#!/usr/local/miniforge3/bin/python

"""
Gathers the data for the ArmLevel Sheet in the CNV report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel", "Anthony Anselmo"]

import os
from os.path import join, dirname, abspath, basename as get_basename
import asyncio
import subprocess
import io
import re
import numpy as np
import polars as pl
import narwhals.stable.v2 as nw
import matplotlib.pyplot as plt

ARM_IMAGE = "aanselmo/cnv_arm:0.1"
ARM_SCRIPT = "/opt/arm_level_calls_from_cnvkit.R"  # script within docker container

CNR_FILES_COL = "cnr_files"
CNS_FILES_COL = "cns_files"
CALL_CNS_FILES_COL = "call_cns_files"
SAMPLE_COL = "sample"


def _sample_table(files: list[str], file_col: str, ext: str) -> pl.DataFrame:
  df = pl.DataFrame({file_col: files}, schema={file_col: pl.String}).with_columns(
      pl.col(file_col).str.replace(r"^.*/", "").str.replace(rf"{re.escape(ext)}$", "").alias(SAMPLE_COL)
  )
  dupes = df.filter(pl.col(SAMPLE_COL).is_duplicated())
  if dupes.height:
    raise ValueError(f"Duplicate sample names in {file_col}:\n{dupes}")
  return df

def organize_files(cnr_files: list[str], cns_files: list[str], call_cns_files: list[str]) -> nw.DataFrame:
  """Organizes all the files by sample/basename into a DataFrame."""
  cnr_df = _sample_table(cnr_files, CNR_FILES_COL, ".cnr")
  cns_df = _sample_table([f for f in cns_files if not f.endswith(".call.cns")], CNS_FILES_COL, ".cns")
  call_df = _sample_table(call_cns_files, CALL_CNS_FILES_COL, ".call.cns")
  all_files = (
      cnr_df.join(cns_df, on=SAMPLE_COL, how="full", coalesce=True)
            .join(call_df, on=SAMPLE_COL, how="full", coalesce=True)
  )
  return nw.from_native(all_files)


async def run_armlevel_one_sample(row: dict, abs_bucket_dir: str, out_dir: str,
                                  sem: asyncio.Semaphore) -> pl.DataFrame:
  sample = row[SAMPLE_COL]
  out_path = join(out_dir, f"{sample}.arm_calls.tsv")
  docker_command = [
      "docker", "run", "--rm",
      #"--user", f"{os.getuid()}:{os.getgid()}",
      "-v", f"{abs_bucket_dir}:{abs_bucket_dir}", "-w", abs_bucket_dir,
      ARM_IMAGE,
      "Rscript", ARM_SCRIPT,
      "--cnr", row[CNR_FILES_COL],
      "--cns", row[CNS_FILES_COL],
      "--call", row[CALL_CNS_FILES_COL],
      "--sample", sample,
      "--cytoband", "/opt/cytoBand.txt.gz",
      "--fraction-by", "bins",
      "--out", out_path,
  ]
  async with sem:
    proc = await asyncio.create_subprocess_exec(
        *docker_command, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
    )
    stdout, stderr = await proc.communicate()
  if proc.returncode != 0:
    raise RuntimeError(f"{sample}: arm-level caller failed (exit {proc.returncode}):\n{stderr.decode()}")
  print(stdout.decode().strip())  # the R script's one-line summary
  return pl.read_csv(out_path, separator="\t", null_values=["NA"])


async def run_armlevel_caller(files: nw.DataFrame, bucket_dir: str, max_parallel: int = 4) -> nw.DataFrame:
  abs_bucket_dir = abspath(bucket_dir)
  out_dir = join(abs_bucket_dir, "arm_calls")
  os.makedirs(out_dir, exist_ok=True)  # write.table won't create the directory
  sem = asyncio.Semaphore(max_parallel)
  results = await asyncio.gather(
      *(run_armlevel_one_sample(row, abs_bucket_dir, out_dir, sem) for row in files.rows(named=True))
  )
  return nw.from_native(pl.concat(results, how="vertical_relaxed"))


async def get_armlevel_calls(cnr_files, cns_files, call_cns_files, bucket_dir="") -> nw.DataFrame:
  if not bucket_dir:
    bucket_dir = dirname(dirname(dirname(cnr_files[0])))
  all_files = organize_files(cnr_files, cns_files, call_cns_files)
  complete = all_files.select(SAMPLE_COL, CNR_FILES_COL, CNS_FILES_COL, CALL_CNS_FILES_COL).drop_nulls()
  armlevel_calls = await run_armlevel_caller(complete, bucket_dir)
  return armlevel_calls.drop("purity")
