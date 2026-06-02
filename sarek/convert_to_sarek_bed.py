#!/usr/local/miniforge3/bin/python

"""
Fixes a bed file:
1. Sorts the chromosomes numerically.
2. Merged overlapping regions.
3. Annotates the regions in a way that CNVKit expects.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
from os.path import basename, join, dirname
import subprocess
from tempfile import TemporaryDirectory
from datetime import date
import polars as pl
import narwhals.stable.v2 as nw
from lib.common import AVAILABLE_CHROMOSOMES
from lib.cnv.chromosome_sheet import CHROMOSOME_COL, START_COL, END_COL
from lib.cnv.gene_sheet import GENE_COL


APPEND_OPTION = "date.bed"
NEW_FILE_ENDING = ".bed"


def setup_arguments():
  parser = argparse.ArgumentParser(description="Sort the chromosomes in the bed file numerically (1,2,3, etc.).")
  parser.add_argument("-i", "--interval-list", dest="interval_list", type=str, required=True,
                      help="The targeted .interval_list for the bait set.")
  parser.add_argument("-r", "--ref-flat-bed", type=str,
                      default="/mnt/storage/centers/ccgd/reference_files/per_tool/cnvkit/20260119_refFlat_no_chr.bed",
                      dest="ref_flat_bed",
                      help=f"The 20260119_refFlat_no_chr.bed used for CNVKit annotations.")
  parser.add_argument("-o", "--output-file", type=str, default=APPEND_OPTION, dest="output_file_name",
                      help=(f"The name of the output file. By default appends {NEW_FILE_ENDING} the original bed in"
                            + " the dir that the original bed is in."))
  return parser.parse_args()


def get_output_file_name(args: argparse.ArgumentParser) -> str:
  output_file_name = args.output_file_name
  today = date.today()
  today_no_dash = today.isoformat().replace("-", "")
  if output_file_name == APPEND_OPTION:
    output_file_name = basename(args.interval_list)
    output_path = dirname(args.interval_list)
    if ".interval_list" in output_file_name:
      output_file_name = output_file_name.replace(".interval_list", f"_{today_no_dash}{NEW_FILE_ENDING}")
    else:
      output_file_name = output_file_name + f"_{today_no_dash}{NEW_FILE_ENDING}"
    output_file_name = join(output_path, output_file_name)
  return output_file_name


def convert_to_bed(tempdir, interval_list: str):
  output_bed_path = join(tempdir, "unsorted.bed")
  convert_command = [
      "/usr/bin/java", "-jar", "/opt/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar",
      "IntervalListToBed", "-I", interval_list,
      "-O", output_bed_path
  ]
  print(f"Running java command: {" ".join(convert_command)}")
  subprocess.run(convert_command, check=True)
  return output_bed_path


def sort_bed(bed_file: str) -> nw.DataFrame:
  bed_df = pl.read_csv(
      bed_file,
      has_header=False,
      separator="\t",
      new_columns=[CHROMOSOME_COL, START_COL, END_COL, GENE_COL],
      schema_overrides={
          CHROMOSOME_COL: pl.String,
          START_COL: pl.Int64,
          END_COL: pl.Int64,
          GENE_COL: pl.String
      }
  )
  bed_df = bed_df.with_columns(
      pl.col(CHROMOSOME_COL).str.replace("chr", "")
        .cast(AVAILABLE_CHROMOSOMES)
        .alias(CHROMOSOME_COL)
  )
  bed_df = bed_df.sort([CHROMOSOME_COL, START_COL, END_COL])
  return nw.from_native(bed_df)


def merge_bed(bed_df: nw.DataFrame) -> nw.DataFrame:
  bed_df = bed_df.with_columns(nw.col(END_COL).shift(1).over(CHROMOSOME_COL).alias("prev_end"))
  bed_df = bed_df.with_columns((nw.col(START_COL) >= nw.col("prev_end")).fill_null(True).alias("is_new_interval"))
  bed_df = bed_df.with_columns(nw.col("is_new_interval").cum_sum().over(CHROMOSOME_COL).alias("interval_group_number"))
  bed_df = bed_df.drop(GENE_COL, "is_new_interval", "prev_end")
  bed_df = bed_df.group_by([CHROMOSOME_COL, "interval_group_number"]).agg(
      nw.col(START_COL).min().alias(START_COL),
      nw.col(END_COL).max().alias(END_COL),
  )
  bed_df = bed_df.select([CHROMOSOME_COL, START_COL, END_COL]).sort([CHROMOSOME_COL, START_COL])
  return bed_df


def annotate_bed_and_write(tempdir, bed_df: nw.DataFrame, ref_flat_bed: str, output_file_name: str):
  input_bed_path = join(tempdir, "unannotated.bed")
  bed_df.to_native().write_csv(input_bed_path, separator="\t", include_header=False)
  cnvkit_annotate_command = [
      "docker", "run",
      "-v", f"{dirname(output_file_name)}:{dirname(output_file_name)}",
      "-v", f"{dirname(ref_flat_bed)}:{dirname(ref_flat_bed)}",
      "-v", f"{tempdir}:{tempdir}",
      "-w", dirname(output_file_name),
      "quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0",
      "cnvkit.py", "target",
      input_bed_path,
      "--annotate", ref_flat_bed,
      "-o", basename(output_file_name)
  ]
  print(f"Running docker command: {" ".join(cnvkit_annotate_command)}")
  subprocess.run(cnvkit_annotate_command, check=True)


def main():
  args = setup_arguments()
  output_file_name = get_output_file_name(args)
  with TemporaryDirectory() as tempdir:
    bed_file = convert_to_bed(tempdir, args.interval_list)
    new_bed = sort_bed(bed_file)
    new_bed = merge_bed(new_bed)
    annotate_bed_and_write(tempdir, new_bed, args.ref_flat_bed, output_file_name)


if __name__ == "__main__":
  main()
