#!/usr/local/miniforge3/bin/python

"""
Collect the data for the "aggregate_summary" table in the aggregate report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import pandas as pd
import polars as pl
import narwhals.stable.v2 as nw
from os.path import basename
from ..common import debug_df


PCT_TARGET_NOT_COV_COL = "% target not covered"


def get_pct_targets(cumcov_dist_file: str, include_not_covered: bool = True):
  """
  Return % targets at 20x, 30x, 50x, 100x
  Optionally include % Target bases not covered
  """

  all_coverages = nw.from_native(pl.read_csv(
      cumcov_dist_file,
      separator="\t",
  ))

  cumcov_df = all_coverages.select(["Sample", "0"])

  max_col = 0
  half_col = 0

  if "20" not in all_coverages.columns:
    number_cols = list(set(all_coverages.columns) - set(["Sample"]))
    number_cols = [int(one_col) for one_col in number_cols]
    max_col = max(number_cols)
    half_col = max_col // 2

    cumcov_df = nw.concat(
        [cumcov_df, all_coverages.select([str(half_col), str(max_col)])],
        how="horizontal"
    )

  for column in ["20", "30", "50", "100"]:
    if column in all_coverages.columns:
      cumcov_df = nw.concat([cumcov_df, all_coverages.select(column)], how="horizontal")
    elif column == "20":
      zero_col = nw.from_native(
          pl.DataFrame({column: all_coverages.shape[0] * [f"({column}, 0)"]})
      )
      cumcov_df = nw.concat([cumcov_df, zero_col], how="horizontal")

  cumcov_df = cumcov_df.filter(nw.col("Sample").str.ends_with(".recal"))

  # --- Start output with Sample column ---
  summarized_pct_targets = cumcov_df.select(
      nw.col("Sample").str.replace(".recal", "")
  )

  # --- Coverage columns ---
  non_zero_coverage_cols = list(set(cumcov_df.columns) - set(["Sample", "0"]))
  non_zero_coverage_cols = sorted([int(c) for c in non_zero_coverage_cols])
  non_zero_coverage_cols = [str(c) for c in non_zero_coverage_cols]

  for cov_column in non_zero_coverage_cols:
    one_pct_target = cumcov_df.select(
        nw.col(cov_column).str.split(", ").list.get(1)
          .str.strip_chars(")")
          .cast(nw.Float64)
          .fill_null(0.0)
        # .cast(nw.Int16)
          .round(2)
          .alias(f"% target bases at {cov_column}x"),
    )

    summarized_pct_targets = nw.concat(
        [summarized_pct_targets, one_pct_target],
        how="horizontal"
    )

  # --- OPTIONAL: Not covered column ---
  if include_not_covered:
    not_covered_col = (
        nw.lit(100)
        - nw.col("0").str.split(", ").list.get(1)
          .str.strip_chars(")")
          .cast(nw.Float64)
        # .cast(nw.Int16)
          .round(2)
    ).alias(PCT_TARGET_NOT_COV_COL)

    summarized_pct_targets = nw.concat(
        [summarized_pct_targets, cumcov_df.select(not_covered_col)],
        how="horizontal"
    )

  # --- Rename columns to match old riker reports ---
  summarized_pct_targets = summarized_pct_targets.rename({
      "Sample": "Sample ID"
  })

  debug_df("Summarized_pct_targets", summarized_pct_targets.to_native(), n=10)

  return summarized_pct_targets


def get_mtc_and_pct_selected(mos_files) -> nw.DataFrame:
  all_mtc_and_pct_selected = nw.from_native(pl.DataFrame({}))
  for one_file in mos_files:
    one_file_metrics = nw.from_native(pl.read_csv(one_file, separator="\t"))
    mtc = one_file_metrics.filter(nw.col("chrom") == nw.lit("total_region")).get_column("mean")[0]
    total_bases = one_file_metrics.filter(nw.col("chrom") == nw.lit("total")).get_column("bases")[0]
    total_bases_on_target = one_file_metrics.filter(nw.col("chrom") == nw.lit("total_region")).get_column("bases")[0]
    pct_selected_bases = total_bases_on_target / total_bases
    sample = basename(one_file).replace(".recal.mosdepth.summary.txt", "")
    one_row = {
        "Sample ID": sample,
        "% Selected bases": pct_selected_bases * 100.0,
        "Mean target Coverage": mtc,
    }
    all_mtc_and_pct_selected = nw.concat([all_mtc_and_pct_selected, nw.from_dict(one_row, backend="polars")],
                                         how="vertical")
  debug_df("All_MTC_and_Pct", all_mtc_and_pct_selected, n=10)
  return all_mtc_and_pct_selected


def get_total_reads(sam_files: list[str]) -> nw.DataFrame:
  all_total_reads = nw.from_native(pl.DataFrame({}))
  for file_path in sam_files:
    sample_name = basename(file_path).replace(".md.cram.stats", "")
    total_reads = 0
    with open(file_path, "r") as one_file:
      for line in one_file:
        if "raw total sequences" in line:
          total_reads = line.split("\t")[2]  # 3rd column
          break
    one_row = {
        "Sample ID": sample_name,
        "Total Reads": total_reads,
    }
    all_total_reads = nw.concat([all_total_reads, nw.from_dict(one_row, backend="polars")], how="vertical")
  debug_df("All total reads", all_total_reads, n=10)
  return all_total_reads


def get_pct_duplication(duplicate_marking_files: list[str]) -> nw.DataFrame:
  all_pct_duplication = nw.from_native(pl.DataFrame({}))
  for file_path in duplicate_marking_files:
    sample_name = basename(file_path).replace(".md.cram.metrics", "")
    pct_dup = 0.0
    start_row = 0
    with open(file_path, "r") as one_file:
      for line_number, line in enumerate(one_file):
        if "PERCENT_DUPLICATION" in line:
          start_row = line_number
          break
    one_file_metrics = nw.from_native(pl.read_csv(file_path, skip_lines=start_row, n_rows=1, separator="\t"))
    pct_dup = one_file_metrics.get_column("PERCENT_DUPLICATION")[0]
    one_row = {
        "Sample ID": sample_name,
        "% Duplication": pct_dup * 100.0,
    }
    all_pct_duplication = nw.concat([all_pct_duplication, nw.from_dict(one_row, backend="polars")], how="vertical")
  return all_pct_duplication


def create_summary_table(
    cumcov_dist_file: str, duplicate_marking_files: list[str], mos_files: list[str], sam_files: list[str]
) -> pd.DataFrame:
  summarized_pct_targets = get_pct_targets(cumcov_dist_file)
  total_reads = get_total_reads(sam_files)
  mtc_and_pct_selected = get_mtc_and_pct_selected(mos_files)
  pct_duplication = get_pct_duplication(duplicate_marking_files)
  summary_table = total_reads.join(mtc_and_pct_selected, on="Sample ID", how="full", suffix="_mos")
  summary_table = summary_table.join(pct_duplication, on="Sample ID", how="full", suffix="_dup")
  summary_table = summary_table.join(summarized_pct_targets, on="Sample ID", how="full", suffix="_target")
  summary_table = summary_table.with_columns(
      nw.coalesce("Sample ID", "Sample ID_mos", "Sample ID_target", "Sample ID_dup", nw.lit("")).alias("Sample ID")
  )
  summary_table = summary_table.drop(["Sample ID_mos", "Sample ID_target", "Sample ID_dup"])
  return summary_table.to_pandas()
