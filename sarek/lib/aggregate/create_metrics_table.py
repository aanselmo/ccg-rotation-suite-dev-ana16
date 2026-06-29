#!/usr/local/miniforge3/bin/python

"""
Organize the data for the "metrics" table and "summary" table in the aggregate report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie, Neil A. Patel"]


import polars as pl
import narwhals.stable.v2 as nw
from ..common import debug_df
from ..setup.google import *


def create_metrics_table(data_lims, data_picard, data_duplicates):
  """
  Comine miniLIMS, picard, and sarek data into format for reporting.
  """
  pl_data_lims = data_lims.to_polars().select(
      "sample",
      pl.col("library").alias("Library"),
      pl.col("barcode_set").alias("Barcode Set"),
      pl.col("barcode_member").alias("Index"),
      pl.col("first_barcode").alias("First barcode"),
      pl.col("second_barcode").alias("Second barcode"),
      pl.col("subject_id").alias("Subject ID"),
      pl.col("capture_name").alias("Capture name")
  )
  pl_data_picard = data_picard.to_polars().with_columns(
      (pl.col("PERCENT_DUPLICATION").alias("% Bases Duplicated"))
  )
  pl_data_duplicates = data_duplicates.to_polars().with_columns(
      (pl.col("% Duplication") / 100).alias("% Reads Duplicated")  # Sarek "pct duplicated"
  ).drop(pl.col("% Duplication"))

  pl_metrics_table = (
      pl_data_lims
      .join(pl_data_picard, how="full", left_on="sample", right_on="sample_collaborator_name")
      .join(pl_data_duplicates, how="full", left_on="sample", right_on="Sample ID")
  ).select(
      ["Sample ID"]
      + pl_data_lims.drop("sample").columns
      + pl_data_picard.drop("sample_collaborator_name", "PERCENT_DUPLICATION", "% Bases Duplicated").columns
      + ["% Bases Duplicated"]
      + pl_data_duplicates.drop("Sample ID").columns
  )
  metrics_table = nw.from_native(pl_metrics_table)
  debug_df("metrics_table", metrics_table.to_native(), n=10)

  pl_summary_table = pl_metrics_table.select(
      "Sample ID", "Library", "Index", "Capture name",
      pl.col("TOTAL_READS").alias("Total Reads"),
      pl.col("PCT_SELECTED_BASES").alias("% Selected bases"),
      pl.col("Mean_full_targets_Coverage").alias("Mean target Coverage"),
      pl.col("% Bases Duplicated"),  # from Picard
      (pl.col("% Reads Duplicated").alias("% Duplication")),  # from Sarek
      pl.col("Pct_full_targets_Bases_20x").alias("% target bases at 20x"),
      pl.col("Pct_full_targets_Bases_30x").alias("% target bases at 30x"),
      pl.col("Pct_full_targets_Bases_50x").alias("% target bases at 50x"),
      pl.col("Pct_full_targets_Bases_100x").alias("% target bases at 100x"),
      pl.col("Pct_full_targets_Not_Covered").alias("% target not covered"),
  )
  summary_table = nw.from_native(pl_summary_table)

  debug_df("summary_table", summary_table.to_native(), n=10)

  return (metrics_table, summary_table)
