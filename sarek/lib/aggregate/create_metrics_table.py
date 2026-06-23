#!/usr/local/miniforge3/bin/python

"""
Collect the data for the "metrics" table and "summary" table in the aggregate report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie, Neil A. Patel"]


from os.path import basename
import pandas as pd
import polars as pl
import narwhals.stable.v2 as nw


def debug_df(name, df, n=10):
    print(f"\n===== {name} (top {n}) =====")
    print(df.head(n))
    print("Shape:", df.shape)


def create_metrics_table(hsmetrics_file: str):
    """
    Return preferred columns from the hsmetrics table.  
    """
    hsmetrics_table = nw.from_native(pl.read_csv(
        hsmetrics_file,
        separator="\t",
    ))

    pl_hsmetrics_table = hsmetrics_table.to_polars()

    pl_hsmetrics_table_reorganized = pl_hsmetrics_table.select(
        pl.col("SAMPLE").alias("sample_collaborator_name"),
        "TOTAL_READS",
        "PF_READS",
        "PCT_PF_READS",
        "PCT_SELECTED_BASES",
        pl.col("MEAN_TARGET_COVERAGE").alias("Mean_full_targets_Coverage"),
        "FOLD_ENRICHMENT",
        pl.col("PCT_EXC_DUPE").alias("PERCENT_DUPLICATION"),
        pl.col("HS_LIBRARY_SIZE").alias("ESTIMATED_LIBRARY_SIZE"), 
        pl.col("PCT_TARGET_BASES_20X").alias("Pct_full_targets_Bases_20x"),
        pl.col("PCT_TARGET_BASES_30X").alias("Pct_full_targets_Bases_30x"),
        pl.col("PCT_TARGET_BASES_50X").alias("Pct_full_targets_Bases_50x"),
        pl.col("PCT_TARGET_BASES_100X").alias("Pct_full_targets_Bases_100x"),
        pl.col("PCT_PF_UQ_READS_ALIGNED"),
        pl.col("ZERO_CVG_TARGETS_PCT").alias("Pct_full_targets_Not_Covered"), 
        pl.col("PCT_USABLE_BASES_ON_BAIT")
    )

    metrics_table = nw.from_native(pl_hsmetrics_table_reorganized)
    debug_df("metrics_table", metrics_table.to_native(), n=10)


    pl_summary_table = pl_hsmetrics_table_reorganized.select(
        pl.col("sample_collaborator_name").alias("Sample ID"),
        pl.col("TOTAL_READS").alias("Total Reads"),
        pl.col("PCT_SELECTED_BASES").alias("% Selected bases"),
        pl.col("Mean_full_targets_Coverage").alias("Mean target Coverage"),
        pl.col("PERCENT_DUPLICATION").alias("% Duplication"),
        pl.col("Pct_full_targets_Bases_20x").alias("% target bases at 20x"),
        pl.col("Pct_full_targets_Bases_30x").alias("% target bases at 30x"),
        pl.col("Pct_full_targets_Bases_50x").alias("% target bases at 50x"),
        pl.col("Pct_full_targets_Bases_100x").alias("% target bases at 100x"),
        pl.col("Pct_full_targets_Not_Covered").alias("% target not covered"),
    )
    summary_table = nw.from_native(pl_summary_table)

    debug_df("summary_table", summary_table.to_native(), n=10)

    return (metrics_table, summary_table)
