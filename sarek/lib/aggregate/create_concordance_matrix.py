#!/usr/local/miniforge3/bin/python

"""
Collect the data for the "FP_concordance_matrix" table in the aggregate report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie"]


import polars as pl
import narwhals.stable.v2 as nw
from ..common import debug_df


def create_concordance_matrix(fingerprint_table):
  """
  Transform the somalier paired data output into a wide matrix.  
  """

  # to make a complete matrix, combine (1) original paired data (2) inverted data (swap sample_a and sample_b
  # and (3) self-self pairs with relatedness value set to 1.0.

  pl_fingerprint_original = (
      fingerprint_table.to_native()
      .select("Sample1", "Sample2", "Relatedness")
      .with_columns(pl.col("Relatedness").cast(pl.Float64, strict=False).fill_nan(None))
  )
  pl_fingerprint_inverted = (pl_fingerprint_original.rename({"Sample1": "Sample2", "Sample2": "Sample1"})
                             .select("Sample1", "Sample2", "Relatedness"))
  fingerprint_samples_unique = (list(set(fingerprint_table["Sample1"].to_list())
                                     .union(set(fingerprint_table["Sample2"].to_list()))))
  pl_fingerprint_self = pl.DataFrame({
      "Sample1": fingerprint_samples_unique,
      "Sample2": fingerprint_samples_unique,
      "Relatedness": 1.0
  })

  pl_fingerprint_tall = pl.concat([pl_fingerprint_original, pl_fingerprint_inverted, pl_fingerprint_self],
                                  how="vertical").sort(by=["Sample1", "Sample2"])

  pl_fingerprint_matrix = (pl_fingerprint_tall.pivot(index="Sample1", on="Sample2", values="Relatedness")
                           .rename({"Sample1": "Somalier Relatedness"}))

  concordance_matrix = nw.from_native(pl_fingerprint_matrix)
  debug_df("concordance_matrix", concordance_matrix.to_native(), n=10)

  return concordance_matrix
