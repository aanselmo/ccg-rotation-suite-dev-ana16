#!/usr/local/miniforge3/bin/python

"""
Creates the "Variant Summary Table" Sheet.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os.path import basename
import polars as pl
import narwhals.stable.v2 as nw
from .snv_indel_sheet import (
    RECOMMENDED_CONSEQUENCE_FILTER,
    RECOMMENDED_FILTER_COL,
    EXISTING_VARIATION_COL,
    SAMPLE_COL,
    VARIANT_CLASS_COL
)

DB_SNP_COL = "DBSNP"
COSMIC_COL = "COSMIC"
BOTH_COL = "BOTH"
NOVEL_COL = "NOVEL"
SNV_COL = "SNV"
INDEL_COL = "Indel"
TOTAL_COL = "Total"


def group_by_existing_variation(filtered_variants: nw.DataFrame) -> nw.DataFrame:
  """
  Group by: DBSNP, COSMIC, BOTH, NOVEL
  Basically, if it's in a database of known variants, let it be known.
  """
  DBSNP_SUBSTRING = "rs"
  COSMIC_SUBSTRING = "COSV"
  filtered_variants = filtered_variants.with_columns(
      nw.when(nw.col(EXISTING_VARIATION_COL).str.contains(DBSNP_SUBSTRING))
      .then(True)
      .otherwise(False)
      .alias(DB_SNP_COL)
  )
  filtered_variants = filtered_variants.with_columns(
      nw.when(nw.col(EXISTING_VARIATION_COL).str.contains(COSMIC_SUBSTRING))
      .then(True)
      .otherwise(False)
      .alias(COSMIC_COL)
  )
  is_dbsnp_and_cosmic = nw.col(DB_SNP_COL) & nw.col(COSMIC_COL)
  filtered_variants = filtered_variants.with_columns(
      nw.when(is_dbsnp_and_cosmic).then(True).otherwise(False).alias(BOTH_COL),
      nw.when(is_dbsnp_and_cosmic).then(False).otherwise(True).alias(DB_SNP_COL),
      nw.when(is_dbsnp_and_cosmic).then(False).otherwise(True).alias(COSMIC_COL)
  )
  filtered_variants = filtered_variants.with_columns(
      nw.when(nw.col(DB_SNP_COL) | nw.col(COSMIC_COL) | nw.col(BOTH_COL))
      .then(False)
      .otherwise(True)
      .alias(NOVEL_COL)
  )
  filtered_variants = filtered_variants.select([SAMPLE_COL, DB_SNP_COL, COSMIC_COL, BOTH_COL, NOVEL_COL])
  grouped_variants = filtered_variants.group_by(SAMPLE_COL).agg(
      pl.col(DB_SNP_COL, COSMIC_COL, BOTH_COL, NOVEL_COL).count()
  )
  return grouped_variants


def group_by_class(filtered_variants: nw.DataFrame) -> nw.DataFrame:
  filtered_variants = filtered_variants.select([SAMPLE_COL, VARIANT_CLASS_COL])
  grouped_variants = filtered_variants.group_by([SAMPLE_COL, VARIANT_CLASS_COL]).agg(
      nw.col(VARIANT_CLASS_COL).count()
  )
  return grouped_variants


def get_variant_summary(all_variants: nw.DataFrame) -> nw.DataFrame:
  filtered_variants = all_variants.filter(nw.col(RECOMMENDED_FILTER_COL) == True)
  variants_grouped_by_db = group_by_existing_variation(filtered_variants)
  variants_grouped_by_class = group_by_class(filtered_variants)
  summary_table = variants_grouped_by_db.join(variants_grouped_by_class, on=SAMPLE_COL, how="full")
  return summary_table
