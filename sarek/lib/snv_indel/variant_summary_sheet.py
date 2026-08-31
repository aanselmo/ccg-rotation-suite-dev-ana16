#!/usr/local/miniforge3/bin/python

"""
Creates the "Variant Summary Table" Sheet.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import polars as pl
import narwhals.stable.v2 as nw
from .vcf2df import (
    BEST_EFFECT_CONSEQUENCE_COL,
    EXISTING_VARIATION_COL,
    COSMIC_v104_COL,
    CSQ_VARIANT_CLASS_COL
)
from .snv_indel_sheet import (
    MUTECT2_FILTER_COL,
    RECOMMENDED_CONSEQUENCE_FILTER,
    REVIEW_REQUIRED_COL,
    RECOMMENDED_FILTER_COL,
    SAMPLE_COL
)

DB_SNP_COL = "DBSNP"
COSMIC_COL = "COSMIC_v104"
BOTH_COL = "BOTH"
NOVEL_COL = "NOVEL"
SNV_COL = "SNV"
INDEL_COL = "Indel"
TOTAL_COL = "Total"
REPEAT_COL_SUFFIX = "_repeat"


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
      nw.when(nw.col(COSMIC_v104_COL).str.contains(COSMIC_SUBSTRING))
      .then(True)
      .otherwise(False)
      .alias(COSMIC_COL)
  )
  is_dbsnp_and_cosmic = nw.col(DB_SNP_COL) & nw.col(COSMIC_COL)
  filtered_variants = filtered_variants.with_columns(
      nw.when(is_dbsnp_and_cosmic).then(True).otherwise(False).alias(BOTH_COL),
      nw.when(is_dbsnp_and_cosmic).then(False).otherwise(nw.col(DB_SNP_COL)).alias(DB_SNP_COL),
      nw.when(is_dbsnp_and_cosmic).then(False).otherwise(nw.col(COSMIC_COL)).alias(COSMIC_COL)
  )
  filtered_variants = filtered_variants.with_columns(
      nw.when(nw.col(DB_SNP_COL) | nw.col(COSMIC_COL) | nw.col(BOTH_COL))
      .then(False)
      .otherwise(True)
      .alias(NOVEL_COL)
  )
  filtered_variants = filtered_variants.select([SAMPLE_COL, DB_SNP_COL, COSMIC_COL, BOTH_COL, NOVEL_COL])
  grouped_variants = filtered_variants.group_by(SAMPLE_COL).agg(
      nw.col(DB_SNP_COL).sum(),
      nw.col(COSMIC_COL).sum(),
      nw.col(BOTH_COL).sum(),
      nw.col(NOVEL_COL).sum(),
  )
  # Blank spaces instead of 0 is easier to read in excel. Otherwise, normally, I don't like None/null values.
  grouped_variants = grouped_variants.with_columns(
      nw.when(nw.col(DB_SNP_COL) == nw.lit(0)).then(None).otherwise(nw.col(DB_SNP_COL)).alias(DB_SNP_COL),
      nw.when(nw.col(COSMIC_COL) == nw.lit(0)).then(None).otherwise(nw.col(COSMIC_COL)).alias(COSMIC_COL),
      nw.when(nw.col(BOTH_COL) == nw.lit(0)).then(None).otherwise(nw.col(BOTH_COL)).alias(BOTH_COL),
      nw.when(nw.col(NOVEL_COL) == nw.lit(0)).then(None).otherwise(nw.col(NOVEL_COL)).alias(NOVEL_COL),
  )
  return grouped_variants


def group_by_class(filtered_variants: nw.DataFrame) -> nw.DataFrame:
  CLASS_COUNT_COL = "CLASS_COUNT"
  filtered_variants = filtered_variants.select([SAMPLE_COL, CSQ_VARIANT_CLASS_COL])
  grouped_variants = filtered_variants.group_by([SAMPLE_COL, CSQ_VARIANT_CLASS_COL]).agg(
      nw.col(CSQ_VARIANT_CLASS_COL).count().alias(CLASS_COUNT_COL)
  )
  all_samples = grouped_variants.get_column(SAMPLE_COL).unique().to_list()
  grouped_partially_transposed = nw.from_native(pl.DataFrame({}))
  for one_sample in all_samples:
    one_group = grouped_variants.filter(nw.col(SAMPLE_COL) == one_sample).drop(SAMPLE_COL)
    one_group = nw.from_native(
        one_group.to_native()
        .transpose(column_names=CSQ_VARIANT_CLASS_COL)
        .with_columns(pl.lit(one_sample).alias(SAMPLE_COL))
    )
    grouped_partially_transposed = nw.concat([grouped_partially_transposed, one_group], how="diagonal")
  return grouped_partially_transposed


def get_variant_summary(all_variants: nw.DataFrame) -> nw.DataFrame:
  filtered_variants = all_variants.filter(nw.col(RECOMMENDED_FILTER_COL) == True)
  variants_grouped_by_db = group_by_existing_variation(filtered_variants)
  variants_grouped_by_class = group_by_class(filtered_variants)
  summary_table = variants_grouped_by_db.join(
      variants_grouped_by_class, on=SAMPLE_COL, how="full", suffix=REPEAT_COL_SUFFIX)
  # Sort columns
  left_columns = [SAMPLE_COL, DB_SNP_COL, COSMIC_COL, BOTH_COL, NOVEL_COL, f"{SAMPLE_COL}{REPEAT_COL_SUFFIX}"]
  right_columns = set(summary_table.columns).difference() - set(left_columns)
  reordered_list = left_columns + list(right_columns)
  summary_table = summary_table.select(reordered_list)
  # Sort rows
  summary_table = summary_table.sort(by=SAMPLE_COL)
  summary_table = summary_table.with_columns(
      nw.lit(
          f"{MUTECT2_FILTER_COL}='PASS';{REVIEW_REQUIRED_COL}='TRUE';{BEST_EFFECT_CONSEQUENCE_COL}={str(RECOMMENDED_CONSEQUENCE_FILTER)}"
      ).alias("ccg_recommended_filter")
  )
  return summary_table
