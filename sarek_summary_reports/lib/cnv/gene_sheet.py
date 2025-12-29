#!/usr/local/miniforge3/bin/python

"""
Gathers the data for the Summary Sheet by Gene in the CNV report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os.path import dirname, basename as get_basename
import asyncio
import polars as pl
import narwhals.stable.v2 as nw
from .qc_sheet import BASENAME_COL
from .variant_sheet import COPY_NUMBER_COL, CHROMOSOME_COL

GENE_COL = "gene"
LOG2_COL = "log2"
DIPLOID = 2.0
GAIN_COL = "No. of Samples with Gains"
LOSS_COL = "No. of Samples with Losses"
GAIN_AND_LOSS_COL = "No. of Samples with Gains & Losses"


def get_one_sample_variants(one_gene_file: str) -> nw.DataFrame:
  one_sample_gene_cnvs = nw.from_native(pl.read_csv(one_gene_file, separator="\t", has_header=True, row_index_offset=0))
  # The last row has the mean target average
  one_sample_gene_cnvs = one_sample_gene_cnvs.with_columns(nw.col(CHROMOSOME_COL).cast(nw.String).alias(CHROMOSOME_COL))
  one_sample_gene_cnvs = one_sample_gene_cnvs.filter((nw.col(CHROMOSOME_COL) != "X") & (nw.col(CHROMOSOME_COL) != "Y"))
  one_sample_gene_cnvs = one_sample_gene_cnvs.with_columns(nw.col(CHROMOSOME_COL).cast(nw.Int8).alias(CHROMOSOME_COL))
  one_sample_gene_cnvs = one_sample_gene_cnvs.with_columns(
      nw.lit(one_gene_file.replace(".genemetrics.tsv", "")).alias(BASENAME_COL)
  )
  one_sample_gene_cnvs = one_sample_gene_cnvs.to_native().with_columns(
      pl.col(BASENAME_COL).map_elements(lambda filename: get_basename(filename)).alias(BASENAME_COL)
  )
  one_sample_gene_cnvs = nw.from_native(one_sample_gene_cnvs)
  return one_sample_gene_cnvs


def calculate_copy_number(original_df: nw.DataFrame, log2col=LOG2_COL) -> nw.DataFrame:
  """
  Take a DataFrame with log2 ratios and calculate the copy numbers.
  :return nw.DataFrame: return an uncollected DataFrame with a new column of copy_number.
  """
  original_df = original_df.with_columns((DIPLOID * (2.0 ** nw.col(log2col))).round(decimals=0).alias(COPY_NUMBER_COL))
  return original_df


def sum_cnv_numbers_by_gene(all_gene_cnvs: nw.DataFrame) -> nw.DataFrame:
  # Log 2 ratio means 2.0 to the power of the log2ratio
  all_gene_cnvs = calculate_copy_number(all_gene_cnvs, log2col=LOG2_COL)
  all_gene_summary = all_gene_cnvs.group_by(GENE_COL, BASENAME_COL).agg(
      nw.when(nw.col(COPY_NUMBER_COL) > DIPLOID).then(1).otherwise(0).sum().alias(GAIN_COL),
      nw.when(nw.col(COPY_NUMBER_COL) < DIPLOID).then(1).otherwise(0).sum().alias(LOSS_COL)
  )
  is_gain_and_loss = (nw.col(GAIN_COL) == 1) & (nw.col(LOSS_COL) == 1)
  all_gene_summary = all_gene_summary.with_columns(
      nw.when(is_gain_and_loss).then(1).otherwise(0).alias(GAIN_AND_LOSS_COL),
      nw.when(is_gain_and_loss).then(0).otherwise(nw.col(GAIN_COL)).alias(GAIN_COL),
      nw.when(is_gain_and_loss).then(0).otherwise(nw.col(LOSS_COL)).alias(LOSS_COL),
  )
  all_gene_summary = all_gene_summary.group_by(GENE_COL).agg(
      nw.col(GAIN_COL).sum(),
      nw.col(LOSS_COL).sum(),
      nw.col(GAIN_AND_LOSS_COL).sum()
  )
  return all_gene_summary


async def summarize_gene_level(gene_files: list[str]) -> nw.DataFrame:
  """
  Summarizes all the calls (gains/losses) by gene.
  :param gene_files list[str]: A list of file paths ending in .genemetrics.tsv
  :return nw.DataFrame: A data frame with the columns: gene, No. of Samples with Gains, No. of Samples with Losses,
      and No. of Samples with Gains & Losses.
  """
  all_gene_cnvs = nw.from_native(pl.DataFrame({}))
  for one_gene_file in gene_files:
    one_sample_variants = get_one_sample_variants(one_gene_file)
    all_gene_cnvs = nw.concat([all_gene_cnvs, one_sample_variants], how="vertical")
  all_gene_summary = sum_cnv_numbers_by_gene(all_gene_cnvs)
  return all_gene_summary
