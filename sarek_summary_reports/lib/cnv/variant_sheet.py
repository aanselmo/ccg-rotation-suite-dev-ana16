#!/usr/local/miniforge3/bin/python

"""
Gathers the data for the Variant Sheet in the CNV report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os.path import dirname
import asyncio
import polars as pl
import narwhals.stable.v2 as nw
from .qc_sheet import BASENAME_COL


TUMOR_CNS_FILES_COL = "tumor_cns_files"
VCF_FILES_COL = "vcf_files"
SKIP_VCF_LINES = 19
CHROMOSOME_COL = "chromosome"
START_COL = "start"
END_COL = "end"
COPY_NUMBER_COL = "copy_number"


def organize_variant_files(vcf_files: list[str], tumor_cns_files: list[str]) -> nw.DataFrame:
  """
  Organizes all the files by sample/basename into a DataFrame.
  """
  vcf_df = nw.from_native(pl.LazyFrame({VCF_FILES_COL: vcf_files}))
  tumor_cns_df = nw.from_native(pl.LazyFrame({TUMOR_CNS_FILES_COL: tumor_cns_files}))
  vcf_df = vcf_df.with_columns(basename=nw.col(vcf_df.col(VCF_FILES_COL).str.replace(".cnvcall.vcf", "")))
  tumor_cns_df = tumor_cns_df.with_columns(basename=nw.col(
      tumor_cns_df.col(TUMOR_CNS_FILES_COL).str.replace(".tumor_only.call.cns", "")))
  variant_files = vcf_df.join(tumor_cns_df, BASENAME_COL, how="full")
  return variant_files.collect()


def get_one_sample_variants(sample_files: tuple[str, str, str]) -> nw.DataFrame:
  (vcf_file, basename, tumor_cns_file) = sample_files
  vcf_df = nw.from_native(pl.scan_csv(vcf_file, separator="\t", has_header=True, skip_lines=SKIP_VCF_LINES))
  tumor_cns_df = nw.from_native(pl.scan_csv(tumor_cns_file, separator="\t", has_header=True))
  vcf_df.rename({"#CHROM": "chromosome", "POS": "start"})
  filtered_tumor_cns_df = tumor_cns_df.join(vcf_df, on=["chromosome", "start"], how="inner")
  filtered_tumor_cns_df = filtered_tumor_cns_df.with_columns(sample=basename)
  filtered_tumor_cns_df = filtered_tumor_cns_df.select((["sample"] + tumor_cns_df.columns))
  filtered_tumor_cns_df = filtered_tumor_cns_df.remove(nw.col("chromosome") == "X" | nw.col("chromosome") == "Y" )
  return filtered_tumor_cns_df.collect()


async def get_cnv_variants(vcf_files: list[str], tumor_cns_files: list[str]) -> nw.DataFrame:
  """
  A dump of all the official CNV variants output by cnvkit.
  :return nw.DataFrame: A DataFrame with the columns: chromsome, start, end, gene, log2, basename, gene, ci_hi, ci_lo,
      copy_number, depth, probes, weight.
  """
  all_variant_files = organize_variant_files(vcf_files, tumor_cns_files)
  rows_all_variant_files = all_variant_files.rows()
  all_variants = nw.from_native(pl.DataFrame({}))
  for row in rows_all_variant_files:
    one_sample_variants = get_one_sample_variants(row)
    all_variants = nw.concat([all_variants, one_sample_variants], how="vertical")
  all_variant_files.rename({"cn" : COPY_NUMBER_COL})
  return all_variants.collect()
