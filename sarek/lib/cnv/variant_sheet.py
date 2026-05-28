#!/usr/local/miniforge3/bin/python

"""
Gathers the data for the Variant Sheet in the CNV report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os.path import basename as get_basename, dirname
from pathlib import Path
import sys
import polars as pl
import narwhals.stable.v2 as nw
script_dir = Path(__file__).resolve().parent
sarek_script_top_dir = dirname(dirname(script_dir))
sys.path.append(sarek_script_top_dir)
from lib.common import AVAILABLE_CHROMOSOMES
from lib.cnv.qc_sheet import SAMPLE_COL


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
  vcf_df = nw.from_native(pl.DataFrame({VCF_FILES_COL: vcf_files}))
  tumor_cns_df = nw.from_native(pl.DataFrame({TUMOR_CNS_FILES_COL: tumor_cns_files}))
  vcf_df = vcf_df.with_columns(nw.col(VCF_FILES_COL).str.replace(".cnvcall.vcf", "").alias(SAMPLE_COL))
  vcf_df = vcf_df.with_columns(nw.col(SAMPLE_COL).str.split("_vs_").list.get(0).alias(SAMPLE_COL))
  tumor_cns_df = tumor_cns_df.with_columns(
      nw.col(TUMOR_CNS_FILES_COL).str.replace(".tumor_only.call.cns", "").alias(SAMPLE_COL)
  )
  tumor_cns_df = tumor_cns_df.with_columns(
      nw.col(SAMPLE_COL).str.replace(".somatic.call.cns", "").alias(SAMPLE_COL)
  )
  tumor_cns_df = tumor_cns_df.with_columns(
      nw.col(SAMPLE_COL).str.replace(".germline.call.cns", "").alias(SAMPLE_COL)
  )
  vcf_df = vcf_df.to_native().with_columns(
      pl.col(SAMPLE_COL).map_elements(lambda name: get_basename(name)).alias(SAMPLE_COL)
  )
  tumor_cns_df = tumor_cns_df.to_native().with_columns(
      pl.col(SAMPLE_COL).map_elements(lambda name: get_basename(name)).alias(SAMPLE_COL)
  )
  variant_files = vcf_df.join(tumor_cns_df, on=SAMPLE_COL, how="inner")
  variant_files = nw.from_native(variant_files)
  return variant_files


def get_one_sample_variants(basename: str, vcf_file: str, tumor_cns_file: str) -> nw.DataFrame:
  vcf_df = nw.read_csv(
    vcf_file, backend="polars", separator="\t", has_header=True, skip_lines=SKIP_VCF_LINES,
    schema_overrides={"#CHROM": AVAILABLE_CHROMOSOMES}
  )
  tumor_cns_df = nw.read_csv(
      tumor_cns_file, backend="polars", separator="\t", has_header=True,
      schema_overrides={CHROMOSOME_COL: AVAILABLE_CHROMOSOMES}
  )
  vcf_df = vcf_df.rename({"#CHROM": CHROMOSOME_COL, "POS": START_COL})
  vcf_df = vcf_df.with_columns(nw.col(CHROMOSOME_COL).cast(nw.String).alias(CHROMOSOME_COL))
  tumor_cns_df = tumor_cns_df.with_columns(nw.col(CHROMOSOME_COL).cast(nw.String).alias(CHROMOSOME_COL))

  # Ensure the start/position column is an integer in both DataFrames
  vcf_df = vcf_df.with_columns(nw.col(START_COL).cast(nw.Int64))
  tumor_cns_df = tumor_cns_df.with_columns(nw.col(START_COL).cast(nw.Int64))

  filtered_tumor_cns_df = tumor_cns_df.join(vcf_df, on=[CHROMOSOME_COL, START_COL], how="inner")

  filtered_tumor_cns_df = filtered_tumor_cns_df.with_columns(nw.lit(basename).alias(SAMPLE_COL))
  columns_to_select = tumor_cns_df.collect_schema().names()
  columns_to_select.append(SAMPLE_COL)
  filtered_tumor_cns_df = filtered_tumor_cns_df.select(columns_to_select)
  filtered_tumor_cns_df = filtered_tumor_cns_df.filter(
      (nw.col(CHROMOSOME_COL) != "X") & (nw.col(CHROMOSOME_COL) != "Y")
  )
  return filtered_tumor_cns_df


async def get_cnv_variants(vcf_files: list[str], tumor_cns_files: list[str]) -> nw.DataFrame:
  """
  A dump of all the official CNV variants output by cnvkit.
  :return nw.DataFrame: A DataFrame with the columns: chromsome, start, end, gene, log2, basename, gene, ci_hi, ci_lo,
      copy_number, depth, probes, weight.
  """
  all_variant_files = organize_variant_files(vcf_files, tumor_cns_files)
  rows_all_variant_files = all_variant_files.rows()
  all_variants = nw.from_native(pl.DataFrame({}))
  for vcf_file, joint_basename, tumor_cns_file in rows_all_variant_files:
    one_sample_variants = get_one_sample_variants(joint_basename, vcf_file, tumor_cns_file)
    all_variants = nw.concat([all_variants, one_sample_variants], how="diagonal")
  all_variants = all_variants.rename({"cn" : COPY_NUMBER_COL})
  # Reorder columns so that sample name comes first.
  sorted_columns = all_variants.columns
  sorted_columns.remove(SAMPLE_COL)
  sorted_columns = [SAMPLE_COL] + sorted_columns
  all_variants = all_variants.select(sorted_columns)
  return all_variants
