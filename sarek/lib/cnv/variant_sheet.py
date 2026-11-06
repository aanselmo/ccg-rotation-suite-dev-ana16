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
SEGMENT_CNS_FILES_COL = "segment_cns_files"
VCF_FILES_COL = "vcf_files"
SKIP_VCF_LINES = 19
CHROMOSOME_COL = "chromosome"
START_COL = "start"
END_COL = "end"
COPY_NUMBER_COL = "copy_number"
CI_HI_COL = "ci_hi"
CI_LO_COL = "ci_lo"
CI_COLS = [CI_HI_COL, CI_LO_COL]
# Canonical column order of a .tumor_only.call.cns file. Used so every sample comes out with the same columns
# in the same order, whichever cnvkit calling mode (tumor_only / somatic / germline) produced its file.
CALL_CNS_COLUMN_ORDER = [
    CHROMOSOME_COL, START_COL, END_COL, "gene", "log2", CI_HI_COL, CI_LO_COL, "cn", "depth", "probes", "weight"
]


def _attach_segment_cns_files(variant_files: pl.DataFrame, cns_files: list[str]) -> pl.DataFrame:
  """
  Adds a column with each sample's unfiltered segment file (<sample>.cns). That file is the source of ci_lo/ci_hi
  (from cnvkit segmetrics --ci) and is used to back-fill them when the call.cns for a sample does not carry them,
  which is the case for germline.call.cns files.
  """
  segment_cns_files = [
      one_file for one_file in cns_files
      if one_file.endswith(".cns") and "call" not in get_basename(one_file) and "bintest" not in get_basename(one_file)
  ]
  if not segment_cns_files:
    return variant_files.with_columns(pl.lit(None, dtype=pl.String).alias(SEGMENT_CNS_FILES_COL))
  segment_cns_df = pl.DataFrame({SEGMENT_CNS_FILES_COL: segment_cns_files}).with_columns(
      pl.col(SEGMENT_CNS_FILES_COL).map_elements(
          lambda name: get_basename(name).removesuffix(".cns"), return_dtype=pl.String
      ).alias(SAMPLE_COL)
  ).unique(subset=SAMPLE_COL, keep="first")
  return variant_files.join(segment_cns_df, on=SAMPLE_COL, how="left")


def organize_variant_files(vcf_files: list[str], tumor_cns_files: list[str],
                           cns_files: list[str] | None = None) -> nw.DataFrame:
  """
  Organizes all the files by sample/basename into a DataFrame.
  Columns: vcf_files, sample, tumor_cns_files, segment_cns_files (segment_cns_files may be null).
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
  variant_files = _attach_segment_cns_files(variant_files, cns_files or [])
  variant_files = variant_files.select([VCF_FILES_COL, SAMPLE_COL, TUMOR_CNS_FILES_COL, SEGMENT_CNS_FILES_COL])
  variant_files = nw.from_native(variant_files)
  return variant_files


def backfill_confidence_intervals(tumor_cns_df: nw.DataFrame, segment_cns_file: str | None) -> nw.DataFrame:
  """
  cnvkit only carries ci_lo/ci_hi through to the call.cns in some calling modes. tumor_only.call.cns keeps them,
  but germline.call.cns does not (germline calling re-merges/filters segments and drops the segmetrics columns).
  When they are missing, take them from the sample's unfiltered <sample>.cns by matching the exact segment
  (chromosome, start, end). Segments that were merged during calling have no exact counterpart and stay null,
  because a CI for a merged segment cannot be derived from the CIs of its parts.
  """
  missing_ci_cols = [col for col in CI_COLS if col not in tumor_cns_df.columns]
  if not missing_ci_cols:
    return tumor_cns_df
  if segment_cns_file:
    segment_cns_df = nw.read_csv(
        segment_cns_file, backend="polars", separator="\t", has_header=True,
        schema_overrides={CHROMOSOME_COL: AVAILABLE_CHROMOSOMES}
    )
    if all(col in segment_cns_df.columns for col in missing_ci_cols):
      segment_cns_df = segment_cns_df.select([CHROMOSOME_COL, START_COL, END_COL] + missing_ci_cols).with_columns(
          nw.col(CHROMOSOME_COL).cast(nw.String),
          nw.col(START_COL).cast(nw.Int64),
          nw.col(END_COL).cast(nw.Int64),
          *[nw.col(col).cast(nw.Float64) for col in missing_ci_cols]
      ).unique(subset=[CHROMOSOME_COL, START_COL, END_COL], keep="first")
      tumor_cns_df = tumor_cns_df.with_columns(nw.col(END_COL).cast(nw.Int64))
      return tumor_cns_df.join(segment_cns_df, on=[CHROMOSOME_COL, START_COL, END_COL], how="left")
  # No source for the CIs: add explicit null columns so the schema is consistent across samples.
  return tumor_cns_df.with_columns(*[nw.lit(None, dtype=nw.Float64).alias(col) for col in missing_ci_cols])


def get_one_sample_variants(basename: str, vcf_file: str, tumor_cns_file: str,
                            segment_cns_file: str | None = None) -> nw.DataFrame:
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
  tumor_cns_df = backfill_confidence_intervals(tumor_cns_df, segment_cns_file)
  # Put the columns in one canonical order (known columns first, anything unexpected after).
  ordered_columns = [col for col in CALL_CNS_COLUMN_ORDER if col in tumor_cns_df.columns]
  ordered_columns += [col for col in tumor_cns_df.columns if col not in ordered_columns]
  tumor_cns_df = tumor_cns_df.select(ordered_columns)

  filtered_tumor_cns_df = tumor_cns_df.join(vcf_df, on=[CHROMOSOME_COL, START_COL], how="inner")

  filtered_tumor_cns_df = filtered_tumor_cns_df.with_columns(nw.lit(basename).alias(SAMPLE_COL))
  columns_to_select = tumor_cns_df.collect_schema().names()
  columns_to_select.append(SAMPLE_COL)
  filtered_tumor_cns_df = filtered_tumor_cns_df.select(columns_to_select)
  filtered_tumor_cns_df = filtered_tumor_cns_df.filter(
      (nw.col(CHROMOSOME_COL) != "X") & (nw.col(CHROMOSOME_COL) != "Y")
  )
  return filtered_tumor_cns_df


async def get_cnv_variants(vcf_files: list[str], tumor_cns_files: list[str],
                           cns_files: list[str] | None = None) -> nw.DataFrame:
  """
  A dump of all the official CNV variants output by cnvkit.
  :param cns_files list[str]: Unfiltered <sample>.cns segment files. Used to back-fill ci_hi/ci_lo for samples whose
      call.cns lacks them (e.g. germline.call.cns).
  :return nw.DataFrame: A DataFrame with the columns: chromsome, start, end, gene, log2, basename, gene, ci_hi, ci_lo,
      copy_number, depth, probes, weight.
  """
  all_variant_files = organize_variant_files(vcf_files, tumor_cns_files, cns_files)
  rows_all_variant_files = all_variant_files.rows()
  all_variants = nw.from_native(pl.DataFrame({}))
  for vcf_file, joint_basename, tumor_cns_file, segment_cns_file in rows_all_variant_files:
    one_sample_variants = get_one_sample_variants(joint_basename, vcf_file, tumor_cns_file, segment_cns_file)
    all_variants = nw.concat([all_variants, one_sample_variants], how="diagonal")
  all_variants = all_variants.rename({"cn" : COPY_NUMBER_COL})
  # Reorder columns so that sample name comes first.
  sorted_columns = all_variants.columns
  sorted_columns.remove(SAMPLE_COL)
  sorted_columns = [SAMPLE_COL] + sorted_columns
  all_variants = all_variants.select(sorted_columns)
  return all_variants
