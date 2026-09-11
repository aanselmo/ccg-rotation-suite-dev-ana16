#!/usr/local/miniforge3/bin/python

"""
Creates the 'SNV/Indel Variants' Sheet.

Most of the filtering and removal of columns happens here and not in vcf2df.py.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import sys
from os.path import basename, dirname
from pathlib import Path
import polars as pl
import narwhals.stable.v2 as nw
script_dir = Path(__file__).resolve().parent
sarek_script_top_dir = dirname(dirname(script_dir))
sys.path.append(sarek_script_top_dir)
from lib.snv_indel.vcf2df import (
    convert_vcf_to_data_frame, BEST_EFFECT_CONSEQUENCE_COL, CSQ_BIOTYPE_COL, BEST_EFFECT_COLS_TO_USE,
    BEST_EFFECT_PREFIX, CANONICAL_PREFIX, CANONICAL_COLS_TO_USE, CSQ_VARIANT_TYPE_COLS, CSQ_POPULATION_COLS,
    COSMIC_v104_COL, CHROM_COL, POSITION_COL
)
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context


VEP_FILE_ENDING = ".mutect2.filtered_VEP.ann.vcf.gz"
SAMPLE_COL = "tumor_sample_name"
REF_COL = "REF"
ALT_COL = "ALT"
RECOMMENDED_FILTER_COL = "ccg_recommended_filter"
REVIEW_REQUIRED_COL = "Review_Required"
MUTECT2_FILTER_COL = "FILTER"
RECOMMENDED_CONSEQUENCE_FILTER = [
    "missense_variant", "frameshift_variant", "inframe_deletion", "inframe_insertion",
    "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant",
    "stop_lost", "start_lost", "stop_gained", "splice_region_variant",
    # These two don't typically appear in SNV/Indels but I'm keeping it here.
    "transcript_ablation", "transcript_amplification"
]


def setup_arguments():
  parser = argparse.ArgumentParser(description=(
      "This snv_indel_sheet.py script was designed to parse 1 file. Otherwise, because the snv_indel_report"
      " is parallelized/multi-process, the stack trace of the error can be obfuscated."
  ))
  parser.add_argument("-o", "--output-file", type=str, default="output.csv",
                      dest="output_file_name", help="The name of the output csv file.")
  parser.add_argument("-v", "--vcf", dest="vcf", type=str, required=True,
                      help=f"A vcf.gz file with the actual variant calls in them.")
  return parser.parse_args()


def add_review_required(all_sample_variants: nw.DataFrame) -> nw.DataFrame:
  sample_variants: pl.DataFrame = all_sample_variants.to_native()

  # gnomADe_AF may be read as String due to "." being used for missing values in VCF
  sample_variants = sample_variants.with_columns(
      pl.col("gnomADe_AF")
      .replace([".", ""], None)
      .cast(pl.Float64)
      .alias("gnomADe_AF")
  )

  sample_variants = sample_variants.with_columns(
      pl.when(
          (pl.col(CSQ_BIOTYPE_COL) == "protein_coding") & (
              (pl.col("gnomADe_AF") <= 0.01)
              | ((pl.col("gnomADe_AF") <= 0.10) & (pl.col(COSMIC_v104_COL).str.contains("COSV")))
              | (pl.col("gnomADe_AF").is_null())
          )
      )
      .then(True)
      .otherwise(False)
      .alias(REVIEW_REQUIRED_COL)
  )
  return nw.from_native(sample_variants)


def add_recommended_consequence_filter(all_sample_variants: nw.DataFrame) -> nw.DataFrame:
  sample_variants: pl.DataFrame = all_sample_variants.to_native()
  sample_variants = sample_variants.with_columns(
      pl.when((pl.col(BEST_EFFECT_CONSEQUENCE_COL)
                 .str.split("&")
                 .list.set_intersection(pl.lit(RECOMMENDED_CONSEQUENCE_FILTER))
                 .list.len() > 0
               ) & (pl.col(MUTECT2_FILTER_COL) == "PASS") & (pl.col(REVIEW_REQUIRED_COL) == True))
      .then(True)
      .otherwise(False)
      .alias(RECOMMENDED_FILTER_COL)
  )
  return nw.from_native(sample_variants)


def process_one_file(one_vcf_file_path: str) -> nw.DataFrame:
  one_sample_variants = convert_vcf_to_data_frame(one_vcf_file_path)
  sample_name = basename(one_vcf_file_path).replace(VEP_FILE_ENDING, "")
  one_sample_variants = one_sample_variants.with_columns(nw.lit(sample_name).alias(SAMPLE_COL))
  one_sample_variants = add_review_required(one_sample_variants)
  one_sample_variants = add_recommended_consequence_filter(one_sample_variants)

  # KEEP ONLY PASSING MUTECT2 VARIANTS
  one_sample_variants = one_sample_variants.filter(
      nw.col(MUTECT2_FILTER_COL) == "PASS"
  )

  # Select only these columns for now. There are a lot more columns. I made an arbitrary call to get it as similar to
  # the riker report as possible. Also, there a lot of gnomAD columns so I tried to summarize them.
  best_effect_cols = [BEST_EFFECT_PREFIX + col for col in BEST_EFFECT_COLS_TO_USE]
  canonical_cols = [CANONICAL_PREFIX + col for col in CANONICAL_COLS_TO_USE]
  one_sample_variants = one_sample_variants.select(
      [SAMPLE_COL, CHROM_COL, POSITION_COL, REF_COL, ALT_COL, MUTECT2_FILTER_COL, REVIEW_REQUIRED_COL,
       RECOMMENDED_FILTER_COL] + CSQ_VARIANT_TYPE_COLS + canonical_cols + best_effect_cols + ["AF", "DP"]
      + CSQ_POPULATION_COLS
  )
  return one_sample_variants


def get_variants(vcf_files: list[str], number_processes=1) -> nw.DataFrame:
  all_sample_variants = nw.from_native(pl.DataFrame({}))
  with ProcessPoolExecutor(max_workers=number_processes, mp_context=get_context("spawn")) as executor:
    list_of_variant_dfs = list(executor.map(process_one_file, vcf_files))
  all_sample_variants = nw.concat(list_of_variant_dfs, how="diagonal")
  all_sample_variants.sort(by=[SAMPLE_COL, CHROM_COL, POSITION_COL])
  return all_sample_variants


def main():
  """
  Main function here was written as a way to quickly test parsing of 1 file. Otherwise, because the snv_indel_report
  is parallelized / multi - process, the stack trace of the error can be obfuscated.
  """
  args = setup_arguments()
  df = process_one_file(args.vcf)
  df.write_csv(args.output_file_name)


if __name__ == "__main__":
  main()
