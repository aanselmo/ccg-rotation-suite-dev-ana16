#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the vcf2df.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
from os.path import join
from pathlib import Path
import polars as pl
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.vcf2df as vcf2df
from sarek.lib.common import AVAILABLE_CHROMOSOMES


EXAMPLE_VCF_PATH = join(Path(__file__).parent, "data", "A008_X108_58797-1-C1.mutect2.filtered_VEP.ann.vcf.gz")


def test_convert_vcf_to_data_frame():
  sample_a008_df = vcf2df.convert_vcf_to_data_frame(EXAMPLE_VCF_PATH, should_filter_mutect2_failed=False)
  # Assert columns that should be further split are indeed split.
  assert "INFO" not in sample_a008_df.columns
  assert "CSQ" not in sample_a008_df.columns
  # Assert 1 variant has the correct alternate allele
  assert sample_a008_df.filter(
      nw.col(vcf2df.POSITION_COL) == 155206037
  ).get_column("ALT")[0] == "A"


def test_split_format_column():
  (header_line, info_col_names, csq_col_names) = vcf2df.parse_vcf_header(EXAMPLE_VCF_PATH)
  sample_a008_df = pl.read_csv(
      EXAMPLE_VCF_PATH,
      skip_lines=header_line,
      separator="\t",
      schema_overrides={vcf2df.CHROM_COL: AVAILABLE_CHROMOSOMES},
  )
  # Data Col of Fromat Header Starts off as the sample name so its not constant.
  sample_a008_df = vcf2df.split_format_column(sample_a008_df)
  assert sample_a008_df.filter(
      pl.col(vcf2df.POSITION_COL) == 135526367
  ).get_column("DP")[0] == "355"


def test_split_csq_column():
  (header_line, info_col_names, csq_col_names) = vcf2df.parse_vcf_header(EXAMPLE_VCF_PATH)
  sample_a008_df = pl.read_csv(
      EXAMPLE_VCF_PATH,
      skip_lines=header_line,
      separator="\t",
      schema_overrides={vcf2df.CHROM_COL: AVAILABLE_CHROMOSOMES},
  )
  info_df = vcf2df.split_info_column(sample_a008_df.select(vcf2df.INFO_COL), info_col_names)
  csq = vcf2df.split_csq_column(csq_col_names, info_df.select(vcf2df.CSQ_COL))
  # No COSMIC v98 rows should exist
  assert csq.filter(pl.col(vcf2df.EXISTING_VARIATION_COL).str.contains("COSV")).is_empty()
