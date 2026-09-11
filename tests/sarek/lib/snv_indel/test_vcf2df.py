#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the vcf2df.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
import sys
from os.path import dirname, abspath, join
from os import makedirs
import polars as pl
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.vcf2df as vcf2df
from sarek.lib.common import AVAILABLE_CHROMOSOMES


GIT_REPO_DIR = abspath(dirname(dirname(dirname(dirname(dirname(__file__))))))
EXAMPLE_VCF_PATH = join(GIT_REPO_DIR, "tests", "data", "A008_X108_58797-1-C1.mutect2.filtered_VEP.ann.vcf.gz")


def test_convert_vcf_to_data_frame():
  sample_a008_df = vcf2df.convert_vcf_to_data_frame(EXAMPLE_VCF_PATH)
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
  data_col_of_format_header = vcf2df.get_format_data_column_name(sample_a008_df)
  format_data = vcf2df.split_format_column(
      sample_a008_df.select(pl.col(vcf2df.FORMAT_COL), pl.col(data_col_of_format_header))
  )
  sample_a008_df = sample_a008_df.with_columns(pl.row_index())
  row = sample_a008_df.filter(
      pl.col(vcf2df.POSITION_COL) == 135526367
  ).get_column("index")[0]
  format_data = format_data.with_columns(pl.row_index())
  assert format_data.filter(pl.col("index") == row).get_column("DP")[0] == "355"
