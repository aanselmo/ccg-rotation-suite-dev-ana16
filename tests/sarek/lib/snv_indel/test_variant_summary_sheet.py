#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the variant_summary_sheet.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
from os.path import join
from pathlib import Path
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.snv_indel_sheet as snv_indel_sheet
import sarek.lib.snv_indel.variant_summary_sheet as variant_summary_sheet


EXAMPLE_VCF_PATH = join(Path(__file__).parent, "data", "A008_X108_58797-1-C1.mutect2.filtered_VEP.ann.vcf.gz")


def test_get_variant_summary():
  all_variant_vcfs = [EXAMPLE_VCF_PATH]
  all_variants = snv_indel_sheet.get_variants(all_variant_vcfs)
  variant_summary_table = variant_summary_sheet.get_variant_summary(all_variants)
  assert variant_summary_table.filter(
      nw.col(snv_indel_sheet.SAMPLE_COL) == "A008_X108_58797-1-C1"
  ).get_column(variant_summary_sheet.BOTH_COL)[0] == 1
