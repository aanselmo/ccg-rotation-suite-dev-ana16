#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the snv_indel_sheet.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
from os.path import join
from pathlib import Path
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.snv_indel_sheet as snv_indel_sheet
from sarek.lib.snv_indel.vcf2df import POSITION_COL
from sarek.lib.common import AVAILABLE_CHROMOSOMES


EXAMPLE_VCF_PATH = join(Path(__file__).parent, "data", "A008_X108_58797-1-C1.mutect2.filtered_VEP.ann.vcf.gz")


def test_get_variants():
  all_variant_vcfs = [EXAMPLE_VCF_PATH]
  all_variants = snv_indel_sheet.get_variants(all_variant_vcfs)
  # filtration columns should have been added.
  assert snv_indel_sheet.REVIEW_REQUIRED_COL in all_variants.columns
  assert snv_indel_sheet.RECOMMENDED_FILTER_COL in all_variants.columns
  # Assert 1 variant has the correct alternate allele
  assert all_variants.filter(
      nw.col(POSITION_COL) == 41266124
  ).get_column("ALT")[0] == "G"
  # Assert that only 1 variant/row makes it through filtration
  assert all_variants.shape[0] == 1
