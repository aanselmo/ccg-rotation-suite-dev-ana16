#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the lego_plot_sheet.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
from os.path import dirname, abspath, join
import polars as pl
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.lego_plot_sheet as lego_plot_sheet
from sarek.lib.snv_indel.vcf2df import POSITION_COL
from sarek.lib.common import AVAILABLE_CHROMOSOMES


GIT_REPO_DIR = abspath(dirname(dirname(dirname(dirname(dirname(__file__))))))
VCF_A016_PATH = join(GIT_REPO_DIR, "tests", "data", "A016_X116_BS07G52324-G1.mutect2.vcf.gz")
EXAMPLE_BED = join(GIT_REPO_DIR, "tests", "data", "POPv3_AllTargetedRegions_20251230.bed")


def test_get_variants():
  all_variant_vcfs = [VCF_A016_PATH]
  plot_memory_location = lego_plot_sheet.graph_lego_plot(all_variant_vcfs, "SCB1118c", "GRCh37", bed_file=EXAMPLE_BED)
  assert plot_memory_location.getbuffer().nbytes > 0
