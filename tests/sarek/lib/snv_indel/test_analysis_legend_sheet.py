#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the anlysis_legen_sheet.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.analysis_legend_sheet as analysis_legend_sheet


def test_get_legend():
  legend_df = analysis_legend_sheet.get_legend("GRCh37")
  assert legend_df.filter(nw.col(analysis_legend_sheet.METRIC_COL) == "REF").get_column(
      analysis_legend_sheet.DEFINITIONS_COL)[0] == "Reference Allele in GRCh37"
