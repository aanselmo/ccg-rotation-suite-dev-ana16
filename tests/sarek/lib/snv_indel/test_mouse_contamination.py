#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the mouse_contamination.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
from os.path import dirname, abspath, join
from pathlib import Path
import polars as pl
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.vcf2df as vcf2df
import sarek.lib.snv_indel.mouse_contamination as mouse_contamination


EXAMPLE_CONTAMINATION_METRICS_FILE = join(Path(__file__).parent, "data", "bbmap-bbsplit_plot.txt")


def test_get_contamination_metrics():
  metrics = mouse_contamination.get_contamination_metrics(EXAMPLE_CONTAMINATION_METRICS_FILE)
  assert (
      metrics.filter(nw.col(mouse_contamination.SAMPLE_COL) == "CEPH1408-H1")
      .get_column(mouse_contamination.PRIMARY_READS_COL)[0] == 27906340
  )
  assert (
      metrics.filter(nw.col(mouse_contamination.SAMPLE_COL) == "CEPH1408-H1")
      .get_column(mouse_contamination.CONTAMINATION_PCT_COL)[0] < 0.01
  )
