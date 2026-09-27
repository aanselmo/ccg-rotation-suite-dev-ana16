#!/usr/local/miniforge3/bin/python

"""
Creates additional tables for PDX samples that have:
1. metrics on mouse contamination.
2. blacklisted Human Aligned Mouse Alleles (HAMAs) that should be filtered out.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import sys
from os.path import basename, dirname
from pathlib import Path
import polars as pl
import narwhals.stable.v2 as nw
script_dir = Path(__file__).resolve().parent
sarek_script_top_dir = dirname(dirname(script_dir))
sys.path.append(sarek_script_top_dir)
from lib.snv_indel.snv_indel_sheet import SAMPLE_COL
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

PRIMARY_READS_COL = "primary reads"
NCBIM37_READS_COL = "NCBIM37 reads"
LAST_ITEM = -1


def get_contamination_metrics(mouse_bbsplit_file: str) -> nw.DataFrame:
  # Number of chunks are define by https://nf-co.re/sarek/3.7.1/docs/usage/#split-fastq-files
  # Typically, there are 12 chunks that a fastq is split into.
  CHUNK_NUM_COL = "chunk_num"
  LANE_COL = "lane"
  contamination = pl.read_csv(mouse_bbsplit_file)
  contamination = contamination.rename({
      "Sample": SAMPLE_COL,
      "primary": PRIMARY_READS_COL,
      "NCBIM37": NCBIM37_READS_COL
  })
  # Example original "Sample" col = "0001.CHGFR0004-PDX-B1-L001"
  # 0001 is the chunk and L001 is the lane entered in the sarek samplesheet.
  contamination = contamination.with_columns(
      pl.col(SAMPLE_COL).str.split_exact(".", 1).struct.rename_fields(CHUNK_NUM_COL, SAMPLE_COL)
  )
  contamination = contamination.with_columns(
      pl.col(SAMPLE_COL).str.split("-")
  )
  contamination = contamination.with_columns(
      pl.col(SAMPLE_COL).list.item(LAST_ITEM).alias(LANE_COL),
      pl.col(SAMPLE_COL).list.head(LAST_ITEM).list.join("-").alias(SAMPLE_COL)
  )
