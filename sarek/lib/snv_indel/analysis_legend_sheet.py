#!/usr/local/miniforge3/bin/python

"""
Creates the 'Legend' Sheet.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


from os.path import join
from pathlib import Path
import polars as pl
import narwhals.stable.v2 as nw


DEFINITIONS_COL = "Definition"
SOURCE_COL = "Source"


def get_legend(reference_genome: str, has_pdx=False) -> nw.DataFrame:
  """
  Get's the 'Legend' Sheet in the snv_indel report as a DataFrame.
  :param reference_genome str: The reference genome used in the report.
  :return nw.DataFrame: The legend sheet as a DataFrame.
  """
  script_dir = Path(__file__).resolve().parent
  legend_df = nw.from_native(pl.read_csv(join(script_dir, "legend.csv")))
  if not has_pdx:
    legend_df.filter(pl.col(SOURCE_COL) != "bbsplit")
  legend_df = legend_df.with_columns(
      nw.col(DEFINITIONS_COL).str.replace("REFERENCE_GENOME", reference_genome).alias(DEFINITIONS_COL)
  )
  return legend_df
