#!/usr/local/miniforge3/bin/python

"""
Create a MultiQC report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


from enum import StrEnum
from pathlib import Path
import pandas as pd
import polars as pl
import narwhals.stable.v2 as nw
import report_creator as rc

RELATEDNESS_COL = "relatedness"
EXPECTED_RELATEDNESS_COL = "expected_relatedness"
MATCH_STATUS_COL = "Match Status"
IBS0_COL = "ibs0"
IBS2_COL = "ibs2"
SAMPLE_A_COL = "#sample_a"
SAMPLE_B_COL = "sample_b"
SAMPLE_PAIR_COL = "sample_pair"


class MATCH_STATUS(StrEnum):
  INAPPROPRIATE_MATCH = "Inappropriate Match"
  CORRECT_MATCH = "OK"
  MISMATCH = "Mismatch"


def parse_somalier_pairs_tsv(tsv_somalier_pairs: str) -> nw.DataFrame:
  somalier_data = pl.read_csv(tsv_somalier_pairs, separator="\t", schema_overrides={RELATEDNESS_COL: pl.Float64})
  somalier_data = somalier_data.select(
      pl.concat_list(SAMPLE_A_COL, SAMPLE_B_COL).alias(SAMPLE_PAIR_COL),
      pl.when((pl.col(RELATEDNESS_COL) > pl.lit(0.9)) & (pl.col(EXPECTED_RELATEDNESS_COL) == pl.lit(-1)))
        .then(pl.lit(MATCH_STATUS.INAPPROPRIATE_MATCH))
        .otherwise(pl.lit(MATCH_STATUS.CORRECT_MATCH)).alias(MATCH_STATUS_COL),
      pl.col(SAMPLE_A_COL, SAMPLE_B_COL, RELATEDNESS_COL, IBS0_COL, IBS2_COL, EXPECTED_RELATEDNESS_COL)
  )
  somalier_data = somalier_data.with_columns(
    pl.when((pl.col(RELATEDNESS_COL) < pl.lit(0.9)) & (pl.col(EXPECTED_RELATEDNESS_COL) == pl.lit(1)))
      .then(pl.lit(MATCH_STATUS.MISMATCH))
      .otherwise(pl.col(MATCH_STATUS_COL)).alias(MATCH_STATUS_COL)
  )
  return nw.from_native(somalier_data)

def get_legend():
  script_dir = Path(__file__).resolve().parent
  legend = nw.read_csv(f"{script_dir}/aggregate_legend.csv", backend="pandas")
  return legend


def get_fingerprinting_block(tsv_somalier_pairs: str):
  somalier_data = parse_somalier_pairs_tsv(tsv_somalier_pairs)
  pd_somalier_data = somalier_data.to_pandas()
  somalier_html_block = rc.Block(
      rc.Text(
          ("This data shows which samples should/shouldn't be matching. However, the coverage"
           + " was not high enough confidently run fingerprinting (with Somalier)."),
          label="Fingerprinting Data",
      ),
  )
  if not pd_somalier_data.loc[:, RELATEDNESS_COL].dropna().empty:
    somalier_html_block = rc.Block(
        rc.Text(
            "This data shows which samples should/shouldn't be matching.",
            label="Fingerprinting Data",
        ),
        rc.DataTable(
            pd_somalier_data.loc[
                pd_somalier_data.loc[:, MATCH_STATUS_COL] == MATCH_STATUS.INAPPROPRIATE_MATCH,
                [SAMPLE_A_COL, SAMPLE_B_COL, RELATEDNESS_COL, IBS0_COL, IBS2_COL]
            ],
            label="Inappropriate Matches",
            index=False,
        ),
        rc.DataTable(
            pd_somalier_data.loc[
                pd_somalier_data.loc[:, MATCH_STATUS_COL] == MATCH_STATUS.MISMATCH,
                [SAMPLE_A_COL, SAMPLE_B_COL, RELATEDNESS_COL, IBS0_COL, IBS2_COL]
            ],
            label="Mismatches",
            index=False,
        ),
        rc.Select(
            blocks=[
                rc.Scatter(
                    pd_somalier_data,
                    x=IBS2_COL,
                    y=RELATEDNESS_COL,
                    label="Scatter Plot: IBS2 v Relatedness",
                    dimension=MATCH_STATUS_COL,
                    hover_name=SAMPLE_PAIR_COL
                ),
                rc.Scatter(
                    pd_somalier_data,
                    x=IBS0_COL,
                    y=RELATEDNESS_COL,
                    dimension=MATCH_STATUS_COL,
                    label="Scatter Plot: IBS0 v Relatedness",
                    hover_name=SAMPLE_PAIR_COL
                ),
            ]
        ),
    )
  return somalier_html_block


def create_html_report(summary_table: pd.DataFrame, tsv_somalier_pairs: str, project_code: str, output_file_name: str):
  somalier_html_block = get_fingerprinting_block(tsv_somalier_pairs)
  legend = get_legend().to_native()
  with rc.ReportCreator(
      title=f"{project_code} Aggregation Metrics Report",
      description="A summary of aggregation metrics and a fingerprinting result report.",
      footer="Created by CCG's Bioinformaticians (Neil A. Patel, & Anthony Anselmo)",
  ) as report:
    view = rc.Block(
        rc.DataTable(
            summary_table,
            label="Aggregation Metrics Summary",
            index=False,
        ),
        rc.Separator(),
        somalier_html_block,  
        rc.Separator(),
        rc.DataTable(
            legend,
            label="Legend",
            index=False,
        )
    )
  report.save(view, output_file_name)
