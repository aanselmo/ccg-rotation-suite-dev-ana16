#!/usr/local/miniforge3/bin/python

"""
Create a MultiQC report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import polars as pl
import narwhals.stable.v2 as nw
import report_creator as rc


RELATEDNESS_COL = "relatedness"
IBS0_COL = "ibs0"
IBS2_COL = "ibs2"
SAMPLE_A_COL = "#sample_a"
SAMPLE_B_COL = "sample_b"
SAMPLE_PAIR_COL = "sample_pair"


def parse_somalier_pairs_tsv(tsv_somalier_pairs: str) -> nw.DataFrame:
  somalier_data = pl.read_csv(tsv_somalier_pairs, separator="\t")
  somalier_data = somalier_data.select(
      pl.concat_list(SAMPLE_A_COL, SAMPLE_B_COL).alias(SAMPLE_PAIR_COL),
      pl.col(SAMPLE_A_COL, SAMPLE_B_COL, RELATEDNESS_COL, IBS0_COL, IBS2_COL)
  )
  return nw.from_native(somalier_data)


def create_html_report(summary_table: nw.DataFrame, tsv_somalier_pairs: str, project_code: str, output_file_name: str):
  somalier_data = parse_somalier_pairs_tsv(tsv_somalier_pairs)
  pd_somalier_data = somalier_data.to_pandas()
  pd_summary_table = summary_table.to_pandas()
  with rc.ReportCreator(
      title=f"{project_code} Aggregation Metrics Report",
      description="A summary of aggregation metrics and a fingerprinting result report.",
      footer="Created by CCG's Bioinformaticians (Neil A. Patel, & Anthony Anselmo)",
  ) as report:
    view = rc.Block(
        rc.DataTable(
            pd_summary_table,
            label="Aggregation Metrics Summary",
            index=True,
        ),
        rc.Separator(),
        rc.DataTable(
            pd_somalier_data.loc[
                pd_somalier_data.loc[:, RELATEDNESS_COL] > 0.9,
                [SAMPLE_A_COL, SAMPLE_B_COL, RELATEDNESS_COL, IBS0_COL, IBS2_COL]
            ],
            label="Matching Samples",
            index=True,
        ),
        rc.Select(
            blocks=[
                rc.Scatter(
                    pd_somalier_data,
                    x=IBS2_COL,
                    y=RELATEDNESS_COL,
                    label="Scatter Plot: IBS2 v Relatedness",
                    hover_name=SAMPLE_PAIR_COL
                ),
                rc.Scatter(
                    pd_somalier_data,
                    x=IBS0_COL,
                    y=RELATEDNESS_COL,
                    label="Scatter Plot: IBS0 v Relatedness",
                    hover_name=SAMPLE_PAIR_COL
                ),
            ]
        ),
    )
  report.save(view, output_file_name)
