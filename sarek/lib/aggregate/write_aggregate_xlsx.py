#!/opt/conda/bin/python
"""
Write the aggregate report data to an excel file.
"""
__docformat__ = "reStructuredText"
__author__ = "Brandon Sie, Neil A. Patel"

import polars as pl
import narwhals as nw
import xlsxwriter
from ..excel_style import get_column_widths, HEADER_STYLE
from .create_html_report import get_legend


def write_xlsx(
    file_name: str,
    metrics: nw.DataFrame,
    aggregate: nw.DataFrame,
    concordance: nw.DataFrame,
    mismatch: nw.DataFrame,
):
  """
  Export the aggregate report as a multi-tab xlsx.  
  """

  # Narwhals does not have a excel engine so this has to be done in polars
  pl_metrics = metrics.to_polars()
  pl_aggregate = aggregate.to_polars()
  pl_concordance = concordance.to_polars()
  pl_mismatch = mismatch.to_polars()
  pl_legend = get_legend().to_polars()

  metrics_pct_columns = [col for col in pl_metrics.columns if
                         ("pct" in col.lower()) or ("percent" in col.lower()) or ("%" in col.lower())]
  summary_pct_columns = [col for col in pl_aggregate.columns if
                         ("pct" in col.lower()) or ("percent" in col.lower()) or ("%" in col.lower())]

  with xlsxwriter.Workbook(file_name) as workbook:

    pl_metrics.write_excel(
        workbook=workbook,
        worksheet="metrics",
        column_widths=get_column_widths(pl_metrics, font_multiplier=8.0, first_row_multiplier=11.0),
        header_format=HEADER_STYLE,
        dtype_formats={pl.Int64: "0;[RED]-0"},
        column_formats={col: "0.00%" for col in metrics_pct_columns},
    )

    pl_aggregate.write_excel(
        workbook=workbook,
        worksheet="aggregate_summary",
        column_widths=get_column_widths(pl_aggregate, font_multiplier=8.0, first_row_multiplier=10.0),
        header_format=HEADER_STYLE,
        dtype_formats={pl.Int64: "0;[RED]-0"},
        column_formats={col: "0.00%" for col in summary_pct_columns},
    )

    fingerprint_issues_sheet = workbook.add_worksheet("fingerprint_issues")
    fingerprint_issues_sheet.write("A1", "Important sample match information:")
    fingerprint_issues_sheet.write("A2", "Relatedness threshold is 0.9")
    pl_mismatch.write_excel(
        workbook=workbook,
        worksheet=fingerprint_issues_sheet,
        column_widths=get_column_widths(pl_mismatch, font_multiplier=8.0, first_row_multiplier=10.0),
        header_format=HEADER_STYLE,
        dtype_formats={pl.Int64: "0;[RED]-0"},
        position="A4",
    )

    pl_concordance.write_excel(
        workbook=workbook,
        worksheet="FP_relatedness_matrix",
        dtype_formats={
            pl.Int64: "0;[RED]-0",
            pl.Float64: "[Red][<0.9]0.000;0.000",
        },
    )

    pl_legend.write_excel(
        workbook=workbook,
        worksheet="Legend",
        column_widths=get_column_widths(pl_legend, font_multiplier=8.0, first_row_multiplier=10.0),
        header_format=HEADER_STYLE,
    )
