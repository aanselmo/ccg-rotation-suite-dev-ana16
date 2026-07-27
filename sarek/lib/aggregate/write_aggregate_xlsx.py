#!/opt/conda/bin/python
"""
Write the aggregate report data to an excel file.
"""
__docformat__ = "reStructuredText"
__author__ = "Brandon Sie, Neil A. Patel"

import pandas as pd
import numpy as np
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
    num_somalier_sites
):
  """
  Export the aggregate report as a multi-tab xlsx.  
  """

  # Narwhals does not have a excel engine so this has to be done in polars
  pl_metrics = metrics.to_polars()
  pl_aggregate = aggregate.to_polars()
  pd_concordance = concordance.to_pandas() # need pandas for this data frame to fill null values with text
  pl_mismatch = mismatch.to_polars()
  pl_legend = get_legend().to_polars()


  metrics_pct_columns = [col for col in pl_metrics.columns if
                         ("pct" in col.lower()) or ("percent" in col.lower()) or ("%" in col.lower())]
  summary_pct_columns = [col for col in pl_aggregate.columns if
                         ("pct" in col.lower()) or ("percent" in col.lower()) or ("%" in col.lower())]
  
  with pd.ExcelWriter(file_name, engine="xlsxwriter") as writer:
    workbook  = writer.book
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
    fingerprint_issues_sheet.write("A2", "Somalier relatedness threshold is 0.9")
    fingerprint_issues_sheet.write("A3", "Informative site calls made with requirement of minimum 20x sequencing depth")
    
    pl_mismatch.write_excel(
        workbook=workbook,
        worksheet=fingerprint_issues_sheet,
        column_widths=get_column_widths(pl_mismatch, font_multiplier=8.0, first_row_multiplier=10.0),
        header_format=HEADER_STYLE,
        dtype_formats={pl.Int64: "0;[RED]-0"},
        position="A5",
    )

    styler_concordance = (
        pd_concordance.style
        .format("{:.3f}")
        .map(pd_color_rule)
        .apply_index(align_left)
    )

    relatedness_matrix_sheet = workbook.add_worksheet("FP_relatedness_matrix")
    matrix_sheet_name = "FP_relatedness_matrix"
    styler_concordance.to_excel(writer, sheet_name=matrix_sheet_name, startrow=2, startcol=0, index=False, na_rep="#NA")  # A3    
    relatedness_matrix_sheet.write("A1", f"Up to {num_somalier_sites} sites used for Somalier fingerprinting")

    pl_legend.write_excel(
        workbook=workbook,
        worksheet="Legend",
        column_widths=get_column_widths(pl_legend, font_multiplier=8.0, first_row_multiplier=10.0),
        header_format=HEADER_STYLE,
    )
  
  print(f"Output written to: {file_name}")


def pd_color_rule(v):
    """
    Color formatting helper for concordance matrix
    """
    if pd.isna(v):
        return "color: black"
    try:
        v = float(v)
    except Exception:
        return "color: black"
    return "color: red" if v < 0.9 else "color: black"


def align_left(x):
    return ['text-align: left' for x in x]