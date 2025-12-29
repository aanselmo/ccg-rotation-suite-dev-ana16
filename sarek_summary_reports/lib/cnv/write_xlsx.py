#!/opt/conda/bin/python
"""
Write the CNV data to an excel file.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"

import polars as pl
import narwhals as nw
import xlsxwriter
from ..excel_style import get_column_widths, HEADER_STYLE
from .qc_sheet import MTC_COL
from .chromosome_sheet import ARM_COL, CHROMOSOME_COL

def write_xlsx(
    file_name: str,
    qc_metrics: nw.DataFrame,
    cnv_variants: nw.DataFrame,
    gene_summary: nw.DataFrame,
    arm_summary: nw.DataFrame
):
  # Narwhals does not have a excel engine so this has to be done in polars
  qc_sheet = qc_metrics.to_native()
  cnv_sheet = cnv_variants.to_native()
  gene_sheet = gene_summary.to_native()
  arm_sheet = arm_summary.to_native()
  qc_sheet_widths = get_column_widths(qc_sheet, first_row_multiplier=13.0)
  qc_sheet_widths[MTC_COL] = qc_sheet_widths[MTC_COL] / 3.0
  arm_sheet_widths = get_column_widths(arm_sheet, first_row_multiplier=6.0)
  arm_sheet_widths[ARM_COL] = arm_sheet_widths[ARM_COL] * 5.0
  arm_sheet_widths[CHROMOSOME_COL] = arm_sheet_widths[CHROMOSOME_COL] * 2.0
  with xlsxwriter.Workbook(file_name) as workbook:
    qc_sheet.write_excel(
        workbook=workbook,
        worksheet="QC Metrics",
        column_widths=qc_sheet_widths,
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )
    cnv_sheet.write_excel(
        workbook=workbook,
        worksheet="Copy Number Variants (CNV)",
        column_widths=get_column_widths(cnv_sheet),
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )
    gene_sheet.write_excel(
        workbook=workbook,
        worksheet="Gene Level Summary",
        column_widths=get_column_widths(gene_sheet, first_row_multiplier=5.0),
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )
    arm_sheet.write_excel(
        workbook=workbook,
        worksheet="Arm Level Summary",
        column_widths=arm_sheet_widths,
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )