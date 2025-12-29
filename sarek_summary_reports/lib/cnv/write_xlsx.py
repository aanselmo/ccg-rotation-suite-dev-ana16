#!/opt/conda/bin/python
"""
Write the CNV data to an excel file.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"

import io
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
    arm_summary: nw.DataFrame,
    mtc_hist_mem_loc: io.BytesIO,
    mad_hist_mem_loc: io.BytesIO,
):
  # Narwhals does not have a excel engine so this has to be done in polars
  pl_qc = qc_metrics.to_native()
  pl_cnv = cnv_variants.to_native()
  pl_gene = gene_summary.to_native()
  pl_arm = arm_summary.to_native()
  qc_sheet_widths = get_column_widths(pl_qc, first_row_multiplier=13.0)
  qc_sheet_widths[MTC_COL] = qc_sheet_widths[MTC_COL] / 3.0
  arm_sheet_widths = get_column_widths(pl_arm, first_row_multiplier=6.0)
  arm_sheet_widths[ARM_COL] = arm_sheet_widths[ARM_COL] * 5.0
  arm_sheet_widths[CHROMOSOME_COL] = arm_sheet_widths[CHROMOSOME_COL] * 2.0
  with xlsxwriter.Workbook(file_name) as workbook:
    qc_sheet = workbook.add_worksheet("QC Metrics")
    pl_qc.write_excel(
        workbook=workbook,
        worksheet=qc_sheet,
        column_widths=qc_sheet_widths,
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )
    pl_cnv.write_excel(
        workbook=workbook,
        worksheet="Copy Number Variants (CNV)",
        column_widths=get_column_widths(pl_cnv),
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )
    pl_gene.write_excel(
        workbook=workbook,
        worksheet="Gene Level Summary",
        column_widths=get_column_widths(pl_gene, first_row_multiplier=5.0),
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )
    pl_arm.write_excel(
        workbook=workbook,
        worksheet="Arm Level Summary",
        column_widths=arm_sheet_widths,
        header_format=HEADER_STYLE,
        dtype_formats={ pl.Int64: "0;[RED]-0" },
    )
    qc_sheet.insert_image("J2", "mtc_graph.png", {"image_data": mtc_hist_mem_loc})
    qc_sheet.insert_image("J30", "mad_graph.png", {"image_data": mad_hist_mem_loc})