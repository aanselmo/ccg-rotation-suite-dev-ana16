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
from .qc_sheet import MTC_COL, THRESHOLD_LEGEND, THRESHOLD_LEGEND2, THRESHOLD_LEGEND3
from .chromosome_sheet import ARM_COL, CHROMOSOME_COL

INT_FORMAT = {pl.Int64: "0;[RED]-0"}


def write_xlsx(
    file_name: str,
    qc_metrics: nw.DataFrame,
    cnv_variants: nw.DataFrame,
    gene_summary: nw.DataFrame,
    arm_summary: nw.DataFrame,
    armlevel_table: nw.DataFrame,
    mtc_hist_mem_loc: io.BytesIO | None = None,
    mad_hist_mem_loc: io.BytesIO | None = None,
):
  """
  Write every CNV sheet to one Excel workbook.

  A sheet is skipped when its DataFrame has no columns (the empty placeholder
  main() uses when that part of the report had no input files). The QC graphs
  are skipped when their image is None.
  """
  # Narwhals does not have a excel engine so this has to be done in polars
  pl_qc = qc_metrics.to_native()
  pl_cnv = cnv_variants.to_native()
  pl_gene = gene_summary.to_native()
  pl_arm = arm_summary.to_native()
  pl_armlevel = armlevel_table.to_native()

  with xlsxwriter.Workbook(file_name) as workbook:
    if pl_qc.width:
      qc_sheet_widths = get_column_widths(pl_qc, first_row_multiplier=13.0)
      qc_sheet_widths[MTC_COL] = qc_sheet_widths[MTC_COL] / 3.0
      qc_sheet = workbook.add_worksheet("QC Metrics")
      pl_qc.write_excel(
          workbook=workbook,
          worksheet=qc_sheet,
          column_widths=qc_sheet_widths,
          header_format=HEADER_STYLE,
          dtype_formats=INT_FORMAT,
      )
      if mtc_hist_mem_loc is not None:
        qc_sheet.insert_image("J2", "mtc_graph.png", {"image_data": mtc_hist_mem_loc})
      if mad_hist_mem_loc is not None:
        qc_sheet.insert_image("J30", "mad_graph.png", {"image_data": mad_hist_mem_loc})
      qc_sheet.write("T2", THRESHOLD_LEGEND)
      qc_sheet.write("T3", THRESHOLD_LEGEND2)
      qc_sheet.write("T4", THRESHOLD_LEGEND3)

    if pl_cnv.width:
      pl_cnv.write_excel(
          workbook=workbook,
          worksheet="Copy Number Variants (CNV)",
          column_widths=get_column_widths(pl_cnv),
          header_format=HEADER_STYLE,
          dtype_formats=INT_FORMAT,
      )

    if pl_gene.width:
      pl_gene.write_excel(
          workbook=workbook,
          worksheet="Gene Level Summary",
          column_widths=get_column_widths(pl_gene, first_row_multiplier=5.0),
          header_format=HEADER_STYLE,
          dtype_formats=INT_FORMAT,
      )

    if pl_arm.width:
      arm_sheet_widths = get_column_widths(pl_arm, first_row_multiplier=6.0)
      arm_sheet_widths[ARM_COL] = arm_sheet_widths[ARM_COL] * 5.0
      arm_sheet_widths[CHROMOSOME_COL] = arm_sheet_widths[CHROMOSOME_COL] * 2.0
      pl_arm.write_excel(
          workbook=workbook,
          worksheet="Arm Level Summary",
          column_widths=arm_sheet_widths,
          header_format=HEADER_STYLE,
          dtype_formats=INT_FORMAT,
      )

    if pl_armlevel.width:
      pl_armlevel.write_excel(
          workbook=workbook,
          worksheet="Arm Level Calls",
          column_widths=get_column_widths(pl_armlevel),
          header_format=HEADER_STYLE,
          dtype_formats=INT_FORMAT,
          freeze_panes=(1, 2),  # keep the header row and the sample/arm columns visible while scrolling
      )
