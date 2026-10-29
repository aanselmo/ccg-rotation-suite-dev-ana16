#!/opt/conda/bin/python
"""
Write the CNV data to an excel file.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"

import io
from pathlib import Path
import polars as pl
import narwhals as nw
import xlsxwriter
from openpyxl import load_workbook
from ..excel_style import get_column_widths, HEADER_STYLE
from .qc_sheet import MTC_COL, THRESHOLD_LEGEND, THRESHOLD_LEGEND2, THRESHOLD_LEGEND3
# from .chromosome_sheet import ARM_COL, CHROMOSOME_COL  # Arm Level Summary tab disabled

INT_FORMAT = {pl.Int64: "0;[RED]-0"}

# openpyxl -> xlsxwriter style name mappings
_H_ALIGN = {"left": "left", "center": "center", "right": "right", "fill": "fill", "justify": "justify",
            "centerContinuous": "center_across", "distributed": "distributed"}
_V_ALIGN = {"top": "top", "center": "vcenter", "bottom": "bottom", "justify": "vjustify",
            "distributed": "vdistributed"}
_BORDER_STYLE = {"thin": 1, "medium": 2, "dashed": 3, "dotted": 4, "thick": 5, "double": 6, "hair": 7,
                 "mediumDashed": 8, "dashDot": 9, "mediumDashDot": 10, "dashDotDot": 11,
                 "mediumDashDotDot": 12, "slantDashDot": 13}


def _rgb(color) -> str | None:
  """
  Convert an openpyxl color to "#RRGGBB". Theme/indexed colors can't be resolved without the
  source theme, so they return None and fall back to Excel's default.
  """
  if color is None or color.type != "rgb" or not isinstance(color.rgb, str):
    return None
  return f"#{color.rgb[-6:]}"


def _cell_format(workbook: xlsxwriter.Workbook, cell, cache: dict):
  """
  Build (or reuse) an xlsxwriter format matching an openpyxl cell's font, fill, alignment,
  borders, and number format. Returns None when the cell has no styling worth copying.
  """
  props = {}
  font = cell.font
  if font.b:
    props["bold"] = True
  if font.i:
    props["italic"] = True
  if font.u:
    props["underline"] = 1
  if font.sz and font.sz != 11:
    props["font_size"] = font.sz
  if font.name and font.name != "Calibri":
    props["font_name"] = font.name
  if (font_color := _rgb(font.color)):
    props["font_color"] = font_color

  if cell.fill.fill_type == "solid" and (bg_color := _rgb(cell.fill.fgColor)):
    props["bg_color"] = bg_color

  alignment = cell.alignment
  if alignment.wrap_text:
    props["text_wrap"] = True
  if alignment.horizontal in _H_ALIGN:
    props["align"] = _H_ALIGN[alignment.horizontal]
  if alignment.vertical in _V_ALIGN:
    props["valign"] = _V_ALIGN[alignment.vertical]

  for side in ("left", "right", "top", "bottom"):
    border_side = getattr(cell.border, side)
    if border_side is not None and border_side.style in _BORDER_STYLE:
      props[side] = _BORDER_STYLE[border_side.style]
      if (border_color := _rgb(border_side.color)):
        props[f"{side}_color"] = border_color

  if cell.number_format and cell.number_format != "General":
    props["num_format"] = cell.number_format

  if not props:
    return None
  key = tuple(sorted(props.items()))
  if key not in cache:
    cache[key] = workbook.add_format(props)
  return cache[key]


def add_reference_sheets(workbook: xlsxwriter.Workbook, reference_files: list[Path | str]) -> None:
  """
  Copy every tab from each reference .xlsx into the report, keeping the tab name, cell values,
  cell formatting, merged cells, column widths, row heights, frozen panes, and filters.

  :param workbook xlsxwriter.Workbook: The open report workbook
  :param reference_files list[Path | str]: Paths to the .xlsx files whose tabs should be copied
  """
  format_cache: dict = {}
  for path in reference_files:
    # data_only=True -> formulas come across as their last calculated values
    src_wb = load_workbook(path, data_only=True)
    for src_ws in src_wb.worksheets:
      dest_ws = workbook.add_worksheet(src_ws.title)

      # Merged ranges are written with merge_range and skipped in the cell loop below
      merged_cells = set()
      for merged_range in src_ws.merged_cells.ranges:
        min_col, min_row, max_col, max_row = merged_range.bounds
        top_left = src_ws.cell(row=min_row, column=min_col)
        dest_ws.merge_range(
            min_row - 1, min_col - 1, max_row - 1, max_col - 1,
            "" if top_left.value is None else top_left.value,
            _cell_format(workbook, top_left, format_cache),
        )
        merged_cells.update(
            (row, col) for row in range(min_row, max_row + 1) for col in range(min_col, max_col + 1)
        )

      for row in src_ws.iter_rows():
        for cell in row:
          if (cell.row, cell.column) in merged_cells:
            continue
          cell_format = _cell_format(workbook, cell, format_cache)
          if cell.value is not None:
            dest_ws.write(cell.row - 1, cell.column - 1, cell.value, cell_format)
          elif cell_format is not None:  # empty but styled (e.g. a colored legend cell)
            dest_ws.write_blank(cell.row - 1, cell.column - 1, None, cell_format)

      for dim in src_ws.column_dimensions.values():
        if dim.width and dim.min and dim.max:
          dest_ws.set_column(dim.min - 1, dim.max - 1, dim.width)

      for row_idx, dim in src_ws.row_dimensions.items():
        if dim.height:
          dest_ws.set_row(row_idx - 1, dim.height)

      if src_ws.freeze_panes:
        dest_ws.freeze_panes(src_ws.freeze_panes)

      if src_ws.auto_filter.ref:  # header-row filter dropdowns
        dest_ws.autofilter(src_ws.auto_filter.ref)


def write_xlsx(
    file_name: str,
    qc_metrics: nw.DataFrame,
    cnv_variants: nw.DataFrame,
    gene_summary: nw.DataFrame,
    # arm_summary: nw.DataFrame,  # Arm Level Summary tab disabled
    armlevel_table: nw.DataFrame,
    mtc_hist_mem_loc: io.BytesIO | None = None,
    mad_hist_mem_loc: io.BytesIO | None = None,
    reference_files: list[Path | str] | None = None,
):
  """
  Write every CNV sheet to one Excel workbook.

  A sheet is skipped when its DataFrame has no columns (the empty placeholder
  main() uses when that part of the report had no input files). The QC graphs
  are skipped when their image is None. Tabs from reference_files are copied in
  after the data sheets, keeping their original names.
  """
  # Narwhals does not have a excel engine so this has to be done in polars
  pl_qc = qc_metrics.to_native()
  pl_cnv = cnv_variants.to_native()
  pl_gene = gene_summary.to_native()
  # pl_arm = arm_summary.to_native()  # Arm Level Summary tab disabled
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

    # Arm Level Summary tab disabled
    # if pl_arm.width:
    #   arm_sheet_widths = get_column_widths(pl_arm, first_row_multiplier=6.0)
    #   arm_sheet_widths[ARM_COL] = arm_sheet_widths[ARM_COL] * 5.0
    #   arm_sheet_widths[CHROMOSOME_COL] = arm_sheet_widths[CHROMOSOME_COL] * 2.0
    #   pl_arm.write_excel(
    #       workbook=workbook,
    #       worksheet="Arm Level Summary",
    #       column_widths=arm_sheet_widths,
    #       header_format=HEADER_STYLE,
    #       dtype_formats=INT_FORMAT,
    #   )

    if pl_armlevel.width:
      pl_armlevel.write_excel(
          workbook=workbook,
          worksheet="Arm Level Calls",
          column_widths=get_column_widths(pl_armlevel),
          header_format=HEADER_STYLE,
          dtype_formats=INT_FORMAT,
          freeze_panes=(1, 2),  # keep the header row and the sample/arm columns visible while scrolling
      )

    if reference_files:
      add_reference_sheets(workbook, reference_files)
