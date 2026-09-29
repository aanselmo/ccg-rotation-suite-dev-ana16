#!/opt/conda/bin/python
"""
Write the SNV Indel data to an excel file.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"

import io
import math
import string
from pip._vendor.rich.progress import Progress
import polars as pl
import narwhals as nw
import xlsxwriter
from .vcf2df import (
    EXISTING_VARIATION_COL, MAX_AF_POPS_COL, CHROM_COL, POSITION_COL, COSMIC_COLS,
    BEST_EFFECT_CONSEQUENCE_COL, BEST_EFFECT_HGVS_C_COL, BEST_EFFECT_HGVS_P_COL,
    CANONICAL_CONSEQUENCE_COL, CANONICAL_HGVS_C_COL, CANONICAL_HGVS_P_COL
)
from .snv_indel_sheet import (
    SAMPLE_COL, REF_COL, ALT_COL, RECOMMENDED_FILTER_COL, MUTECT2_FILTER_COL
)
from .mouse_contamination import CONTAMINATION_PCT_COL
from .variant_summary_sheet import REPEAT_COL_SUFFIX
from ..excel_style import get_column_widths, HEADER_STYLE, PERCENT_FORMAT


MAX_ROWS_EXCEL = 1048576
SNV_INDEL_SHEET = "SNVs and Indels"
VARIANT_SUMMARY_SHEET = "Variant Summary Table"
LEGEND_SHEET = "Legend"


def write_by_cell(df: pl.DataFrame, worksheet: xlsxwriter.worksheet, column_formats: dict):
  """
  This method is about 13% slower than writing by row. However, the advantage of this is that you can add column/cell
  formats.
  """
  with Progress() as progress:
    write_progress_bar = progress.add_task(f"Writing worksheet {worksheet.get_name()}...", total=df.shape[0])
    rows_for_1_percent = math.ceil(df.shape[0] / 100)
    for row_number, row in enumerate(df.iter_rows(named=True), start=1):
      for column_number, column in enumerate(row.keys(), start=0):
        if column in column_formats:
          worksheet.write(row_number, column_number, row[column], column_formats[column])
        else:
          worksheet.write(row_number, column_number, row[column])
      if (row_number % rows_for_1_percent) == 0:
        progress.update(write_progress_bar, completed=row_number)


def write_by_row(df: pl.DataFrame, worksheet: xlsxwriter.worksheet):
  """
  This method is 13% faster than writing by cell. However, you cannot use column/cell formats.
  """
  with Progress() as progress:
    write_progress_bar = progress.add_task(f"Writing worksheet {worksheet.get_name()}...", total=df.shape[0])
    rows_for_1_percent = math.ceil(df.shape[0] / 100)
    for row_number, row in enumerate(df.iter_rows(), start=1):
      worksheet.write_row(row_number, 0, row)
      if (row_number % rows_for_1_percent) == 0:
        progress.update(write_progress_bar, completed=row_number)


def write_df_to_excel(df: pl.DataFrame, worksheet: xlsxwriter.worksheet, column_widths: dict, header_format,
                      user_autofilter=True, column_formats={}):
  """
  Because I'm using {"constant_memory": True} as an option in xlsxwriter, I cannot use polars.write_excel.
  I have to do it via loop. I'm doing this to support over 1M rows in excel.
  https://xlsxwriter.readthedocs.io/working_with_memory.html
  """
  for column_number, one_width in enumerate(column_widths.values()):
    worksheet.set_column_pixels(column_number, column_number, one_width)
  worksheet.write_row(0, 0, df.columns, cell_format=header_format)
  if column_formats:
    write_by_cell(df, worksheet, column_formats)
  else:
    write_by_row(df, worksheet)
  if user_autofilter:
    worksheet.autofilter(0, 0, df.shape[0] - 1, df.shape[1] - 1)


def get_samples_per_sheet(pl_variants: pl.DataFrame) -> tuple[list[list[str]], bool]:
  samples_to_split = []
  samples_per_sheet = [[]]
  rows_in_sheet = 1  # start with 1 for header
  sheet_number = 0
  COUNT_COL = "count"
  rows_per_sample = pl_variants.group_by(pl.col(SAMPLE_COL)).agg(pl.count(POSITION_COL).alias(COUNT_COL))
  rows_per_sample = rows_per_sample.sort(by=SAMPLE_COL)
  samples_to_split = rows_per_sample.filter(pl.col(COUNT_COL) >= MAX_ROWS_EXCEL).get_column(SAMPLE_COL).to_list()
  for row in rows_per_sample.iter_rows(named=True):
    rows_in_sheet = rows_in_sheet + row[COUNT_COL]
    if rows_in_sheet <= MAX_ROWS_EXCEL:
      samples_per_sheet[sheet_number].append(row[SAMPLE_COL])
    else:
      if len(samples_per_sheet[sheet_number]) > 0:
        sheet_number = sheet_number + 1
      rows_in_sheet = row[COUNT_COL] + 1
      if (len(samples_per_sheet) - 1) == sheet_number:
        samples_per_sheet[sheet_number].append(row[SAMPLE_COL])
      else:
        samples_per_sheet.append([row[SAMPLE_COL]])
  return (samples_per_sheet, samples_to_split)


def write_all_snv_indel_sheets(workbook: xlsxwriter.Workbook, pl_variants: pl.DataFrame, sheet_basename: str,
                               header_format):
  pl_variants = pl_variants.sort(by=[SAMPLE_COL, CHROM_COL, POSITION_COL])
  pl_variants_widths = get_column_widths(pl_variants, first_row_multiplier=13.0)
  for one_col in (COSMIC_COLS + [
      REF_COL, ALT_COL, RECOMMENDED_FILTER_COL, MUTECT2_FILTER_COL,
      EXISTING_VARIATION_COL, MAX_AF_POPS_COL
  ]):
    pl_variants_widths[one_col] = 130.0
  for one_col in [
      BEST_EFFECT_CONSEQUENCE_COL, BEST_EFFECT_HGVS_C_COL, BEST_EFFECT_HGVS_P_COL,
      CANONICAL_CONSEQUENCE_COL, CANONICAL_HGVS_C_COL, CANONICAL_HGVS_P_COL,
  ]:
    pl_variants_widths[one_col] = 250.0
  (samples_per_sheet, samples_to_split) = get_samples_per_sheet(pl_variants)
  for group_number, sample_group in enumerate(samples_per_sheet):
    pl_one_group = pl_variants.filter(pl.col(SAMPLE_COL).is_in(sample_group))
    if len(samples_per_sheet) == 1:
      one_sheet_name = sheet_basename
    else:
      one_sheet_name = f"{sheet_basename} {group_number}"
    if len(sample_group) == 1 and sample_group[0] in samples_to_split:
      max_one_group_sheets = math.ceil(pl_one_group.shape[0] / MAX_ROWS_EXCEL)
      for one_group_sheet in range(max_one_group_sheets):
        # append a, b, c, ..., z
        subsheet_name = one_sheet_name + f"_{string.ascii_lowercase[one_group_sheet]}"
        pl_group_slice = pl_one_group.slice(one_group_sheet * MAX_ROWS_EXCEL, (one_group_sheet + 1) * MAX_ROWS_EXCEL)
        pl_group_slice = pl_group_slice.sort([SAMPLE_COL, CHROM_COL, POSITION_COL])
        write_df_to_excel(pl_group_slice, workbook.add_worksheet(subsheet_name), pl_variants_widths, header_format)
    else:
      pl_one_group = pl_one_group.sort([SAMPLE_COL, CHROM_COL, POSITION_COL])
      write_df_to_excel(pl_one_group, workbook.add_worksheet(one_sheet_name), pl_variants_widths, header_format)


def write_summary_sheet(workbook: xlsxwriter.Workbook, pl_summary: pl.DataFrame, sheet_name: str, header_format):
  pl_summary_widths = get_column_widths(pl_summary)
  not_sample_name_cols = pl_summary_widths.keys() - set([SAMPLE_COL, f"{SAMPLE_COL}{REPEAT_COL_SUFFIX}"])
  not_sample_name_cols = list(not_sample_name_cols)
  for one_col in not_sample_name_cols:
    pl_summary_widths[one_col] = 100.0
  write_df_to_excel(pl_summary, workbook.add_worksheet(sheet_name), pl_summary_widths, header_format)


def write_lego_plot_sheet(workbook: xlsxwriter.Workbook, lego_plot: io.BytesIO, sheet_name: str):
  lego_plot_sheet = workbook.add_worksheet(sheet_name)
  if lego_plot.getbuffer().nbytes > 0:
    lego_plot_sheet.insert_image("B2", "lego_plot_graph.png", {"image_data": lego_plot})
  else:
    lego_plot_sheet.write("B2", "No lego plot provided.")


def write_legend_sheet(workbook: xlsxwriter.Workbook, pl_legend: pl.DataFrame, sheet_name: str, header_format):
  pl_legend_widths = get_column_widths(pl_legend)
  write_df_to_excel(pl_legend, workbook.add_worksheet(sheet_name),
                    pl_legend_widths, header_format, user_autofilter=False)


def write_mouse_contamination_table(workbook: xlsxwriter.Workbook, pl_mouse_contamination: pl.DataFrame,
                                    sheet_name: str, header_format, percent_format):
  widths = get_column_widths(pl_mouse_contamination)
  write_df_to_excel(pl_mouse_contamination, workbook.add_worksheet(sheet_name),
                    widths, header_format, user_autofilter=True,
                    column_formats={CONTAMINATION_PCT_COL: percent_format})  


def write_all_hamas(workbook: xlsxwriter.Workbook, pl_all_hamas: pl.DataFrame, sheet_name: str,
                    header_format):
  widths = get_column_widths(pl_all_hamas)
  write_df_to_excel(pl_all_hamas, workbook.add_worksheet(sheet_name),
                    widths, header_format, user_autofilter=True)


def write_xlsx(
    file_name: str,
    all_variants: nw.DataFrame,
    variant_summary_table: nw.DataFrame,
    mouse_contamination_table: nw.DataFrame,
    all_hamas: nw.DataFrame,
    legend: nw.DataFrame,
    lego_plot: io.BytesIO,
):
  """
  :param nw.DataFrame all_variants: Asumes that all_variants are ordered by sample and chromosome position.
  :param lego_plot io.BytesIO: lego_plot.getbuffer().nbytes == 0 when there is no plot.
  """
  # Narwhals does not have a excel engine so this has to be done in polars
  pl_variants = all_variants.to_native()
  pl_summary = variant_summary_table.to_native()
  pl_legend = legend.to_native()
  pl_mouse_contamination = mouse_contamination_table.to_native()
  pl_all_hamas = all_hamas.to_native()
  with xlsxwriter.Workbook(file_name, options={"constant_memory": True, "use_zip64": True}) as workbook:
    header_format = workbook.add_format(HEADER_STYLE)
    percent_format = workbook.add_format(PERCENT_FORMAT)
    write_all_snv_indel_sheets(workbook, pl_variants, SNV_INDEL_SHEET, header_format)
    write_summary_sheet(workbook, pl_summary, VARIANT_SUMMARY_SHEET, header_format)
    if not pl_mouse_contamination.is_empty():
      write_mouse_contamination_table(workbook, pl_mouse_contamination, "Mouse Contamination", header_format,
                                      percent_format)
    if not pl_all_hamas.is_empty():
      write_all_hamas(workbook, pl_all_hamas, "Mouse Alleles (HAMAs)", header_format)
    write_lego_plot_sheet(workbook, lego_plot, "Mutation Lego Plot")
    write_legend_sheet(workbook, pl_legend, LEGEND_SHEET, header_format)
