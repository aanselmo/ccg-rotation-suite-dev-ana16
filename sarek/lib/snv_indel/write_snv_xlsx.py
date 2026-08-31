#!/opt/conda/bin/python
"""
Write the SNV Indel data to an excel file.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"

import io
import polars as pl
import narwhals as nw
import xlsxwriter
from .vcf2df import (
    BEST_EFFECT_AMINO_ACIDS_COL, BEST_EFFECT_CONSEQUENCE_COL, EXISTING_VARIATION_COL, MAX_AF_POPS_COL
)
from .snv_indel_sheet import (
    SAMPLE_COL, REF_COL, ALT_COL, RECOMMENDED_FILTER_COL, MUTECT2_FILTER_COL
)
from .variant_summary_sheet import REPEAT_COL_SUFFIX
from ..excel_style import get_column_widths, HEADER_STYLE


MAX_ROWS_EXCEL = 1048576


def get_samples_per_sheet(pl_variants: pl.DataFrame):
  rows_per_sample = pl_variants.group_by(pl.col(SAMPLE_COL)).len()
  samples_per_sheet = [[]]
  rows_in_sheet = 1  # start with 1 for header
  sheet_number = 0
  for row in rows_per_sample.iter_rows(named=True):
    rows_in_sheet = rows_in_sheet + row["len"]
    if rows_in_sheet <= MAX_ROWS_EXCEL:
      samples_per_sheet[sheet_number] = row[SAMPLE_COL]
    else:
      sheet_number = sheet_number + 1
      rows_in_sheet = row["len"] + 1
      samples_per_sheet[sheet_number] = row[SAMPLE_COL]
  return samples_per_sheet


def write_xlsx(
    file_name: str,
    all_variants: nw.DataFrame,
    variant_summary_table: nw.DataFrame,
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
  pl_variants_widths = get_column_widths(pl_variants, first_row_multiplier=13.0)
  for one_col in [
      REF_COL, ALT_COL, RECOMMENDED_FILTER_COL, MUTECT2_FILTER_COL, BEST_EFFECT_CONSEQUENCE_COL,
      BEST_EFFECT_AMINO_ACIDS_COL, EXISTING_VARIATION_COL, MAX_AF_POPS_COL
  ]:
    pl_variants_widths[one_col] = 100.0
  pl_summary_widths = get_column_widths(pl_summary)
  not_sample_name_cols = pl_summary_widths.keys() - set([SAMPLE_COL, f"{SAMPLE_COL}{REPEAT_COL_SUFFIX}"])
  not_sample_name_cols = list(not_sample_name_cols)
  for one_col in not_sample_name_cols:
    pl_summary_widths[one_col] = 100.0
  with xlsxwriter.Workbook(file_name) as workbook:
    samples_per_sheet = get_samples_per_sheet(pl_variants)
    if len(samples_per_sheet) == 1:
      pl_variants.write_excel(
          workbook=workbook,
          worksheet="SNVs and Indels",
          column_widths=pl_variants_widths,
          header_format=HEADER_STYLE,
          dtype_formats={pl.Int64: "0;[RED]-0"},
      )
    else:
      for group_number, sample_group in enumerate(samples_per_sheet):
        pl_one_group = pl_variants.filter(pl.col(SAMPLE_COL).is_in(sample_group))
        pl_one_group.write_excel(
            workbook=workbook,
            worksheet=f"SNVs and Indels {group_number}",
            column_widths=pl_variants_widths,
            header_format=HEADER_STYLE,
            dtype_formats={pl.Int64: "0;[RED]-0"},
        )
    pl_summary.write_excel(
        workbook=workbook,
        worksheet="Variant Summary Table",
        column_widths=pl_summary_widths,
        header_format=HEADER_STYLE,
        dtype_formats={pl.Int64: "0;[RED]-0"},
    )
    lego_plot_sheet = workbook.add_worksheet("Mutation Lego Plot")
    if lego_plot.getbuffer().nbytes > 0:
      lego_plot_sheet.insert_image("B2", "lego_plot_graph.png", {"image_data": lego_plot})
    else:
      lego_plot_sheet.write("B2", "No lego plot provided.")

    pl_legend.write_excel(
        workbook=workbook,
        worksheet="Legend",
        column_widths=get_column_widths(pl_legend),
        header_format=HEADER_STYLE,
        dtype_formats={pl.Int64: "0;[RED]-0"},
    )
