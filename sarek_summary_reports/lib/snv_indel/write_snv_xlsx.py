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
from .snv_indel_sheet import (
    SAMPLE_COL, REF_COL, ALT_COL, RECOMMENDED_FILTER_COL, CONSEQUENCE_COL, AMINO_ACIDS_COL, EXISTING_VARIATION_COL,
    MAX_AF_POPS_COL,
)
from .variant_summary_sheet import REPEAT_COL_SUFFIX
from ..excel_style import get_column_widths, HEADER_STYLE


def write_xlsx(
    file_name: str,
    all_variants: nw.DataFrame,
    variant_summary_table: nw.DataFrame,
    legend: nw.DataFrame,
    lego_plot: io.BytesIO,
):
  # Narwhals does not have a excel engine so this has to be done in polars
  pl_variants = all_variants.to_native()
  pl_summary = variant_summary_table.to_native()
  pl_legend = legend.to_native()
  pl_variants_widths = get_column_widths(pl_variants, first_row_multiplier=13.0)
  for one_col in [
      REF_COL, ALT_COL, RECOMMENDED_FILTER_COL, CONSEQUENCE_COL, AMINO_ACIDS_COL, EXISTING_VARIATION_COL,
      MAX_AF_POPS_COL
  ]:
    pl_variants_widths[one_col] = 100.0
  pl_summary_widths = get_column_widths(pl_summary)
  not_sample_name_cols = pl_summary_widths.keys() - set([SAMPLE_COL, f"{SAMPLE_COL}{REPEAT_COL_SUFFIX}"])
  not_sample_name_cols = list(not_sample_name_cols)
  for one_col in not_sample_name_cols:
    pl_summary_widths[one_col] = 100.0
  with xlsxwriter.Workbook(file_name) as workbook:
    pl_variants.write_excel(
        workbook=workbook,
        worksheet="SNVs and Indels",
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
    lego_plot_sheet.insert_image("B2", "lego_plot_graph.png", {"image_data": lego_plot})
    pl_legend.write_excel(
        workbook=workbook,
        worksheet="Legend",
        column_widths=get_column_widths(pl_legend),
        header_format=HEADER_STYLE,
        dtype_formats={pl.Int64: "0;[RED]-0"},
    )
