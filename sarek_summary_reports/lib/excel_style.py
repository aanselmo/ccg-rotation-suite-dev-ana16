#!/opt/conda/bin/python
"""
Some helper functions to add style to worksheets in openpyxl.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel, Julius Nevin"

import polars as pl
# Narwhals has no excel writing engine so it's excluded from imports


# Create 'header' format to style page header
HEADER_STYLE = {
    "bold": True,
    "text_wrap": True,
    "align": "left",
}


def get_column_widths(df: pl.DataFrame, font_multiplier=10.0, first_row_multiplier=0) -> dict:
  """
  Get the column withs of a worksheet based on how many letters there are in a cell.
  :param font_multiplier float: a multiplier to increase if the letters are bigger. Default is 1.8
  :param first_row_multiplier float: Use a different multiplier for the first row (i.e. the header). Zero means use the
      same as the rest of the data.
  :return dict: a dictionary of column names to their column widths.
  """
  column_to_width = {}
  max_width = df.select(pl.all().cast(pl.Utf8).str.len_chars().max() * font_multiplier)
  if not first_row_multiplier:
    first_row_multiplier = font_multiplier
  name_width = dict(map(lambda column_name: (column_name, len(column_name) * first_row_multiplier), df.columns))
  max_width.extend(pl.DataFrame(name_width))
  max_width = max_width.max()
  column_to_width = max_width.row(0, named=True)
  return column_to_width
