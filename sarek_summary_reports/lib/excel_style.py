#!/opt/conda/bin/python
"""
Some helper functions to add style to worksheets in openpyxl.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel, Julius Nevin"

import xlsxwriter

def as_text(value):
  if value is None:
    return ""
  return str(value)


def adjust_column_widths(worksheet, font_multiplier=1.8, first_row_multiplier=0):
  """
  Adjust the column with of a worksheet based on how many letters there are in a cell.
  :param font_multiplier float: a multiplier to increase if the letters are bigger. Default is 1.8
  :param first_row_multiplier float: Use a different multiplier for the first row. Zero means use the same as
      the rest. 
  """
  for column_cells in worksheet.columns:
    length = 0
    for row_num, cell in enumerate(column_cells):
      if row_num == 0 and first_row_multiplier:
        length = len(as_text(cell.value)) * first_row_multiplier
      else:
        length = max(len(as_text(cell.value)) * font_multiplier, length)
    worksheet.column_dimensions[column_cells[0].column_letter].width = length
  return worksheet


def get_range_of_column(column_cells, header_row=0):
  low_range = 100000000
  high_range = -100000000
  for row_num, cell in enumerate(column_cells):
    if row_num > header_row:
      if cell.value > high_range:
        high_range = cell.value
      elif cell.value < low_range:
        low_range = cell.value
  return (low_range, high_range)


def convert_to_numerical(df_in: pd.DataFrame) -> pd.DataFrame:
  df_out = df_in.copy()
  for col_index, column_name in enumerate(df_in.columns):
    if pd.api.types.is_string_dtype(df_in.iloc[:, col_index]):
      if re.fullmatch(r"^\d+$", df_in.iat[0, col_index]):
        df_out.iloc[:, col_index] = df_in.iloc[:, col_index].astype("int64")
      elif re.fullmatch(r"^\d+\.\d+$", df_in.iat[0, col_index]):
        df_out.iloc[:, col_index] = df_in.iloc[:, col_index].astype("float64")
      elif re.fullmatch(r"^.*%$", df_in.iat[0, col_index]):
        df_out.iloc[:, col_index] = df_in.iloc[:, col_index].str.replace("%", "").astype("float64") / 100
  return df_out


def format_numerical_data(worksheet, exclude_columns=[], format_map={}, header_row=0):
  """
  :param header_row int: The row where the column names are.
  :param exclude_columns list[str]: A list of column names in the header_row that you want to exclude from this
      format_numerical_data formatting.
  :param format_map dict{ str : str }: A map of the column name to the format you want to use. This is if you have a
      desired numerical format for a column that you for sure know you want to use.
  """
  for column_cells in worksheet.columns:
    header_val = column_cells[header_row].value
    if header_val in format_map.keys():
      mapped_column_format = format_map[header_val]
      for cell in column_cells[header_row+1:]:
        cell.number_format = mapped_column_format
    elif header_val not in exclude_columns:
      if isinstance(column_cells[header_row+1].value, float) or isinstance(column_cells[header_row+1].value, int):
        (low_range, high_range) = get_range_of_column(column_cells, header_row)
        if low_range == 0 and high_range == 0:
          for cell in column_cells[header_row+1:]:
            cell.number_format = "#,###0"
        elif low_range >= 0 and high_range <= 1:
          for cell in column_cells[header_row+1:]:
            cell.number_format = numbers.FORMAT_PERCENTAGE_00
        elif high_range >= 1000:
          for cell in column_cells[header_row+1:]:
            if isinstance(cell.value, float):
              cell.number_format = numbers.FORMAT_NUMBER_COMMA_SEPARATED1
            else:
              cell.number_format = "#,###0"


# Create 'header' format to style page header
HEADER_STYLE = NamedStyle(name="header")
HEADER_STYLE.font = Font(name='Arial', bold=True, underline='single', size="12")
HEADER_STYLE.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

INDEX_COL_STYLE = NamedStyle(name="index_col")
INDEX_COL_STYLE.font = Font(name='Arial', size="12")
INDEX_COL_STYLE.alignment = Alignment(horizontal="left", vertical="bottom")

# Create 'data' format to style tool statistic columns
DATA_STYLE = NamedStyle(name="data")
DATA_STYLE.font = Font(name='Arial', size="12")
DATA_STYLE.alignment = Alignment(horizontal="right", vertical="bottom")


def style_worksheets(output_file_name: str, sheets_to_style: list = [-1]):
  """
  :param sheets_to_style list: A list of sheet names to style. By default if [-1] is passed, it styles all of them.
  """
  workbook = load_workbook(filename=output_file_name, data_only=False)
  if sheets_to_style == [-1]:
    sheets_to_style = workbook.sheetnames
  for sheet_name in sheets_to_style:
    sheet = workbook[sheet_name]
    for cell in sheet[1]:
      cell.style = HEADER_STYLE
    sheet.row_dimensions[1].height = 45
    for row_num in range(2, sheet.max_row + 1):
      sheet.cell(row=row_num, column=1).style = INDEX_COL_STYLE
      for col_num in range(2, sheet.max_column + 1):
        sheet.cell(row=row_num, column=col_num).style = DATA_STYLE
    adjust_column_widths(sheet, font_multiplier=1.7, first_row_multiplier=0.58)
    format_numerical_data(sheet)
  return workbook
