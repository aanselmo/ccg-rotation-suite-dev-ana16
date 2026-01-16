
#!/usr/local/miniforge3/bin/python

"""
Used to convert a .vcf.gz file to a narwhals dataframe.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import gzip
import polars as pl
import narwhals.stable.v2 as nw


INFO_COL = "INFO"
CSQ_COL = "CSQ"
FORMAT_COL = "FORMAT"


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  parser.add_argument("-o", "--output-file", type=str, default="output.csv",
                      dest="output_file_name", help="The name of the output csv file.")
  parser.add_argument("-v", "--vcf", dest="vcf", type=str, required=True,
                      help=f"A vcf.gz file with the actual variant calls in them.")
  return parser.parse_args()


def find_header_line(vcf_path: str) -> int:
  header_line = 0
  with gzip.open(vcf_path, "rt") as vcf_file:
    for line_number, line in enumerate(vcf_file):
      if line.startswith("#CHROM"):
        header_line = line_number
        break
  return header_line


def get_csq_header(vcf_path: str) -> list[str]:
  """
  :return list: Will return empty list if no CSQ header is found.
  """
  header = []
  with gzip.open(vcf_path, "rt") as vcf_file:
    for line in vcf_file:
      if line.startswith("##INFO=<ID=CSQ"):
        header = line.strip().split("|")
        last_whitespace = header[0].rfind(" ")
        header[0] = header[0][last_whitespace + 1:]
        header[-1] = header[-1].rstrip('">')
  if header:
    af_index = header.index("AF")
    header[af_index] = "gnomAD_AF"
  return header


def get_format_data_column_name(df: pl.DataFrame):
  format_data_col_index = df.get_column_index(FORMAT_COL) + 1
  return df.columns[format_data_col_index]


def drop_null_columns(df: pl.DataFrame) -> pl.DataFrame:
  """
  Drop columns that are completely null.
  :param df pl.DataFrame: A dataframe with potentially some null columns.
  :return pl.DataFrame: A dataframe with null columns dropped.
  """
  columns_to_drop = []
  for column_name in df.columns:
    if df.get_column(column_name).null_count() == df.height:
      columns_to_drop.append(column_name)
  df = df.drop(columns_to_drop)
  return df


def parse_info_section_headers(df: pl.DataFrame) -> pl.DataFrame:
  """
  The previous str.split.list.to_struct operation leaves columns like:
  |field_0              | field_1                    |
  ====================================================
  |AS_FilterStatus=SITE | AS_SB_TABLE=264,236|99,106 |

  This parses out the AS_FilterStatus and assigns it to the column name properly.
  :param df pl.DataFrame: DataFrame to be parsed.
  :return pl.DataFrame: DataFrame formatted properly.
  """
  columns_to_format = [column for column in df.columns if column.startswith("field")]
  # Get unique column names
  data_split_as_fields = df.select(columns_to_format).cast(pl.Utf8)
  unique_column_names = data_split_as_fields.select(
      pl.all().str.split("=").list.get(0)
  ).unique().to_numpy().reshape(-1).tolist()
  unique_column_names = list(set(unique_column_names))
  unique_column_names = list(filter(lambda x: x is not None, unique_column_names))
  corrected_info_section = pl.DataFrame({})
  # transpose for the filter operation (polars ONLY filters rows)
  data_split_as_fields = data_split_as_fields.transpose()
  for new_column_name in unique_column_names:
    print(new_column_name)
    # filter out starts with column name
    rows_with_new_column_name = (
        (pl.all() == new_column_name)
        | (pl.all().str.starts_with(f"{new_column_name}="))
    )
    # Filter any rows with the column name
    new_column = data_split_as_fields.filter(pl.any_horizontal(rows_with_new_column_name))
    # Filter any cells with the column name
    new_column = new_column.select(pl.when(rows_with_new_column_name).then(pl.all()).otherwise(pl.lit("")))
    # Verify that there is only 1 real value per column
    num_values = new_column.select(rows_with_new_column_name).sum()
    num_values = num_values.transpose()
    if not num_values.filter(pl.nth(0) > 1).is_empty():
      raise IndexError(f"There are multiple values in a row to flatten into 1 column ({new_column_name})")
    # transpose the filtered rows again
    new_column = new_column.transpose()
    # concat_str to flatten into 1 column
    new_column = new_column.select(pl.concat_str(pl.all()).alias(new_column_name))
    # Remove repeated column name in all rows.
    if "=" in str(new_column.item(0, 0)):
      new_column = new_column.with_columns(pl.col(new_column_name).str.split("=").list.get(1).alias(new_column_name))
    corrected_info_section = pl.concat([corrected_info_section, new_column], how="horizontal")
  df = df.drop(columns_to_format)
  df = pl.concat([df, corrected_info_section], how="horizontal")
  df = df.rename({"DP": "INFO_DP"})
  return df


def split_format_column(df: pl.DataFrame) -> pl.DataFrame:
  """
  Take a Dataframe with the format column (i.e. header column names) and it's data column and split the data into
  it's own column. The format column should be first, followed by the data column.
  :param df pl.DataFrame: A dataframe with 2 columns only.
  :return pl.DataFrame: A dataframe with those 2 columns split.
  """
  new_df = pl.DataFrame({})
  for one_row in df.rows():
    FIRST_FORMAT_COL = 0
    SECOND_DATA_COL = 1
    header = one_row[FIRST_FORMAT_COL].split(":")
    data = one_row[SECOND_DATA_COL].split(":")
    combined = pl.DataFrame(dict(zip(header, data)))
    new_df = pl.concat([new_df, combined], how="diagonal_relaxed")
  return new_df


def convert_vcf_to_data_frame(vcf_path: str) -> nw.DataFrame:
  header_line = find_header_line(vcf_path)
  csq_header = get_csq_header(vcf_path)
  df = pl.read_csv(
      vcf_path,
      skip_lines=header_line,
      separator="\t",
      schema_overrides={"#CHROM": pl.String},
  )
  # Data Col of Fromat Header Starts off as the same name so it's not constant.
  data_col_of_format_header = get_format_data_column_name(df)
  format_data = split_format_column(df.select(pl.col(FORMAT_COL), pl.col(data_col_of_format_header)))
  df = df.with_columns(pl.col(INFO_COL).str.split(";").list.to_struct(upper_bound=20).struct.unnest())
  # The to_struct() expression creates null columns.
  df = drop_null_columns(df)
  df = df.drop(INFO_COL)
  df = parse_info_section_headers(df)
  if csq_header:
    df = df.with_columns(pl.col(CSQ_COL).str.split("|").list.to_struct(fields=csq_header).struct.unnest())
    df = df.drop(CSQ_COL)
  df = df.drop([FORMAT_COL, data_col_of_format_header])
  df = pl.concat([df, format_data], how="horizontal")
  return nw.from_native(df)


def main():
  args = setup_arguments()
  df = convert_vcf_to_data_frame(args.vcf)
  df.write_csv(args.output_file_name)


if __name__ == "__main__":
  main()
