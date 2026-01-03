
#!/usr/local/miniforge3/bin/python

"""
Used to convert a .vcf.gz file to a narwhals dataframe.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import gzip
import polars as pl
import narwhals.stable.v2 as nw


INFO_COL = "INFO"
CSQ_COL = "CSQ"
FORMAT_COL = "FORMAT"


def find_header_line(vcf_path: str) -> int:
  header_line = 0
  with gzip.open(vcf_path, "rt") as vcf_file:
    for line_number, line in enumerate(vcf_file):
      if line.startswith("#CHROM"):
        header_line = line_number
        break
  return header_line


def get_csq_header(vcf_path: str) -> list[str]:
  header = []
  with gzip.open(vcf_path, "rt") as vcf_file:
    for line in vcf_file:
      if line.startswith("##INFO=<ID=CSQ"):
        header = line.strip().split("|")
        last_whitespace = header[0].rfind(" ")
        header[0] = header[0][last_whitespace + 1:]
        header[-1] = header[-1].rstrip('">')
  af_index = header.index("AF")
  header[af_index] = "CSQ_AF"
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
  columns_to_format = df.columns
  columns_to_format = filter(lambda column: column.startswith("field"), columns_to_format)
  for column_name in columns_to_format:
    one_cell: str = df.item(0, column_name)
    new_column_name = one_cell.split("=")[0]
    df = df.rename({column_name: new_column_name})
    df = df.with_columns(pl.col(new_column_name).str.split("=").list.get(1).alias(new_column_name))
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
  df = df.with_columns(pl.col(CSQ_COL).str.split("|").list.to_struct(fields=csq_header).struct.unnest())
  df = df.drop(CSQ_COL)
  df = df.drop([FORMAT_COL, data_col_of_format_header])
  df = pl.concat([df, format_data], how="horizontal")
  return nw.from_native(df)


def main():
  df = convert_vcf_to_data_frame(
      "/home/nap24/Documents/Pipeline/bulk_dev/Sarek/20251218_CCG1200c_ponv2/annotation/mutect2/A008_X108_58797-1-C1/A008_X108_58797-1-C1.mutect2.filtered.bcftools_filtered_VEP.ann.vcf.gz")
  df.write_csv("/home/nap24/bin/ccgd-rotation-suite/scratch/test_vcf_converter.csv")


if __name__ == "__main__":
  main()