
#!/usr/local/miniforge3/bin/python

"""
Used to convert a .vcf.gz file to a narwhals dataframe.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import gzip
from os.path import dirname
from pathlib import Path
import sys
import polars as pl
import narwhals.stable.v2 as nw
script_dir = Path(__file__).resolve().parent
sarek_script_top_dir = dirname(dirname(script_dir))
sys.path.append(sarek_script_top_dir)
from lib.common import AVAILABLE_CHROMOSOMES

INFO_COL = "INFO"
CSQ_COL = "CSQ"
FORMAT_COL = "FORMAT"
VARIANT_CLASS_COL = "VARIANT_CLASS"
CONSEQUENCE_COL = "Consequence"
AMINO_ACIDS_COL = "Amino_acids"
EXISTING_VARIATION_COL = "Existing_variation"
MAX_AF_POPS_COL = "MAX_AF_POPS"
BEST_EFFECT_COLS_TO_USE = [
    CONSEQUENCE_COL, VARIANT_CLASS_COL, "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "Gene",
    "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "cDNA_position", "CDS_position",
    "Protein_position", AMINO_ACIDS_COL
]
CANONICAL_COLS_TO_USE = BEST_EFFECT_COLS_TO_USE
GNOMAD_1000_GENOMES_AF_COL = "gnomAD_AF"
GNOMAD_COLS = [GNOMAD_1000_GENOMES_AF_COL, "gnomADe_AF", "gnomADg_AF", "MAX_AF", MAX_AF_POPS_COL]
COSMIC_COLS = [
    "COSMIC_v104", "COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT", "COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT",
    "COSMIC_v80", "COSMIC_v80_CNT"
]
CSQ_COLS_TO_USE = [EXISTING_VARIATION_COL] + COSMIC_COLS + GNOMAD_COLS


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
    header[af_index] = GNOMAD_1000_GENOMES_AF_COL
  return header


def get_csq_columns(csq_header: list[str], df: pl.DataFrame):
  csq = df.select(pl.col(CSQ_COL).str.split(",").list.eval(
      pl.element().str.split("|").list.to_struct(fields=csq_header)
  ))
  # best_effect = csq.select(pl.col(CSQ_COL).list.get(0).struct.unnest())
  best_effect = csq.select(
      pl.col(CSQ_COL).list.filter(pl.element().struct.field("PICK") == "1")
      .list.first()
      .struct.unnest()
  )
  canonical = csq.select(
      pl.col(CSQ_COL).list.filter(pl.element().struct.field("CANONICAL") == "YES")
      .list.first()
      .struct.unnest()
  )
  breakpoint()
  csq = best_effect.select(CSQ_COLS_TO_USE)
  csq = csq.with_columns(
      pl.col(EXISTING_VARIATION_COL).str.split("&")
      .list.filter(~pl.element().str.contains("COSV"))
      .list.join("&").alias(EXISTING_VARIATION_COL)
  )
  best_effect = best_effect.select(pl.col(BEST_EFFECT_COLS_TO_USE).name.prefix("BestEffect_"))
  canonical = canonical.select(pl.col(CANONICAL_COLS_TO_USE).name.prefix("Canonical_"))
  csq = pl.concat([canonical, best_effect, csq], how="horizontal")
  return csq


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


def parse_info_section_headers(all_variants: pl.DataFrame) -> pl.DataFrame:
  """
  The previous str.split.list.to_struct operation leaves columns like:
  |field_0              | field_1                    |
  ====================================================
  |AS_FilterStatus=SITE | AS_SB_TABLE=264,236|99,106 |

  This parses out the AS_FilterStatus and assigns it to the column name properly.
  :param df pl.DataFrame: DataFrame to be parsed.
  :return pl.DataFrame: DataFrame formatted properly.
  """
  columns_to_format = [column for column in all_variants.columns if column.startswith("field")]
  # Get unique column names
  info_sec_as_fields = all_variants.select(columns_to_format).cast(pl.Utf8)
  info_sec_as_fields = info_sec_as_fields.select(pl.all().str.split("="))
  info_sec_headers = info_sec_as_fields.select(pl.all().list.get(0))
  info_sec_data = pl.DataFrame({})
  for column in info_sec_as_fields.columns:
    one_column = info_sec_as_fields.get_column(column).map_elements(
        lambda cell_list: cell_list[1] if (len(cell_list) > 1) else cell_list[0],
        return_dtype=pl.String
    )
    info_sec_data = pl.concat([info_sec_data, one_column.to_frame()], how="horizontal")
  corrected_info_section = []
  for row_index in range(info_sec_as_fields.shape[0]):
    headers = info_sec_headers.row(row_index)
    values = info_sec_data.row(row_index)
    head_to_val_map = dict(zip(headers, values))
    head_to_val_map = {k: v for k, v in head_to_val_map.items() if k is not None}
    corrected_info_section.append(head_to_val_map)
  all_variants = all_variants.drop(columns_to_format)
  df_corrected_info_section = pl.DataFrame(corrected_info_section)
  all_variants = pl.concat([all_variants, df_corrected_info_section], how="horizontal")
  all_variants = all_variants.rename({"DP": "INFO_DP"})
  return all_variants


def split_format_column(format_and_data_cols: pl.DataFrame) -> pl.DataFrame:
  """
  Take a Dataframe with the format column (i.e. header column names) and its data column and split the data into
  its own column. The format column should be first, followed by the data column.
  :param df pl.DataFrame: A dataframe with 2 columns only.
  :return pl.DataFrame: A dataframe with those 2 columns split.
  """
  FIRST_FORMAT_COL = 0
  SECOND_DATA_COL = 1
  split_column = format_and_data_cols.select(
      pl.nth(FIRST_FORMAT_COL).str.split(":"),
      pl.nth(SECOND_DATA_COL).str.split(":")
  )
  split_column = [dict(zip(row[FIRST_FORMAT_COL], row[SECOND_DATA_COL])) for row in split_column.iter_rows()]
  split_column = pl.DataFrame(split_column)
  return split_column


def convert_vcf_to_data_frame(vcf_path: str) -> nw.DataFrame:
  print(vcf_path)
  header_line = find_header_line(vcf_path)
  csq_header = get_csq_header(vcf_path)
  df = pl.read_csv(
      vcf_path,
      skip_lines=header_line,
      separator="\t",
      schema_overrides={"#CHROM": AVAILABLE_CHROMOSOMES},
  )
  # Data Col of Fromat Header Starts off as the same name so its not constant.
  data_col_of_format_header = get_format_data_column_name(df)
  format_data = split_format_column(df.select(pl.col(FORMAT_COL), pl.col(data_col_of_format_header)))
  df = df.with_columns(pl.col(INFO_COL).str.split(";").list.to_struct(upper_bound=20).struct.unnest())
  # The to_struct() expression creates null columns.
  df = drop_null_columns(df)
  df = df.drop(INFO_COL)
  df = parse_info_section_headers(df)
  if csq_header:
    csq_df = get_csq_columns(csq_header, df)
  df = df.drop([FORMAT_COL, data_col_of_format_header])
  df = pl.concat([df, format_data, csq_df], how="horizontal")
  return nw.from_native(df)


def main():
  args = setup_arguments()
  df = convert_vcf_to_data_frame(args.vcf)
  df.write_csv(args.output_file_name)


if __name__ == "__main__":
  main()
