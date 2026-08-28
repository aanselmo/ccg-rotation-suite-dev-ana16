
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
BEST_EFFECT_PREFIX = "BestEffect_"
CANONICAL_PREFIX = "Canonical_"
CSQ_VARIANT_CLASS_COL = "VARIANT_CLASS"
BEST_EFFECT_VARIANT_CLASS_COL = BEST_EFFECT_PREFIX + CSQ_VARIANT_CLASS_COL
CSQ_CONSEQUENCE_COL = "Consequence"
EXISTING_VARIATION_COL = "Existing_variation"
CSQ_BIOTYPE_COL = "BIOTYPE"
BEST_EFFECT_BIOTYPE_COL = BEST_EFFECT_PREFIX + CSQ_BIOTYPE_COL
BEST_EFFECT_CONSEQUENCE_COL = BEST_EFFECT_PREFIX + CSQ_CONSEQUENCE_COL
CSQ_AMINO_ACIDS = "Amino_acids"
BEST_EFFECT_AMINO_ACIDS_COL = BEST_EFFECT_PREFIX + CSQ_AMINO_ACIDS
BEST_EFFECT_COLS_TO_USE = [
    CSQ_CONSEQUENCE_COL, CSQ_VARIANT_CLASS_COL, "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "Gene",
    "Feature_type", "Feature", CSQ_BIOTYPE_COL, "EXON", "INTRON", "cDNA_position", "CDS_position",
    "Protein_position", CSQ_AMINO_ACIDS
]
CANONICAL_COLS_TO_USE = BEST_EFFECT_COLS_TO_USE
GNOMAD_1000_GENOMES_AF_COL = "gnomAD_AF"
MAX_AF_POPS_COL = "MAX_AF_POPS"
GNOMAD_COLS = [GNOMAD_1000_GENOMES_AF_COL, "gnomADe_AF", "gnomADg_AF", "MAX_AF", MAX_AF_POPS_COL]
COSMIC_v104_COL = "COSMIC_v104"
COSMIC_COLS = [
    COSMIC_v104_COL, "COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT", "COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT",
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


def parse_vcf_header(vcf_path: str) -> list[str]:
  info_col_names = []
  csq_col_names = []
  header_line = 0
  with gzip.open(vcf_path, "rt") as vcf_file:
    for line_number, line in enumerate(vcf_file):
      if line.startswith("##INFO=<ID="):
        split_line = line.split(",")
        one_column = split_line[0].replace("##INFO=<ID=", "")
        info_col_names.append(one_column)
        if line.startswith("##INFO=<ID=CSQ"):
          csq_col_names = line.strip().split("|")
          last_whitespace = csq_col_names[0].rfind(" ")
          csq_col_names[0] = csq_col_names[0][last_whitespace + 1:]
          csq_col_names[-1] = csq_col_names[-1].rstrip('">')
      elif line.startswith("#CHROM"):
        header_line = line_number
        break
  if csq_col_names:
    af_index = csq_col_names.index("AF")
    csq_col_names[af_index] = GNOMAD_1000_GENOMES_AF_COL
    # Apparently, CSQ columns are noted in the header the same way as INFO columns sometimes.
    info_col_names = list(set(info_col_names) - set(csq_col_names))
  return (header_line, info_col_names, csq_col_names)


def split_csq_column(csq_header: list[str], csq_column: pl.DataFrame):
  csq = csq_column.select(pl.col(CSQ_COL).str.split(",").list.eval(
      pl.element().str.split("|").list.to_struct(fields=csq_header)
  ))
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
  csq = best_effect.select(CSQ_COLS_TO_USE)
  csq = csq.with_columns(
      pl.col(EXISTING_VARIATION_COL).str.split("&")
      .list.filter(~pl.element().str.contains("COSV"))
      .list.join("&").alias(EXISTING_VARIATION_COL)
  )
  best_effect = best_effect.select(pl.col(BEST_EFFECT_COLS_TO_USE).name.prefix(BEST_EFFECT_PREFIX))
  canonical = canonical.select(pl.col(CANONICAL_COLS_TO_USE).name.prefix(CANONICAL_PREFIX))
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


def split_info_column(info_column: pl.DataFrame, column_names: list[str]):
  info_df = info_column.with_columns([
      pl.col("INFO")
      .str.extract(fr"{one_col}=([^;]*)", 1)
      .alias(one_col)
      for one_col in column_names
  ])
  info_df = info_df.drop(INFO_COL)
  info_df = info_df.rename({"DP": "INFO_DP"})
  return info_df


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
  (header_line, info_col_names, csq_col_names) = parse_vcf_header(vcf_path)
  df = pl.read_csv(
      vcf_path,
      skip_lines=header_line,
      separator="\t",
      schema_overrides={"#CHROM": AVAILABLE_CHROMOSOMES},
  )
  # Data Col of Fromat Header Starts off as the same name so its not constant.
  data_col_of_format_header = get_format_data_column_name(df)
  format_data = split_format_column(df.select(pl.col(FORMAT_COL), pl.col(data_col_of_format_header)))
  info_df = split_info_column(df.select(INFO_COL), info_col_names)
  if csq_col_names:
    csq_df = split_csq_column(csq_col_names, info_df.select(CSQ_COL))
  info_df = info_df.drop(CSQ_COL)
  df = df.drop([FORMAT_COL, data_col_of_format_header, INFO_COL])
  df = pl.concat([df, format_data, info_df, csq_df], how="horizontal")
  return nw.from_native(df)


def main():
  args = setup_arguments()
  df = convert_vcf_to_data_frame(args.vcf)
  df.write_csv(args.output_file_name)


if __name__ == "__main__":
  main()
