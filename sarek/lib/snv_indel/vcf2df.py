
#!/usr/local/miniforge3/bin/python

"""
Used to convert a .vcf.gz file to a narwhals dataframe.

There is are some columns filtered/excluded in this module, but for the most part I try to preserve everything for
later filtration. I intended this module to just convert the formats. The data that are excluded are the following:
    1. Any CSQ (AKA Consequences or AKA Transcripts) columns besides Best Effect and Canonical.
    2. Old COSMIC v98 annotations from the VEP cache in Existing_Variation column.
    3. Only use the variants that PASS the Mutect2 FILTER column.
    4. Allele fractions that are <2%

Things this module does do:
    1. splitting columns into simpler, unique columns.

Things this module does NOT do (to keep things organized - put these in snv_indel_sheet.py):
    1. Add extra filter columns.
    2. Combining columns.
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

CHROM_COL = "#CHROM"
POSITION_COL = "POS"
REF_COL = "REF"
ALT_COL = "ALT"
INFO_COL = "INFO"
CSQ_COL = "CSQ"
FORMAT_COL = "FORMAT"
AF_COL = "AF"
DP_COL = "DP"
MUTECT2_FILTER_COL = "FILTER"
BEST_EFFECT_PREFIX = "BestEffect_"
CANONICAL_PREFIX = "Canonical_"
CSQ_VARIANT_CLASS_COL = "VARIANT_CLASS"
CSQ_CONSEQUENCE_COL = "Consequence"
EXISTING_VARIATION_COL = "Existing_variation"
CSQ_BIOTYPE_COL = "BIOTYPE"
BEST_EFFECT_CONSEQUENCE_COL = BEST_EFFECT_PREFIX + CSQ_CONSEQUENCE_COL
CSQ_HGVS_C = "HGVSc"
CSQ_HGVS_P = "HGVSp"
BEST_EFFECT_HGVS_C_COL = BEST_EFFECT_PREFIX + CSQ_HGVS_C
BEST_EFFECT_HGVS_P_COL = BEST_EFFECT_PREFIX + CSQ_HGVS_P
BEST_EFFECT_COLS_TO_USE = [
    CSQ_CONSEQUENCE_COL, "EXON", "INTRON", "cDNA_position", CSQ_HGVS_C, CSQ_HGVS_P
]
CANONICAL_COLS_TO_USE = BEST_EFFECT_COLS_TO_USE
CANONICAL_CONSEQUENCE_COL = CANONICAL_PREFIX + CSQ_CONSEQUENCE_COL
CANONICAL_HGVS_C_COL = CANONICAL_PREFIX + CSQ_HGVS_C
CANONICAL_HGVS_P_COL = CANONICAL_PREFIX + CSQ_HGVS_P
GNOMAD_1000_GENOMES_AF_COL = "gnomAD_AF"
GNOMAD_EXOMES_AF_COL = "gnomADe_AF"
GNOMAD_WHOLE_GENOMES_AF_COL = "gnomADg_AF"
MAX_AF_COL = "MAX_AF"
MAX_AF_POPS_COL = "MAX_AF_POPS"
GNOMAD_COLS = [GNOMAD_1000_GENOMES_AF_COL, "gnomADe_AF", "gnomADg_AF", "MAX_AF", MAX_AF_POPS_COL]
COSMIC_v104_COL = "COSMIC_v104"
COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT_COL = "COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT"
COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT_COL = "COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT"
COSMIC_v80_COL = "COSMIC_v80"
COSMIC_v80_CNT_COL = "COSMIC_v80_CNT"
COSMIC_COLS = [
    COSMIC_v104_COL, COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT_COL, COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT_COL,
    COSMIC_v80_COL, COSMIC_v80_CNT_COL
]
CSQ_VARIANT_TYPE_COLS = [CSQ_VARIANT_CLASS_COL, "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "Gene", CSQ_BIOTYPE_COL]
CSQ_POPULATION_COLS = [EXISTING_VARIATION_COL] + COSMIC_COLS + GNOMAD_COLS


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process 1 vcf to turn it into a csv. Useful for testing.")
  parser.add_argument("-o", "--output-file", type=str, default="output.csv",
                      dest="output_file_name", help="The name of the output csv file.")
  parser.add_argument("-v", "--vcf", dest="vcf", type=str, required=True,
                      help=f"A vcf.gz file with the actual variant calls in them.")
  parser.add_argument("-f", "--filter-mutect-fail", dest="filter_mutect_fail", action="store_true",
                      help=f"Keep only the PASS calls for Mutect2 FILTER column.")
  return parser.parse_args()


def parse_vcf_header(vcf_path: str) -> tuple[int, list[str], list[str]]:
  """
  Parse out information from the .vcf header (header section lines begin with a hashtag sign i.e. #)

  :param vcf_path str: The file path to open and parse the header from.
  :return tuple[int, list[str], list[str]]: (1) The int is the header_line where you can read the column names for the
      table. (2) The first list[str] is the col_names in the info section of the dataframe.
      (3) The second list[str] is the csq_column_names or the column names parsed out from the CSQ section of the .vcf
      that gets put into the DataFrame.
  """
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
      elif line.startswith(CHROM_COL):
        # This is the last piece of info needed to grabbed from the header so you can exit the loop after this.
        header_line = line_number
        break
  if csq_col_names:
    af_index = csq_col_names.index("AF")
    csq_col_names[af_index] = GNOMAD_1000_GENOMES_AF_COL
    # Apparently, CSQ columns are noted in the header the same way as INFO columns sometimes.
    info_col_names = list(set(info_col_names) - set(csq_col_names))
  return (header_line, info_col_names, csq_col_names)


def split_csq_column(csq_header: list[str], csq_column: pl.DataFrame):
  """
  Take the CSQ column and split it into many unique, smaller columns.
  :param csq_header list[str]: the column names that the CSQ section will be split into.
  :param csq_column pl.DataFrame: The original CSQ column to be split into smaller columns.
  :return csq pl.DataFrame: Many columns after split.
  """
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
  csq = best_effect.select(CSQ_VARIANT_TYPE_COLS + CSQ_POPULATION_COLS)
  # Removing the COSV annotations from the existing variation column because it's the old v98 ones
  # We have COSMIC_v104 columns in CSQ with updated information.
  csq = csq.with_columns(
      pl.col(EXISTING_VARIATION_COL).str.split("&")
      .list.filter(~pl.element().str.contains("COSV"))
      .list.join("&").alias(EXISTING_VARIATION_COL)
  )
  csq = csq.with_columns(
      pl.when(pl.col(COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT_COL).str.len_chars() == 0)
      .then(pl.lit(None))
      .otherwise(pl.col(COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT_COL))
      .alias(COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT_COL)
  )
  csq = csq.with_columns(
      pl.when(pl.col(COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT_COL).str.len_chars() == 0)
      .then(pl.lit(None))
      .otherwise(pl.col(COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT_COL))
      .alias(COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT_COL)
  )
  csq = csq.with_columns(
      pl.col(COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT_COL).str.to_integer(),
      pl.col(COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT_COL).str.to_integer()
  )
  best_effect = best_effect.select(pl.col(BEST_EFFECT_COLS_TO_USE).name.prefix(BEST_EFFECT_PREFIX))
  canonical = canonical.select(pl.col(CANONICAL_COLS_TO_USE).name.prefix(CANONICAL_PREFIX))
  csq = pl.concat([canonical, best_effect, csq], how="horizontal")
  return csq


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
  """
  This splits the single info_column into multiple columns with unique information.
  :param info_column pl.DataFrame: the original info_column to be split.
  :param column_names list[str]: The column names in the INFO section. 
  """
  info_df = info_column.with_columns([
      pl.col("INFO")
      .str.extract(fr"{one_col}=([^;]*)", 1)
      .alias(one_col)
      for one_col in column_names
  ])
  info_df = info_df.drop(INFO_COL)
  info_df = info_df.rename({"DP": "INFO_DP"})
  return info_df


def split_format_column(df: pl.DataFrame) -> pl.DataFrame:
  """
  Take a Dataframe with the FORMAT column (i.e. header column names) and its data column and split the data into
  its own column.
  e.g.
  1. If it's tumor-only mutect2, then it's right after the FORMAT column.
     FORMAT  01-16_TR00066-G9
     GT:AD:AF:DP:F1R2:F2R1:FAD:SB    0/1:0,1040:0.999:1040:0,370:0,362:0,850:0,0,530,510
  2. If it's matched tumor-normal mutect2, then it's the 2nd one after the FORMAT column.
     FORMAT  CHGFR0086_CHGFR0086-Germline-A4 CHGFR0086_CHGFR0086-T1-Org-5050-G3
     GT:AD:AF:DP:F1R2:F2R1:FAD:SB    0/0:795,3:5.085e-03:798:387,1:365,2:768,3:415,380,1,2   0/1:695,4:7.416e-03:699:346,1:295,3:660,4:324,371,1,3
  :param df pl.DataFrame: The partially converted DataFrame from the original vcf.
  :return pl.DataFrame: The original DataFrame with the 2 columns split.
  """
  FIRST_FORMAT_COL = 0
  SECOND_DATA_COL = 1
  index_format_col = df.get_column_index(FORMAT_COL)
  index_format_data_col = -1   # i.e. last column
  format_data_col = df.columns[index_format_data_col]
  index_format_data_col = df.get_column_index(format_data_col)
  if index_format_data_col != (index_format_col + 1):
    # then it's matched tumor-normal. Drop the normal data.
    df = df.drop(df.columns[index_format_col + 1])
  split_column = df.select(
      pl.col(FORMAT_COL).str.split(":"),
      pl.col(format_data_col).str.split(":")
  )
  split_column = [dict(zip(row[FIRST_FORMAT_COL], row[SECOND_DATA_COL])) for row in split_column.iter_rows()]
  split_column = pl.DataFrame(split_column, schema_overrides={AF_COL: pl.Float32, DP_COL: pl.Int16})
  df = df.drop([FORMAT_COL, format_data_col])
  df = pl.concat([df, split_column], how="horizontal")
  return df


def set_dtypes(df: pl.DataFrame) -> pl.DataFrame:
  dtypes = {
    POSITION_COL: pl.Int64,
    REF_COL: pl.String,
    ALT_COL: pl.String,
    MUTECT2_FILTER_COL: pl.String,
    CSQ_VARIANT_CLASS_COL: pl.String,
    "SYMBOL": pl.String,
    "SYMBOL_SOURCE": pl.String,
    "HGNC_ID": pl.String,
    "Gene": pl.String,
    CSQ_BIOTYPE_COL: pl.String,
    CANONICAL_CONSEQUENCE_COL: pl.String,
    "Canonical_EXON": pl.String,
    "Canonical_INTRON": pl.String,
    "Canonical_cDNA_position": pl.String,
    CANONICAL_HGVS_C_COL: pl.String,
    CANONICAL_HGVS_P_COL: pl.String,
    BEST_EFFECT_CONSEQUENCE_COL: pl.String,
    "BestEffect_EXON": pl.String,
    "BestEffect_INTRON": pl.String,
    "BestEffect_cDNA_position": pl.String,
    BEST_EFFECT_HGVS_C_COL: pl.String,
    BEST_EFFECT_HGVS_P_COL: pl.String,
    EXISTING_VARIATION_COL: pl.String,
    COSMIC_v104_COL: pl.String,
    COSMIC_v104_GENOME_SCREEN_SAMPLE_COUNT_COL: pl.Int64,
    COSMIC_v104_TARGETED_SCREEN_SAMPLE_COUNT_COL: pl.Int64,
    COSMIC_v80_COL: pl.String,
    COSMIC_v80_CNT_COL: pl.String,
    GNOMAD_1000_GENOMES_AF_COL: pl.Float64,
    GNOMAD_EXOMES_AF_COL: pl.Float64,
    GNOMAD_WHOLE_GENOMES_AF_COL: pl.Float64,
    MAX_AF_COL: pl.Float64,
    MAX_AF_POPS_COL: pl.String
  }
  if GNOMAD_1000_GENOMES_AF_COL in df.columns:
    df = df.with_columns(
        pl.col(GNOMAD_1000_GENOMES_AF_COL).replace([".", ""], None),
        pl.col(GNOMAD_EXOMES_AF_COL).replace([".", ""], None),
        pl.col(GNOMAD_WHOLE_GENOMES_AF_COL).replace([".", ""], None),
        pl.col(MAX_AF_COL).replace([".", ""], None)
    )
  for column, column_dtype in dtypes.items():
    if column in df.columns:
      df = df.cast({column: column_dtype})
  return df


def convert_vcf_to_data_frame(vcf_path: str, should_filter_mutect2_failed=True) -> nw.DataFrame:
  """
  Convert a .vcf file into a nw.DataFrame. Does not do any filtering for the most part. That's not part of this module.
  :param vcf_path str: the path to the file of .vcf that needs to be parsed.
  :return nw.DataFrame: A narwhals DataFrame that the .vcf converts to.
  """
  print(vcf_path)
  csq_df = pl.DataFrame({})
  (header_line, info_col_names, csq_col_names) = parse_vcf_header(vcf_path)
  df = pl.read_csv(
      vcf_path,
      skip_lines=header_line,
      separator="\t",
      schema_overrides={CHROM_COL: AVAILABLE_CHROMOSOMES},
  )
  if should_filter_mutect2_failed:
    # KEEP ONLY PASSING MUTECT2 VARIANTS
    df = df.filter(pl.col(MUTECT2_FILTER_COL) == "PASS")
  if not df.is_empty():
    df = split_format_column(df)
    df = df.filter(pl.col(AF_COL) >= 0.02)
    info_df = split_info_column(df.select(INFO_COL), info_col_names)
    if csq_col_names:
      csq_df = split_csq_column(csq_col_names, info_df.select(CSQ_COL))
      info_df = info_df.drop(CSQ_COL)
    df = df.drop([INFO_COL])
    df = pl.concat([df, info_df, csq_df], how="horizontal")
    df = set_dtypes(df)    
  return nw.from_native(df)


def main():
  """
  Main function here was written as a way to quickly test parsing of 1 file. Otherwise, because the snv_indel_report
  is parallelized/multi-process, the stack trace of the error can be obfuscated.
  """
  args = setup_arguments()
  df = convert_vcf_to_data_frame(args.vcf, should_filter_mutect2_failed=args.filter_mutect_fail)
  df.write_csv(args.output_file_name)


if __name__ == "__main__":
  main()
