#!/usr/local/miniforge3/bin/python

"""
Gathers the data for the Summary Sheet by Gene in the CNV report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import sys
from os.path import dirname, basename as get_basename, dirname
from pathlib import Path
import polars as pl
import narwhals.stable.v2 as nw
script_dir = Path(__file__).resolve().parent
sarek_script_top_dir = dirname(dirname(script_dir))
sys.path.append(sarek_script_top_dir)
from lib.common import AVAILABLE_CHROMOSOMES
from .qc_sheet import SAMPLE_COL
from .variant_sheet import COPY_NUMBER_COL, CHROMOSOME_COL

GENE_COL = "gene"
LOG2_COL = "log2"
DIPLOID = 2.0
GAIN_COL = "No. of Samples with Gains"
LOSS_COL = "No. of Samples with Losses"
GAIN_AND_LOSS_COL = "No. of Samples with Gains & Losses"



def get_one_sample_variants(one_gene_file: str) -> nw.DataFrame:
    # Setup metadata
    base = get_basename(one_gene_file)
    sample_name = base.replace(".genemetrics.tsv", "")
    
    # Standard 12 columns from CNVkit genemetrics + our Sample column
    EXPECTED_COLUMNS = [
        GENE_COL, CHROMOSOME_COL, "start", "end", LOG2_COL, 
        "depth", "weight", "ci_hi", "ci_lo", "probes", 
        "segment_weight", "segment_probes", SAMPLE_COL
    ]

    try:
        # Attempt to read the file
        raw_df = pl.read_csv(
            one_gene_file,
            separator="\t",
            has_header=True,
            schema_overrides={ 
                CHROMOSOME_COL: AVAILABLE_CHROMOSOMES,
                "segment_weight": pl.Float64 
            }
        )

        # Validation: Handle empty or 5-column "no result" files
        if raw_df.width < 12 or raw_df.height == 0:
            print(f"{sample_name}: No data or partial header. Adding dummy row.")
    
            # Create a dummy row
            dummy_df = nw.from_dict({
                col: [sample_name if col == SAMPLE_COL else None] 
                for col in EXPECTED_COLUMNS
            }, backend="polars")

            # Force types to match dfs generated from genemetrics files with actual data
            # This prevents data type compatibility errors during vertical concat
            return dummy_df.with_columns([
                nw.col(GENE_COL).cast(nw.String),
                nw.col(CHROMOSOME_COL).cast(nw.Int8),
                nw.col("start").cast(nw.Int64),
                nw.col("end").cast(nw.Int64),
                nw.col(LOG2_COL).cast(nw.Float64),
                nw.col("depth").cast(nw.Float64),
                nw.col("weight").cast(nw.Float64),
                nw.col("ci_hi").cast(nw.Float64),
                nw.col("ci_lo").cast(nw.Float64),
                nw.col("probes").cast(nw.Int64),
                nw.col("segment_weight").cast(nw.Float64),
                nw.col("segment_probes").cast(nw.Int64),
                nw.col(SAMPLE_COL).cast(nw.String)
            ])

        # Standard Processing for genemetrics files with actual data
        df = nw.from_native(raw_df)
        
        # Cast and Filter Chromosomes (XY removal)
        # Important: CHOMOSOME_COL as Int8 data type
        df = df.with_columns(nw.col(CHROMOSOME_COL).cast(nw.String))
        df = df.filter(~nw.col(CHROMOSOME_COL).is_in(["X", "Y"]))
        df = df.with_columns(nw.col(CHROMOSOME_COL).cast(nw.Int8))
        
        # Add Sample Column
        df = df.with_columns(nw.lit(sample_name).alias(SAMPLE_COL))
        
        # Ensure we return exactly the columns needed for vertical concat
        return df.select(EXPECTED_COLUMNS)

    except Exception as e:
        print(f"Error reading {sample_name}: {e}")
        # Return a dummy row so the entire report doesn't fail
        return nw.from_dict({
            col: [sample_name if col == SAMPLE_COL else None] 
            for col in EXPECTED_COLUMNS
        }, backend="polars")


def calculate_copy_number(original_df: nw.DataFrame, log2col=LOG2_COL) -> nw.DataFrame:
  """
  Take a DataFrame with log2 ratios and calculate the copy numbers.
  :return nw.DataFrame: return an uncollected DataFrame with a new column of copy_number.
  """
  original_df = original_df.with_columns((DIPLOID * (2.0 ** nw.col(log2col))).round(decimals=0).alias(COPY_NUMBER_COL))
  return original_df


def sum_cnv_numbers_by_gene(all_gene_cnvs: nw.DataFrame) -> nw.DataFrame:
  # Log 2 ratio means 2.0 to the power of the log2ratio
  all_gene_cnvs = calculate_copy_number(all_gene_cnvs, log2col=LOG2_COL)
  all_gene_summary = all_gene_cnvs.group_by(GENE_COL, SAMPLE_COL).agg(
      nw.when(nw.col(COPY_NUMBER_COL) > DIPLOID).then(1).otherwise(0).sum().alias(GAIN_COL),
      nw.when(nw.col(COPY_NUMBER_COL) < DIPLOID).then(1).otherwise(0).sum().alias(LOSS_COL)
  )
  is_gain_and_loss = (nw.col(GAIN_COL) == 1) & (nw.col(LOSS_COL) == 1)
  all_gene_summary = all_gene_summary.with_columns(
      nw.when(is_gain_and_loss).then(1).otherwise(0).alias(GAIN_AND_LOSS_COL),
      nw.when(is_gain_and_loss).then(0).otherwise(nw.col(GAIN_COL)).alias(GAIN_COL),
      nw.when(is_gain_and_loss).then(0).otherwise(nw.col(LOSS_COL)).alias(LOSS_COL),
  )
  all_gene_summary = all_gene_summary.group_by(GENE_COL).agg(
      nw.col(GAIN_COL).sum(),
      nw.col(LOSS_COL).sum(),
      nw.col(GAIN_AND_LOSS_COL).sum()
  )
  return all_gene_summary


async def summarize_gene_level(gene_files: list[str]) -> nw.DataFrame:
  """
  Summarizes all the calls (gains/losses) by gene.
  :param gene_files list[str]: A list of file paths ending in .genemetrics.tsv
  :return nw.DataFrame: A data frame with the columns: gene, No. of Samples with Gains, No. of Samples with Losses,
      and No. of Samples with Gains & Losses.
  """
  all_gene_cnvs = nw.from_native(pl.DataFrame({}))
  for one_gene_file in gene_files:
    one_sample_variants = get_one_sample_variants(one_gene_file)
    all_gene_cnvs = nw.concat([all_gene_cnvs, one_sample_variants], how="vertical")
  all_gene_summary = sum_cnv_numbers_by_gene(all_gene_cnvs)
  return all_gene_summary
