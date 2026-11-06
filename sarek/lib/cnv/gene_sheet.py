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
# log2 cutoffs of cnvkit's `call --method threshold` (defaults -1.1,-0.25,0.2,0.7): log2 <= -0.25 is a loss
# (cn <= 1) and log2 > 0.2 is a gain (cn >= 3). These are the cutoffs behind the cn column of the call.cns files
# used for the CNV tab (and the same cutoffs used for the arm-level calls and the CNV plots), so using them here keeps
# the Gene Level Summary consistent with the CNV tab.
# Do NOT derive copy number as round(2 * 2**log2): that rule needs log2 >= 0.32 for a gain and <= -0.42 for a loss,
# so modest gains/losses that cnvkit calls (e.g. log2 0.20-0.32) were counted as neutral.
GAIN_LOG2_CUTOFF = 0.2
LOSS_LOG2_CUTOFF = -0.25



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
  # Classify each gene row with the same log2 cutoffs cnvkit's threshold calling uses (see GAIN_LOG2_CUTOFF).
  is_gain = nw.col(LOG2_COL) > GAIN_LOG2_CUTOFF
  is_loss = nw.col(LOG2_COL) <= LOSS_LOG2_CUTOFF
  # A gene can be listed more than once for one sample when a breakpoint falls inside it (one row per segment it
  # overlaps). Count each sample at most once per gene: use max (i.e. "any row is a gain/loss"), not sum, otherwise
  # such a sample is counted twice and the "Gains & Losses" test (== 1) misses it.
  all_gene_summary = all_gene_cnvs.group_by(GENE_COL, SAMPLE_COL).agg(
      nw.when(is_gain).then(1).otherwise(0).max().alias(GAIN_COL),
      nw.when(is_loss).then(1).otherwise(0).max().alias(LOSS_COL)
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


SAMPLES_WITH_GAINS_COL = "Samples with Gains"
SAMPLES_WITH_LOSSES_COL = "Samples with Losses"
SAMPLE_LIST_SEPARATOR = ", "
# cnvkit's placeholder for off-target (antitarget) regions in the call.cns gene column; not a gene.
NOT_A_GENE = "-"


def add_sample_lists(gene_summary: nw.DataFrame, cnv_variants: nw.DataFrame) -> nw.DataFrame:
  """
  Adds two columns to the gene summary listing, for each gene, the samples that have it in a gain row / a loss row of
  the CNV tab (the comma-separated gene column of the call.cns variants). A sample whose CNV-tab rows put the gene in
  both a gain and a loss (gene split across two called segments) is listed in both columns; such samples are the ones
  counted under "No. of Samples with Gains & Losses".
  :param gene_summary nw.DataFrame: output of summarize_gene_level().
  :param cnv_variants nw.DataFrame: output of get_cnv_variants() (the CNV tab).
  :return nw.DataFrame: gene_summary with SAMPLES_WITH_GAINS_COL and SAMPLES_WITH_LOSSES_COL appended
      (null where no sample qualifies).
  """
  pl_summary = gene_summary.to_native()
  if not pl_summary.width:
    return gene_summary
  pl_variants = cnv_variants.to_native()
  needed = {GENE_COL, SAMPLE_COL, COPY_NUMBER_COL}
  if not pl_variants.width or not needed.issubset(pl_variants.columns):
    return nw.from_native(pl_summary.with_columns(
        pl.lit(None, dtype=pl.String).alias(SAMPLES_WITH_GAINS_COL),
        pl.lit(None, dtype=pl.String).alias(SAMPLES_WITH_LOSSES_COL),
    ))

  genes_by_sample = (
      pl_variants
      .select(
          pl.col(GENE_COL).cast(pl.String).str.split(",").alias(GENE_COL),
          pl.col(SAMPLE_COL).cast(pl.String),
          pl.col(COPY_NUMBER_COL).cast(pl.Float64),
      )
      .explode(GENE_COL)
      .with_columns(pl.col(GENE_COL).str.strip_chars())
      .filter(pl.col(GENE_COL).is_not_null() & (pl.col(GENE_COL) != "") & (pl.col(GENE_COL) != NOT_A_GENE))
  )

  def samples_where(condition: pl.Expr, out_col: str) -> pl.DataFrame:
    return (
        genes_by_sample.filter(condition)
        .select(GENE_COL, SAMPLE_COL).unique()
        .group_by(GENE_COL)
        .agg(pl.col(SAMPLE_COL).sort().str.join(SAMPLE_LIST_SEPARATOR).alias(out_col))
    )

  gains = samples_where(pl.col(COPY_NUMBER_COL) > DIPLOID, SAMPLES_WITH_GAINS_COL)
  losses = samples_where(pl.col(COPY_NUMBER_COL) < DIPLOID, SAMPLES_WITH_LOSSES_COL)
  pl_summary = (
      pl_summary.with_columns(pl.col(GENE_COL).cast(pl.String))
      .join(gains, on=GENE_COL, how="left")
      .join(losses, on=GENE_COL, how="left")
  )
  return nw.from_native(pl_summary)
