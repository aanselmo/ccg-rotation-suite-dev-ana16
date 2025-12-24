
#!/usr/local/miniforge3/bin/python

"""
Gathers the data for the Summary Sheet by Chromosomal Arm in the CNV report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os.path import dirname, basename
import asyncio
import polars as pl
import narwhals.stable.v2 as nw
from qc_sheet import BASENAME_COL
from variant_sheet import CHROMOSOME_COL, START_COL, END_COL, COPY_NUMBER_COL
from gene_sheet import GAIN_AND_LOSS_COL, GAIN_COL, LOSS_COL, DIPLOID

ARM_COL = "arm"


async def summarize_chromosome_arm_level(chromosomal_arms: nw.LazyFrame, cnv_variants: nw.LazyFrame) -> nw.LazyFrame:
  """
  Summarizes all the calls (gains/losses) by chromosomal arm.
  :param cnv_variants nw.LazyFrame: All the variants in all samples excluding those on the X & Y chromosome.
  """
  summary_chromosomal_arm = chromosomal_arms.with_columns(
      nw.col(0).alias(GAIN_COL),
      nw.col(0).alias(LOSS_COL),
      nw.col(0).alias(GAIN_AND_LOSS_COL),
  )
  all_samples = cnv_variants.select(nw.col(BASENAME_COL)).group_by(BASENAME_COL)
  summary_chromosomal_arm = summary_chromosomal_arm.join(other=all_samples, how="cross")
  for variant in cnv_variants.iter_rows():  # Or similar function
    (chromosome, start, end, gene, log2, basename, gene, ci_hi, ci_lo,
     copy_number, depth, probes, weight) = variant
    column_to_add_to = ""
    if copy_number > DIPLOID:
      column_to_add_to = GAIN_COL
    if copy_number < DIPLOID:
      column_to_add_to = LOSS_COL
    # Because of half-open coordinates:
    # https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
    # +1 if cnv_variant.start is in between (>= or <) start & end of p arm
    # +1 if cnv_variant.end is in between (>= or <) start & end of q arm
    arm_to_add_to = ""
    chromosomal_arms.filter((nw.col(CHROMOSOME_COL) == chromosome) & (nw.col(CHROMOSOME_COL) == chromosome))
    summary_chromosomal_arm.with_columns(
        pl.when(
            (pl.col(BASENAME_COL) == basename)
            & (pl.col(CHROMOSOME_COL) == chromosome)
            & (pl.col(ARM_COL) == )
        ).then(pl.col(column_to_add_to) + 1)
        .otherwise(pl.col(column_to_add_to))
        .alias(column_to_add_to)
    )
