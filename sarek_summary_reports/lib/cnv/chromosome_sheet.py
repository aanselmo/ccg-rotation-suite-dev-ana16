
#!/usr/local/miniforge3/bin/python

"""
Gathers the data for the Summary Sheet by Chromosomal Arm in the CNV report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import gzip
import polars as pl
import narwhals.stable.v2 as nw
from .qc_sheet import BASENAME_COL
from .variant_sheet import CHROMOSOME_COL, START_COL, END_COL, COPY_NUMBER_COL
from .gene_sheet import GAIN_AND_LOSS_COL, GAIN_COL, LOSS_COL, DIPLOID

ARM_COL = "arm"

def derive_chromosomal_arms(genome_version: str) -> nw.DataFrame:
  """
  Derives the chromosomal arms from USCS cytoband file. Downloaded from:
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
  UCSC uses a half open system for coordinates for 0 based counting:
  https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
  :return nw.DataFrame: A dataframe with the columns "chromosome", "arm", "start", "end"
  """
  l_genome = genome_version.lower()
  cytoband_data = nw.from_native(pl.DataFrame({}))
  with gzip.open(f'{l_genome}CytoBand.txt.gz', 'rb') as cytoband_file:
    cytoband_data = pl.scan_csv(
        cytoband_file,
        separator="\t",
        has_header=False,
        new_columns=[CHROMOSOME_COL, START_COL, END_COL, "band", "stain"]
    )
    # Aggregate by arm substring
    cytoband_data = nw.from_native(cytoband_data)
    cytoband_data = cytoband_data.with_columns(arm=nw.col("band").str.slice(offset=0, length=1))  # first letter: p/q
    cytoband_data = cytoband_data.drop("band", "stain")
    cytoband_data = cytoband_data.group_by(CHROMOSOME_COL, ARM_COL).agg(nw.col(START_COL).min(), nw.col(END_COL).max())
  return cytoband_data


def reduce_data_memory(chromosomal_arms: nw.DataFrame, cnv_variants: nw.DataFrame) -> tuple[nw.DataFrame, nw.DataFrame]:
  """
  Drop some columns and rows to save memory in the next step (a matrix cross product)
  """
  chromosomal_arms = chromosomal_arms.filter((nw.col(CHROMOSOME_COL) != "chrX") & (nw.col(CHROMOSOME_COL) != "chrY"))
  chromosomal_arms = chromosomal_arms.with_columns(
      nw.col(CHROMOSOME_COL).str.replace("chr", "").cast(nw.Int8).alias(CHROMOSOME_COL)
  )
  cnv_variants = cnv_variants.select([
      nw.col(CHROMOSOME_COL),
      nw.col(START_COL),
      nw.col(END_COL),
      nw.col(COPY_NUMBER_COL),
      nw.col(BASENAME_COL)
  ])
  return (chromosomal_arms, cnv_variants)


def assign_variants_to_arms(chromosomal_arms: nw.DataFrame, cnv_variants: nw.DataFrame) -> nw.DataFrame:
  """
  Takes each variant and assigns it to chromosomal arm(s).
  """
  # The cross product is memory intensive, but computationally not. Easy for polars to do matrix math.
  # I think it'll take less than 5 MB of memory
  variants_per_arm = (cnv_variants.join(chromosomal_arms, how="cross"))
  # Categorize gain or loss to sum later.
  # Because of half-open coordinates:
  # https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
  # +1 if cnv_variant.start is in between (>= or <) start & end of p arm
  # +1 if cnv_variant.end is in between (>= or <) start & end of q arm
  variants_per_arm = variants_per_arm.with_columns(
      nw.when(
          (nw.col(COPY_NUMBER_COL) > DIPLOID)
          & (nw.col(CHROMOSOME_COL) == nw.col("chromosome_right"))
          & (nw.col(START_COL) >= nw.col("start_right"))
          & (nw.col(START_COL) < nw.col("end_right"))
      ).then(1).otherwise(0).alias(GAIN_COL),
      nw.when(
          (nw.col(COPY_NUMBER_COL) < DIPLOID)
          & (nw.col(CHROMOSOME_COL) == nw.col("chromosome_right"))
          & (nw.col(END_COL) >= nw.col("start_right"))
          & (nw.col(END_COL) < nw.col("end_right"))
      ).then(1).otherwise(0).alias(LOSS_COL)
  )
  return variants_per_arm


def sum_variants_per_arm(variants_per_arm: nw.DataFrame) -> nw.DataFrame:
  """
  Sum to create summary table.
  """
  variants_per_arm = variants_per_arm.group_by(
      BASENAME_COL, "chromosome_right", ARM_COL, "start_right", "end_right"
  ).agg(
      nw.col(GAIN_COL).sum(),
      nw.col(LOSS_COL).sum()
  )
  variants_per_arm = variants_per_arm.with_columns(
      nw.when(nw.col(GAIN_COL) > 0).then(1).otherwise(0).alias(GAIN_COL),
      nw.when(nw.col(LOSS_COL) > 0).then(1).otherwise(0).alias(LOSS_COL),
  )
  is_gain_and_loss = (nw.col(GAIN_COL) == 1) & (nw.col(LOSS_COL) == 1)
  variants_per_arm = variants_per_arm.with_columns(
      nw.when(is_gain_and_loss).then(1).otherwise(0).alias(GAIN_AND_LOSS_COL),
      nw.when(is_gain_and_loss).then(0).otherwise(nw.col(GAIN_COL)).alias(GAIN_COL),
      nw.when(is_gain_and_loss).then(0).otherwise(nw.col(LOSS_COL)).alias(LOSS_COL),
  )
  variants_per_arm = variants_per_arm.group_by("chromosome_right", ARM_COL, "start_right", "end_right").agg(
      nw.col(GAIN_COL).sum(),
      nw.col(LOSS_COL).sum(),
      nw.col(GAIN_AND_LOSS_COL).sum(),
  )
  variants_per_arm = variants_per_arm.rename({
      "chromosome_right": CHROMOSOME_COL,
      "start_right": START_COL,
      "end_right": END_COL
  })
  variants_per_arm = variants_per_arm.sort(by=[CHROMOSOME_COL, ARM_COL])
  return variants_per_arm


async def summarize_chromosome_arm_level(genome_version: str, cnv_variants: nw.DataFrame) -> nw.DataFrame:
  """
  Summarizes all the calls (gains/losses) by chromosomal arm.
  :param cnv_variants nw.DataFrame: All the variants in all samples excluding those on the X & Y chromosome.
  """
  chromosomal_arms = derive_chromosomal_arms(genome_version)
  (chromosomal_arms, cnv_variants) = reduce_data_memory(chromosomal_arms.lazy(), cnv_variants)
  variants_per_arm = assign_variants_to_arms(chromosomal_arms, cnv_variants)
  variants_per_arm = sum_variants_per_arm(variants_per_arm)
  return variants_per_arm
