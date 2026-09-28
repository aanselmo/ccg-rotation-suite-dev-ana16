#!/usr/local/miniforge3/bin/python

"""
Creates additional tables for PDX samples that have:
1. metrics on mouse contamination.
2. blacklisted Human Aligned Mouse Alleles (HAMAs) that should be filtered out.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import sys
from os.path import basename, dirname
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import polars as pl
import narwhals.stable.v2 as nw
script_dir = Path(__file__).resolve().parent
sarek_script_top_dir = dirname(dirname(script_dir))
sys.path.append(sarek_script_top_dir)
from lib.snv_indel.snv_indel_sheet import SAMPLE_COL, REF_COL, ALT_COL
from lib.snv_indel.vcf2df import convert_vcf_to_data_frame, CHROM_COL, POSITION_COL
from lib.snv_indel.lego_plot_sheet import UNFILTERED_VCF_FILE_ENDING


PRIMARY_READS_COL = "primary reads"
NCBIM37_READS_COL = "NCBIM37 reads"
TOTAL_READS_COL = "total reads"
CONTAMINATION_PCT_COL = "% mouse reads"
MOUSE_POSITIVE_CONTROL_THRESHOLD = .995

SEEN_IN_MOUSE_COL = "Seen_In_Mouse"


def get_contamination_metrics(mouse_bbsplit_file: str) -> nw.DataFrame:
  """
  Get the amount of mouse contamination in each sample.
  """
  # Number of chunks are define by https://nf-co.re/sarek/3.7.1/docs/usage/#split-fastq-files
  # Typically, there are 12 chunks that a fastq is split into.
  CHUNK_NUM_COL = "chunk_num"
  LANE_COL = "lane"
  contamination = pl.read_csv(mouse_bbsplit_file, separator="\t")
  contamination = contamination.rename({
      "Sample": SAMPLE_COL,
      "primary": PRIMARY_READS_COL,
      "NCBIM37": NCBIM37_READS_COL
  })
  contamination = contamination.cast({
    PRIMARY_READS_COL: pl.Int64,
    NCBIM37_READS_COL: pl.Int64
  })
  # Example original "Sample" col = "0001.CHGFR0004-PDX-B1-L001"
  # 0001 is the chunk and L001 is the lane entered in the sarek samplesheet.
  contamination = contamination.with_columns(
      pl.col(SAMPLE_COL).str.split_exact(".", 1)
      .struct.rename_fields([CHUNK_NUM_COL, SAMPLE_COL])
      .struct.unnest()
  )
  contamination = contamination.with_columns(
      pl.col(SAMPLE_COL).str.split("-")
  )
  contamination = contamination.with_columns(
      # only the last item (-1)
      pl.col(SAMPLE_COL).list.get(-1).alias(LANE_COL),
      # everything but the last item (-1)
      pl.col(SAMPLE_COL).list.head(pl.col(SAMPLE_COL).list.len()-1).list.join("-").alias(SAMPLE_COL)
  )
  contamination = contamination.group_by(SAMPLE_COL).agg(
      pl.col(LANE_COL),
      pl.col(CHUNK_NUM_COL),
      pl.col(PRIMARY_READS_COL).sum(),
      pl.col(NCBIM37_READS_COL).sum(),
  )
  contamination = contamination.with_columns(
      pl.col(PRIMARY_READS_COL).add(pl.col(NCBIM37_READS_COL)).alias(TOTAL_READS_COL)
  )
  contamination = contamination.with_columns(
      pl.col(NCBIM37_READS_COL).truediv(pl.col(TOTAL_READS_COL))
      .cast(pl.Float64).alias(CONTAMINATION_PCT_COL)
  )
  contamination = contamination.drop([CHUNK_NUM_COL, LANE_COL])
  return nw.from_native(contamination)


def get_mouse_pos_control_vcfs(
    mouse_positive_vcfs: list[str], unfiltered_vcf_files: list[str], contamination_table: nw.DataFrame
) -> tuple[list[str], list[str]]:
  """
  If mouse_positive_vcfs is empty list, dynamically determine which samples should be used as mouse positive controls
  based on contamination levels.
  """
  mouse_pos_controls = []
  if not mouse_positive_vcfs:
    mouse_pos_controls = contamination_table.filter(
        nw.col(CONTAMINATION_PCT_COL) > MOUSE_POSITIVE_CONTROL_THRESHOLD
    ).get_column(SAMPLE_COL).to_list()
    for vcf_path in unfiltered_vcf_files:
      for sample_name in mouse_pos_controls:
        if sample_name in vcf_path:
          mouse_positive_vcfs.append(vcf_path)
    print("Warning! No Mouse Positive controls supplied with --mouse-positive-vcfs flag.")
    print((f"Using the following samples (with >{MOUSE_POSITIVE_CONTROL_THRESHOLD} contamination)"
          + f" as mouse positive controls: {mouse_pos_controls}."))
  else:
    for one_vcf_file_path in unfiltered_vcf_files:
      sample_name = basename(one_vcf_file_path).replace(UNFILTERED_VCF_FILE_ENDING, "")
      mouse_pos_controls.append(sample_name)
  return (mouse_pos_controls, mouse_positive_vcfs)


def process_one_mouse_file(one_vcf_file_path: str) -> nw.DataFrame:
  """
  Run vcf2df for one Human Aligned Mouse Alleles (HAMAs) unfiltered vcf file.
  """
  one_sample_variants = convert_vcf_to_data_frame(one_vcf_file_path, should_filter_mutect2_failed=False)
  if not one_sample_variants.is_empty():
    sample_name = basename(one_vcf_file_path).replace(UNFILTERED_VCF_FILE_ENDING, "")
    one_sample_variants = one_sample_variants.with_columns(nw.lit(sample_name).alias(SAMPLE_COL))
  return one_sample_variants


def get_human_aligned_mouse_alleles(
    mouse_positive_vcfs: list[str], unfiltered_vcf_files: list[str], contamination_table: nw.DataFrame,
    number_processes=1
) -> tuple[nw.DataFrame, list[str]]:
  """
  Get a table of Human Aligned Mouse Alleles (HAMAs) that make it past bbsplit to use as a list of blacklisted variants.
  """
  (mouse_pos_controls, mouse_positive_vcfs) = get_mouse_pos_control_vcfs(
      mouse_positive_vcfs, unfiltered_vcf_files, contamination_table
  )
  with ProcessPoolExecutor(max_workers=number_processes, mp_context=get_context("spawn")) as executor:
    list_of_variant_dfs = list(executor.map(process_one_mouse_file, mouse_positive_vcfs))
  all_hamas = nw.concat(list_of_variant_dfs, how="diagonal")
  all_hamas.sort(by=[SAMPLE_COL, CHROM_COL, POSITION_COL])
  return (all_hamas, mouse_pos_controls)


def remove_human_aligned_mouse_alleles(
    all_variants: nw.DataFrame, all_hamas: nw.DataFrame, mouse_pos_controls: list[str]
) -> nw.DataFrame:
  """
  Remove blacklisted Human Aligned Mouse Alleles (HAMAs) from all_variants and update the summary table.
  """
  VARIANT_ID_COLUMNS = [CHROM_COL, POSITION_COL, REF_COL, ALT_COL]
  all_hamas = all_hamas.select(VARIANT_ID_COLUMNS)
  all_hamas = all_hamas.with_columns(nw.lit(True).alias(SEEN_IN_MOUSE_COL))
  all_variants = all_variants.filter(~nw.col(SAMPLE_COL).is_in(mouse_pos_controls))
  all_variants = all_variants.join(all_hamas, on=VARIANT_ID_COLUMNS, how="left")
  all_variants = all_variants.with_columns(
      nw.when(nw.col(SEEN_IN_MOUSE_COL).is_null())
      .then(False)
      .otherwise(nw.col(SEEN_IN_MOUSE_COL))
      .alias(SEEN_IN_MOUSE_COL)
  )
  return all_variants
