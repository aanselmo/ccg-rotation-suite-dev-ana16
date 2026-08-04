#!/usr/local/miniforge3/bin/python

"""
Download sarek, picard, and somalier data from google cloud and organize for reporting.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie"]


import re
import narwhals as nw
import polars as pl
from os.path import basename
from ..setup.google import *
from .create_html_report import parse_somalier_pairs_tsv
from .somalier_processing import rename_somalier_samples
from ..common import debug_df

RELATEDNESS_COL = "relatedness"


def get_google_picard(storage_client, picard_uri: list[str]):
  """
  Download picard hsmetrics data from google cloud and organize for reporting.
  """
  sample_names = [""] * len(picard_uri)
  picard_outputs = [""] * len(picard_uri)

  for this_picard_file in picard_uri:
    this_blob = get_blob_from_gs_uri(storage_client, this_picard_file)
    this_sample_name = re.sub(f".*/picard/(.*)/hsmetrics.*", r"\1", this_picard_file)

    row_index = picard_uri.index(this_picard_file)
    sample_names[row_index] = this_sample_name
    this_picard = pl.read_csv(this_blob.download_as_string(), separator="\t", skip_lines=6, n_rows=1).select(
        "TOTAL_READS",
        "PF_READS",
        "PCT_PF_READS",
        "PCT_SELECTED_BASES",
        pl.col("MEAN_TARGET_COVERAGE"),
        "FOLD_ENRICHMENT",
        pl.col("PCT_EXC_DUPE").alias("PERCENT_DUPLICATION"),
        pl.col("HS_LIBRARY_SIZE").alias("ESTIMATED_LIBRARY_SIZE"),
        pl.col("PCT_TARGET_BASES_20X").alias("Pct_full_targets_Bases_20x"),
        pl.col("PCT_TARGET_BASES_30X").alias("Pct_full_targets_Bases_30x"),
        pl.col("PCT_TARGET_BASES_50X").alias("Pct_full_targets_Bases_50x"),
        pl.col("PCT_TARGET_BASES_100X").alias("Pct_full_targets_Bases_100x"),
        pl.col("PCT_PF_UQ_READS_ALIGNED"),
        pl.col("ZERO_CVG_TARGETS_PCT").alias("Pct_full_targets_Not_Covered"),
        pl.col("PCT_USABLE_BASES_ON_BAIT"))
    picard_outputs[row_index] = this_picard

  picard_outputs_table = pl.concat(picard_outputs, how='vertical')

  picard_outputs_table_names = picard_outputs_table.with_columns(pl.Series("sample_collaborator_name", sample_names))
  cols_ordered = ["sample_collaborator_name"] + picard_outputs_table.columns
  picard_outputs_table_sorted = picard_outputs_table_names.select(cols_ordered)
  picard_metrics = nw.from_native(picard_outputs_table_sorted)
  debug_df("picard_metrics", picard_metrics.to_native(), n=10)
  return (picard_metrics)


def get_google_somalier_pairs(storage_client, somalier_pairs_uri, data_lims):
  """
  Download somalier pairs data from google cloud and organize for reporting
  """

  somalier_blob = get_blob_from_gs_uri(storage_client, somalier_pairs_uri)
  somalier_pairs_raw = pl.read_csv(
      somalier_blob.download_as_string(), separator="\t", schema_overrides={RELATEDNESS_COL: pl.Float64})
  somalier_data = parse_somalier_pairs_tsv(somalier_pairs_raw)

  pl_somalier_pairs = (
      nw.to_native(somalier_data)
      .select('#sample_a', 'sample_b', 'Match Status', 'relatedness', 'ibs0', 'ibs2', 'expected_relatedness')
      .rename({"#sample_a": "Sample1", "sample_b": "Sample2", "Match Status": "Type",
               "relatedness": "Relatedness", "expected_relatedness": "Expected Relatedness"}))

  if data_lims:
    pl_somalier_pairs = rename_somalier_samples(nw.from_native(pl_somalier_pairs), data_lims).to_polars()

  somalier_pairs = nw.from_native(pl_somalier_pairs)
  debug_df("somalier_pairs", somalier_pairs.to_native(), n=10)
  pl_somalier_mismatch = pl_somalier_pairs.filter(pl.col("Type") != "OK")
  somalier_mismatch = nw.from_native(pl_somalier_mismatch)
  debug_df("somalier_mismatch", somalier_mismatch.to_native(), n=10)
  return (somalier_pairs, somalier_mismatch)


def get_google_somalier_samples(storage_client, somalier_samples_uri, data_lims):
  """
  Download somalier samples data from google cloud and organize for reporting
  """
  
  somalier_samples_blob = get_blob_from_gs_uri(storage_client, somalier_samples_uri)
  somalier_samples_raw = pl.read_csv(somalier_samples_blob.download_as_string(), separator="\t")

  pl_somalier_samples = (
    somalier_samples_raw
      .select('sample_id', 'gt_depth_mean', 'n_hom_ref', 'n_het', 'n_hom_alt', 'n_unknown')
      .with_columns(n_sites = pl.col('n_hom_ref') + pl.col('n_het') + pl.col('n_hom_alt'))
      .with_columns(n_possible_sites = pl.col('n_sites') + pl.col('n_unknown'))
  )
  if data_lims:
    pl_somalier_samples = rename_somalier_samples(nw.from_native(pl_somalier_samples), data_lims).to_polars()
  return(nw.from_native(pl_somalier_samples))


def get_google_duplicates(storage_client, duplicate_marking_uri: list[str]) -> nw.DataFrame:
  """
  Download markduplicates md.cram.metrics data from google cloud and organize for reporting
  """

  all_pct_duplication = nw.from_native(pl.DataFrame({}))

  for file_path in duplicate_marking_uri:
    sample_name = basename(file_path).replace(".md.cram.metrics", "")
    pct_dup = 0.0
    duplicate_blob = get_blob_from_gs_uri(storage_client, file_path)
    one_file_metrics = pl.read_csv(
        duplicate_blob.download_as_string(), separator="\t", skip_lines=6, n_rows=1)
    pct_dup = one_file_metrics.get_column("PERCENT_DUPLICATION")[0]
    one_row = {
        "Sample ID": sample_name,
        "% Duplication": pct_dup * 100.0,
    }
    all_pct_duplication = nw.concat([all_pct_duplication, nw.from_dict(one_row, backend="polars")], how="vertical")
  debug_df("all_pct_duplication", all_pct_duplication)
  return all_pct_duplication


def get_google_lanes(storage_client, samplesheet_uri):
  """
  Download sample sheet and count number of lanes per sample
  """
  samplesheet_blob = get_blob_from_gs_uri(storage_client, samplesheet_uri)
  pl_samplesheet = pl.read_csv(samplesheet_blob.download_as_string())
  pl_lanes_per_sample = pl_samplesheet.group_by("sample").len().rename({"len": "aggregate_count"})
  debug_df("lanes_per_sample", pl_lanes_per_sample)
  lanes_per_sample = nw.from_native(pl_lanes_per_sample)
  return(lanes_per_sample)
