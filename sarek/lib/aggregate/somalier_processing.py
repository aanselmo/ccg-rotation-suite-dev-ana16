#!/usr/local/miniforge3/bin/python

"""
Combine and format some somalier outputs for reporting.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie"]

import polars as pl
import narwhals as nw

def append_somalier_sample_info_to_mismatch(somalier_mismatch, somalier_samples):
  somalier_samples_select = somalier_samples.to_polars().select('sample_id', pl.col('gt_depth_mean').alias("mean_depth"), 'n_sites')
  somalier_mismatch_annotated = (
    somalier_mismatch.to_polars()
      .join(somalier_samples_select.rename(lambda col: f"Sample1_{col}"), 
            left_on = "Sample1", right_on = "Sample1_sample_id")
      .join(somalier_samples_select.rename(lambda col: f"Sample2_{col}"), 
            left_on = "Sample2", right_on = "Sample2_sample_id")
  )
  return(nw.from_native(somalier_mismatch_annotated))


def rename_somalier_samples(somalier_data, data_lims):
  """
  Somalier sample names are a concatenation of LIMS sample name and collaborator sample name. For similarity to Riker,
  the order of these names is being flipped.  
  """
  data_lims = data_lims.with_columns(
    somalier_paste = nw.col('patient') + '_' + nw.col('sample'),
    somalier_paste_invert = nw.col('sample') + '_' + nw.col('patient_original'),
  )
  pl_somalier_data = somalier_data.to_polars()
  pl_data_lims = data_lims.to_polars()

  if 'Sample1' in pl_somalier_data.columns:
      pl_somalier_data = (
            pl_somalier_data
            .join(
            pl_data_lims.select('somalier_paste', pl.col('somalier_paste_invert').alias('Sample1_rename')),
            left_on='Sample1', right_on='somalier_paste', how='left'
            )
            .join(
            pl_data_lims.select('somalier_paste', pl.col('somalier_paste_invert').alias('Sample2_rename')),
            left_on='Sample2', right_on='somalier_paste', how='left'
            )
            .drop('Sample1', 'Sample2')
            .rename({"Sample1_rename": "Sample1", "Sample2_rename": "Sample2"})
            .select(
                  "Sample1",
                  "Sample2",
                  *[c for c in pl_somalier_data.columns if c not in ("Sample1", "Sample2")]
            )
      )

  if 'sample_id' in pl_somalier_data.columns:
      pl_somalier_data = (
            pl_somalier_data
            .join(
            pl_data_lims.select('somalier_paste', pl.col('somalier_paste_invert').alias('sample_id_rename')),
            left_on='sample_id', right_on='somalier_paste', how='left'
            )
            .drop('sample_id')
            .rename({"sample_id_rename": "sample_id"})
            .select(
                  "sample_id",
                  *[c for c in pl_somalier_data.columns if c not in ("sample_id")]
            )
      )
  return nw.from_native(pl_somalier_data)

