#!/usr/local/miniforge3/bin/python

"""
Creates the "SNV/Indel Variants" Sheet.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import polars as pl
import narwhals.stable.v2 as nw
from .vcf2df import convert_vcf_to_data_frame


def gather_variants(vcf_files: list[str]) -> nw.DataFrame:
  all_sample_variants = nw.from_native(pl.DataFrame({}))
  for one_file in vcf_files:
    one_sample_variants = convert_vcf_to_data_frame(one_file)
    all_sample_variants = nw.concat([all_sample_variants, one_sample_variants], how="diagonal")
  return all_sample_variants
