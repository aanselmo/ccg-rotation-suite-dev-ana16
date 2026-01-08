#!/usr/local/miniforge3/bin/python

"""
Creates the 'SNV/Indel Variants' Sheet.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os.path import basename
import polars as pl
import narwhals.stable.v2 as nw
from .vcf2df import convert_vcf_to_data_frame


VEP_FILE_ENDING = ".mutect2.filtered.bcftools_filtered_VEP.ann.vcf.gz"
SAMPLE_COL = "tumor_sample_name"
REF_COL = "REF"
ALT_COL = "ALT"
RECOMMENDED_FILTER_COL = "recommended_filter"
CONSEQUENCE_COL = "Consequence"
AMINO_ACIDS_COL = "Amino_acids"
EXISTING_VARIATION_COL = "Existing_variation"
VARIANT_CLASS_COL = "VARIANT_CLASS"
MAX_AF_POPS_COL = "MAX_AF_POPS"
RECOMMENDED_CONSEQUENCE_FILTER = [
    "missense_variant", "frameshift_variant", "inframe_deletion", "inframe_insertion",
    "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant",
    "stop_lost", "start_lost", "stop_gained", "splice_region_variant",
    "transcript_ablation", "transcript_amplification"
]


def add_recommended_consequence_filter(all_sample_variants: nw.DataFrame) -> nw.DataFrame:
  all_sample_variants = all_sample_variants.with_columns(
      nw.when(nw.col(CONSEQUENCE_COL).is_in(RECOMMENDED_CONSEQUENCE_FILTER))
      .then(True)
      .otherwise(False)
      .alias(RECOMMENDED_FILTER_COL)
  )
  return all_sample_variants


def get_variants(vcf_files: list[str]) -> nw.DataFrame:
  all_sample_variants = nw.from_native(pl.DataFrame({}))
  for one_file in vcf_files:
    one_sample_variants = convert_vcf_to_data_frame(one_file)
    sample_name = basename(one_file).replace(VEP_FILE_ENDING, "")
    one_sample_variants = one_sample_variants.with_columns(nw.lit(sample_name).alias(SAMPLE_COL))
    all_sample_variants = nw.concat([all_sample_variants, one_sample_variants], how="diagonal")
  all_sample_variants = add_recommended_consequence_filter(all_sample_variants)
  # Select only these columns for now. There are a lot more columns. I made an arbitrary call to get it as similar to
  # the riker report as possible. Also, there a lot of gnomAD columns so I tried to summarize them.
  all_sample_variants = all_sample_variants.select([
      SAMPLE_COL, "#CHROM", "POS", REF_COL, ALT_COL, RECOMMENDED_FILTER_COL, CONSEQUENCE_COL, "CANONICAL",
      VARIANT_CLASS_COL, "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "Gene",
      "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "cDNA_position", "CDS_position",
      "Protein_position", AMINO_ACIDS_COL, EXISTING_VARIATION_COL, "AF", "DP",
      "gnomAD_AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF", MAX_AF_POPS_COL,
  ])
  return all_sample_variants
