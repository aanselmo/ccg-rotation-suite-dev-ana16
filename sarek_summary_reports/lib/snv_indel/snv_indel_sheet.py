#!/usr/local/miniforge3/bin/python

"""
Creates the "SNV/Indel Variants" Sheet.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os.path import basename
import polars as pl
import narwhals.stable.v2 as nw
from .vcf2df import convert_vcf_to_data_frame


VEP_FILE_ENDING = ".mutect2.filtered.bcftools_filtered_VEP.ann.vcf"
CONSEQUENCE_COL = "Consequence"
RECOMMENDED_FILTER_COL = "recommended_filter"
EXISTING_VARIATION_COL = "Existing_variation"
VARIANT_CLASS_COL = "VARIANT_CLASS"
MAX_GNOMAD_EXOME_AF_COL = "MAX_gnomADe_AF"
MAX_GNOMAD_EXOME_POP_COL = "MAX_gnomADe_AF_POPS"
MAX_GNOMAD_GENOME_AF_COL = "MAX_gnomADg_AF"
MAX_GNOMAD_GENOME_POP_COL = "MAX_gnomADg_AF_POPS"
SAMPLE_COL = "tumor_sample_name"
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


def find_gnomad_max_pop(all_sample_variants: nw.DataFrame) -> nw.DataFrame:
  """
  Find the gnomAD max allele fraction populations for the exome and genome cohorts. For some reason, VEP does this
  for the 1000 genomes cohort, but not the other 2.
  
  :param all_sample_variants nw.DataFrame: the variants from the VEP annotated vcf files in a DataFrame.
  :return nw.DataFrame: The same all_sample_variants, but with 4 columns added on:
      MAX_gnomADe_AF -> The max allele fraction in the exome cohort.
      MAX_gnomADe_AF_POPS -> What population the above allele fraction belongs to.
      MAX_gnomADg_AF -> The max allele fraction in the genome cohort.
      MAX_gnomADg_AF_POPS -> What population the above allele fraction belongs to.
  """
  GNOMAD_EXOME_COLS = [
      "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF",
      "gnomADe_MID_AF", "gnomADe_NFE_AF", "gnomADe_REMAINING_AF", "gnomADe_SAS_AF"
  ]
  GNOMAD_GENOME_COLS = [
      "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF", "gnomADg_EAS_AF",
      "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_REMAINING_AF", "gnomADg_SAS_AF"
  ]
  all_sample_variants = all_sample_variants.with_columns(
      nw.max_horizontal(GNOMAD_EXOME_COLS).alias(MAX_GNOMAD_EXOME_AF_COL)
  )
  all_sample_variants = all_sample_variants.with_columns(
    nw.struct(GNOMAD_EXOME_COLS)
      .map_elements(
          lambda s: max(s, key=s.get),
          return_dtype=nw.Utf8
      )
      .alias(MAX_GNOMAD_EXOME_POP_COL)
  )
  all_sample_variants = all_sample_variants.with_columns(
      nw.max_horizontal(GNOMAD_GENOME_COLS).alias(MAX_GNOMAD_GENOME_AF_COL)
  )
  all_sample_variants = all_sample_variants.with_columns(
    nw.struct(GNOMAD_GENOME_COLS)
      .map_elements(
          lambda s: max(s, key=s.get),
          return_dtype=nw.Utf8
      )
      .alias(MAX_GNOMAD_GENOME_POP_COL)
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
  all_sample_variants = find_gnomad_max_pop(all_sample_variants)
  # Select only these columns for now. There are a lot more columns. I made an arbitrary call to get it as similar to
  # the riker report as possible. Also, there a lot of gnomAD columns so I tried to summarize them.
  all_sample_variants = all_sample_variants.select([
      "#CHROM", "POS", "REF", "ALT", RECOMMENDED_FILTER_COL, CONSEQUENCE_COL, VARIANT_CLASS_COL, "SYMBOL",
      "SYMBOL_SOURCE", "HGNC_ID", "Gene",
      "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
      "Protein_position", "Amino_acids", "Codons", EXISTING_VARIATION_COL, "AF", "DP",
      "gnomAD_AF", "MAX_AF", "MAX_AF_POPS", "gnomADe_AF", "gnomADg_AF",
      MAX_GNOMAD_EXOME_AF_COL, MAX_GNOMAD_EXOME_POP_COL, MAX_GNOMAD_GENOME_AF_COL, MAX_GNOMAD_GENOME_POP_COL
  ])
  return all_sample_variants
