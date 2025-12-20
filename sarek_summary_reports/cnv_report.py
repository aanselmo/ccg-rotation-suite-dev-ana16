#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the CNV report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import gzip
from os import makedirs
from datetime import date
from os.path import dirname, abspath, basename, join
from glob import glob
from multiprocessing import get_context
import xlsxwriter
import polars as pl
import narwhals.stable.v2 as nw


AVAILABLE_GENOME_VERSIONS = ["GRCh37", "GRCh38"]


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  parser.add_argument("-o", "--output-file", type=str, default="output.xlsx",
                      dest="output_file_name", help="The name of the output file. Overrides project_code argument.")
  parser.add_argument("-p", "--project", dest="project_code", type=str, default="",
                      help="The project code for the report. Used to generate filename with the project_code in it.")
  parser.add_argument("-r", "--reference", "genome_reference", type=str, default="GRCh37",
                      choices=AVAILABLE_GENOME_VERSIONS,
                      help="The Genome Version. Options are 'GRCh37' or 'GRCh38'.")
  parser.add_argument("-b" "--bucket-folder", dest="bucket_folder", default="", type=str,
                      help=("The sarek base folder with all of the other files within it. You can either enter this"
                            + "or pass all the files individually."))
  # These 3 types of files are for the QC metrics sheet
  parser.add_argument("-m" "--mos-files", "mos_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .recal.mosdepth.summary.txt files with mean coverage in them.")
  parser.add_argument("-r" "--cnr-files", "cnr_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .cnr files with normalized log2 ratios in them.")
  parser.add_argument("-s" "--cns-files", "cns_files", default=[], action="extend", nargs="+",
                      help=("A list of paths to .cns files with segment calls in them. Use the most unfiltered ones"
                            + " i.e. the .cns files not .call.cns or .tumor_only.call.cns"))
  # Gene summary sheet
  parser.add_argument("-g" "--gene-files", "gene_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .genemetrics.tsv files with gene summaries in them.")
  # Variant sheet
  parser.add_argument("-v" "--vcf-files", "vcf_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .vcf files with the actual variant calls in them.")
  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if args.bucket_folder:
    if (args.mos_files or args.cnr_files or args.cns_files, args.gene_files, args.vcf_files):
      parser.error("No files should be specified if --bucket-folder is specified.")


def gather_input_files(args: argparse.Namespace) -> tuple[list[str], list[str], list[str], list[str], list[str]]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns tuple[5 lists]: 5 lists of file paths (mos_files, cnr_files, cns_files, gene_files, vcf_files) 
  """
  mos_files = []
  cnr_files = []
  cns_files = []
  gene_files = []
  vcf_files = []
  if not args.mos_files:
    mos_files = glob(f"{args.bucket_folder}/reports/mosdepth/*/*.recal.mosdepth.summary.txt")
  else:
    mos_files = args.mos_files

  if not args.cnr_files:
    cnr_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cnr")
  else:
    cnr_files = args.cnr_files

  if not args.cns_files:
    cns_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cns")
    # Remove call.cns files which are basically duplicates. In this case, we want the ones without call in filenames
    # because they have more segment calls which will lower the median absolute deviation (MAD) scores.
    for one_file in cns_files:
      if "call" in one_file:
        cns_files.remove(one_file)
  else:
    cns_files = args.cns_files

  if not args.gene_files:
    gene_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.genemetrics.tsv")
  else:
    gene_files = args.gene_files

  if not args.vcf_files:
    vcf_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cnvcall.vcf")
  else:
    vcf_files = args.vcf_files
  return (mos_files, cnr_files, cns_files, gene_files, vcf_files)


def derive_chromosomal_arms(genome_version: str) -> nw.LazyFrame:
  """
  Derives the chromosomal arms from USCS cytoband file. Downloaded from:
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
  """
  l_genome = genome_version.lower()
  cytoband_data = nw.from_native(pl.LazyFrame({}))
  with gzip.open(f'{l_genome}CytoBand.txt.gz', 'rb') as cytoband_file:
    cytoband_data = pl.scan_csv(
        cytoband_file,
        separator="\t",
        has_header=False,
        new_columns=["chromosome", "start", "end", "band", "stain"]
    )
    # Aggregate by arm substring
    cytoband_data = nw.from_native(cytoband_data)
    cytoband_data = cytoband_data.with_columns(arm=nw.col("band").str.slice(offset=0, length=1))  # first letter: p/q
    cytoband_data = cytoband_data.drop("band", "stain")
    cytoband_data = cytoband_data.group_by("chromosome", "arm").agg(nw.col("start").min(), nw.col("end").max())
  return cytoband_data.collect()


def organize_qc_files(mos_files: list[str], cnr_files: list[str], cns_files: list[str]) -> nw.LazyFrame:
  """
  Organizes all the files by sample/basename into a DataFrame.
  """
  mos_df = nw.from_native(pl.LazyFrame({"mos_files": mos_files}))
  cnr_df = nw.from_native(pl.LazyFrame({"cnr_files": cnr_files}))
  cns_df = nw.from_native(pl.LazyFrame({"cns_files": cns_files}))
  mos_df = mos_df.with_columns(basename=nw.col(mos_df.col("mos_files").str.replace(".recal.mosdepth.summary.txt", "")))
  cnr_df = cnr_df.with_columns(basename=nw.col(cnr_df.col("cnr_files").str.replace(".cnr", "")))
  cns_df = cns_df.with_columns(basename=nw.col(cns_df.col("cns_files").str.replace(".cns", "")))
  all_qc_files = mos_df.join(cnr_df, "basename", how="full")
  all_qc_files = all_qc_files.join(cnr_df, "basename", how="full")
  return all_qc_files.collect()


def get_qc_metrics(mos_files: list[str], cnr_files: list[str], cns_files: list[str]) -> nw.LazyFrame:
  """
  Gather all the QC Metrics for the CNV report. Some of these will be generated by a docker command `cnvkit.py metrics`
  """
  all_qc_files = organize_qc_files(mos_files, cnr_files, cns_files)
  # qc_metrics = columns(sample_name, MTC, cnvkit.py metrics)
  for 


def graph_qc_metrics():
  """
  Create 2 QC graphs:
  1. Mean target coverage (MTC) histogram
  2. Median absolute deviation (MAD) histogram.
  """


def get_cnv_variants() -> nw.LazyFrame:
  """
  A dump of all the official CNV variants output by cnvkit.
  """


def summarize_gene_level() -> nw.LazyFrame:
  """
  Summarizes all the calls (gains/losses) by gene.
  """


def summarize_chromosome_arm_level() -> nw.LazyFrame:
  """
  Summarizes all the calls (gains/losses) by chromosomal arm.
  """


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  output_file_name = args.output_file_name
  today = date.today().strftime('%Y%m%d')
  if output_file_name == "output.xlsx" and args.project_code:
    project_report_dir = f"{args.project_code}_Reports"
    makedirs(project_report_dir, exist_ok=True)
    output_file_name = join(project_report_dir, f"{today}_{args.project_code.lower()}_cnvkit_report.xlsx")
  chromosomal_arms = derive_chromosomal_arms(args.genome_version)
  (mos_files, cnr_files, cns_files, gene_files, vcf_files) = gather_input_files(args)
  qc_metrics = nw.from_native(pl.LazyFrame({}))
  cnv_variants = nw.from_native(pl.LazyFrame({}))
  gene_summary = nw.from_native(pl.LazyFrame({}))
  arm_summary = nw.from_native(pl.LazyFrame({}))
  if mos_files and cnr_files and cns_files:
    qc_metrics = get_qc_metrics(mos_files, cnr_files, cns_files)
    graph_qc_metrics()
  if vcf_files:
    cnv_variants = get_cnv_variants(vcf_files)
  if gene_files:
    gene_summary = summarize_gene_level(gene_files)
  if cnv_variants:
    arm_summary = summarize_chromosome_arm_level(chromosomal_arms, cnv_variants)


if __name__ == "__main__":
  main()
