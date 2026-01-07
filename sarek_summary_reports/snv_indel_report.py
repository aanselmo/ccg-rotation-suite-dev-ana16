#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the SNV/Indel report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
from os import makedirs
from os.path import join, basename
from glob import glob
from datetime import date
import asyncio
import polars as pl
import narwhals.stable.v2 as nw
from lib.common_args import append_common_args, append_reference_arg, append_vcf_files_arg
from lib.snv_indel.lego_plot_sheet import graph_lego_plot
from lib.snv_indel.snv_indel_sheet import get_variants, VEP_FILE_ENDING
from lib.snv_indel.variant_summary_sheet import get_variant_summary
from lib.snv_indel.analysis_legend_sheet import get_legend


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  append_common_args(parser)
  append_reference_arg(parser)
  parser.add_argument("-v", "--vcf-files", dest="vcf_files", default=[], action="extend", nargs="+",
                      help=f"A list of paths to {VEP_FILE_ENDING} files with the actual variant calls in them.")
  parser.add_argument("-e", "--bed-file", dest="bed_file", type=str, required=True,
                      help="The targeted bed file for the bait set. If it's Whole genome sequencing, pass 'WGS'.")
  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if args.bucket_folder:
    if (args.vcf_files):
      parser.error("No files should be specified if --bucket-folder is specified.")
  else:
    if not (args.vcf_files):
      parser.error("Please specify a --bucket-folder or individual files (--vcf-files).")


def gather_input_files(args: argparse.Namespace) -> list[str]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns list[str]: A list of strings for vcf_files paths
  """
  vcf_files = []

  if not args.vcf_files:
    vcf_files = glob(f"{args.bucket_folder}/annotation/mutect2/*/*{VEP_FILE_ENDING}")
  else:
    vcf_files = args.vcf_files
  return vcf_files


async def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  output_file_name = args.output_file_name
  today = date.today().strftime('%Y%m%d')
  project_code = "test"
  if output_file_name == "output.xlsx" and args.project_code:
    project_code = args.project_code
    project_report_dir = f"{project_code}_Reports"
    makedirs(project_report_dir, exist_ok=True)
    output_file_name = join(project_report_dir, f"{today}_{project_code.lower()}_cnvkit_report.xlsx")
  vcf_files = gather_input_files(args)
  bed_file = None
  if args.bed_file != "WGS":
    bed_file = args.bed_file
  graph_task = asyncio.create_task(graph_lego_plot(vcf_files, project_code, args.reference_genome, bed_file))
  all_variants = get_variants(vcf_files)
  variant_summary_table = get_variant_summary(all_variants)
  legend = get_legend(args.reference_genome)
  await graph_task


if __name__ == "__main__":
  asyncio.run(main())
