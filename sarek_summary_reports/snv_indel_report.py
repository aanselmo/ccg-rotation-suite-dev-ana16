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
import polars as pl
import narwhals.stable.v2 as nw
from lib.snv_indel.snv_indel_sheet import get_variants, VEP_FILE_ENDING
from lib.snv_indel.variant_summary_sheet import get_variant_summary


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  parser.add_argument("-o", "--output-file", type=str, default="output.xlsx",
                      dest="output_file_name", help="The name of the output file. Overrides project_code argument.")
  parser.add_argument("-p", "--project", dest="project_code", type=str, default="",
                      help="The project code for the report. Used to generate filename with the project_code in it.")
  parser.add_argument("-b", "--bucket-folder", dest="bucket_folder", default="", type=str,
                      help=("The sarek base folder with all of the other files within it. You can either enter this"
                            + "or pass all the files individually."))
  parser.add_argument("-v", "--vcf-files", dest="vcf_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .vcf files with the actual variant calls in them.")
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
  vcf_files = gather_input_files(args)
  all_variants = get_variants(vcf_files)
  variant_summary_table = get_variant_summary(all_variants)
