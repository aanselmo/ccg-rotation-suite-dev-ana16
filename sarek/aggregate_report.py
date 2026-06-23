#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the aggregate report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import polars as pl
import narwhals as nw
import re
from glob import glob
from lib.common import (append_common_args, append_bucket_arg, get_output_file_name, 
                        DEFAULT_HTML_OUTPUT, DEFAULT_EXCEL_OUTPUT)
from lib.aggregate.create_metrics_table import create_metrics_table
from lib.aggregate.create_concordance_matrix import *
from lib.aggregate.create_html_report import create_html_report, parse_somalier_pairs_tsv
from lib.aggregate.write_aggregate_xlsx import write_xlsx

DEFAULT_AGG_HTML_REPORT_FILE_ENDING = "aggregation_metrics_report.html"
DEFAULT_AGG_EXCEL_REPORT_FILE_ENDING = "aggregation_metrics_report.xlsx"


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  append_common_args(parser, default_output_file=DEFAULT_HTML_OUTPUT)
  append_bucket_arg(parser)
  parser.add_argument("-c", "--cumcov-dist-file", dest="cumcov_dist_file", type=str,
                      nargs="+", help="One path to a mosdepth-cumcoverage-dist-id.txt file.")
  parser.add_argument("-d", "--duplicate-marking-files", dest="duplicate_marking_files", default=[], action="extend",
                      nargs="+", help="A list of paths to .md.cram.metrics files PERCENT_DUPLICATION in them.")
  parser.add_argument("-m", "--mos-files", dest="mos_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .recal.mosdepth.summary.txt files with mean coverage in them.")
  parser.add_argument("-s", "--samtool-files", dest="samtool_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .md.cram.stats files with total reads in them.")
  parser.add_argument("-t", "--tsv-somalier-pairs", dest="tsv_somalier_pairs", default="", type=str,
                      help="A list of paths to .md.cram.stats files with total reads in them.")
  parser.add_argument("-e", "--hsmetrics-file", dest="hsmetrics_file", default="", type=str,
                      nargs="+", help="One path to an hsmetrics_report.csv file.")

  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if args.bucket_folder:
    if (args.cumcov_dist_file or args.duplicate_marking_files or args.mos_files):
      parser.error("No files should be specified if --bucket-folder is specified.")
  else:
    if not (args.cumcov_dist_file or args.duplicate_marking_files or args.mos_files):
      parser.error("Please specify a --bucket-folder or individual files (--cumcov-dist-file, "
                   + "--duplicate-marking-files, --mos-files, --gene-files, --vcf-files).")


def gather_input_files(args: argparse.Namespace) -> tuple[list[str], list[str], list[str], list[str], list[str], list[str]]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns tuple[6 lists]: 6 lists of file paths (mos_files, cnr_files, cns_files, gene_files, vcf_files, hsmetrics_file) 
  """
  cumcov_dist_file = ""
  duplicate_marking_files = []
  mos_files = []
  samtool_files = []
  tsv_somalier_pairs = ""
  hsmetrics_file = ""

  if not args.cumcov_dist_file:
    cumcov_dist_file = f"{args.bucket_folder}/multiqc/multiqc_data/mosdepth-cumcoverage-dist-id.txt"
  else:
    cumcov_dist_file = args.cumcov_dist_file

  if not args.duplicate_marking_files:
    duplicate_marking_files = glob(f"{args.bucket_folder}/reports/markduplicates/*/*.md.cram.metrics")
  else:
    duplicate_marking_files = args.duplicate_marking_files

  if not args.mos_files:
    mos_files = glob(f"{args.bucket_folder}/reports/mosdepth/*/*.recal.mosdepth.summary.txt")
  else:
    mos_files = args.mos_files

  if not args.samtool_files:
    samtool_files = glob(f"{args.bucket_folder}/reports/samtools/*/*.md.cram.stats")
  else:
    samtool_files = args.samtool_files

  if not args.tsv_somalier_pairs:
    tsv_somalier_pairs = glob(f"{args.bucket_folder}/reports/somalier/*.pairs.tsv")[0]
  else:
    tsv_somalier_pairs = args.tsv_somalier_pairs

  if not args.hsmetrics_file:
    hsmetrics_file = resolve_hsmetrics_file(args)
  else:
    hsmetrics_file = args.hsmetrics_file

  return (cumcov_dist_file, duplicate_marking_files, mos_files, samtool_files, tsv_somalier_pairs, hsmetrics_file)


def resolve_hsmetrics_file(args: argparse.Namespace):
  """
  Identify proper hsmetrics file to use. If there is more than one *_hsmetrics_reports.csv in the reports folder, 
  then choose the one with the most recent date.  
  :param args argparse.Namespace: All arguments for the script
  :returns string: Path to the most recent hsmetrics file
  """
  
  reports_folder_pattern = f"{args.bucket_folder}/reports/ccg_*_reports"
  reports_folders = glob(reports_folder_pattern)
  if len(reports_folders) != 1:
      raise ValueError(f"expected exactly 1 reports folder, found {len(reports_folders)}")
  reports_folder = reports_folders[0]

  hsmetrics_pattern = f"{reports_folder}/*_hsmetrics_report.csv"
  hsmetrics_matches = glob(hsmetrics_pattern)
  if not hsmetrics_matches:
      raise FileNotFoundError(f"No file matching pattern: {hsmetrics_pattern}")

  # Choose match with latest date

  hsmetrics_dates = [(re.search("\\d{8}", path).group(0), path) for path in hsmetrics_matches]


  latest_path = max(hsmetrics_dates)[1]
  return(latest_path)


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (project_code, html_output_file_name) = get_output_file_name(args, DEFAULT_AGG_HTML_REPORT_FILE_ENDING, 
                                                               DEFAULT_HTML_OUTPUT)
  
  args_xlsx = args
  args_xlsx.default_output_file = 'output.xlsx'
  args_xlsx.output_file = 'output.xlsx'
  args_xlsx.output_file_name = 'output.xlsx'
  (project_code, xlsx_output_file_name) = get_output_file_name(args, DEFAULT_AGG_EXCEL_REPORT_FILE_ENDING,
                                                               DEFAULT_EXCEL_OUTPUT)

  (cumcov_dist_file, duplicate_marking_files, mos_files, 
   samtool_files, tsv_somalier_pairs, hsmetrics_file) = gather_input_files(args)

  (metrics_table, summary_table) = create_metrics_table(hsmetrics_file)
  concordance_table = create_concordance_matrix(tsv_somalier_pairs)


  somalier_pairs = parse_somalier_pairs_tsv(tsv_somalier_pairs)
  pl_somalier_pairs = (
    nw.to_native(somalier_pairs)
    .select('#sample_a', 'sample_b', 'Match Status', 'relatedness', 'ibs0', 'ibs2', 'expected_relatedness')
    .rename({"#sample_a": "Sample1", "sample_b": "Sample2", "Match Status": "Type", 
             "relatedness": "Relatedness", "expected_relatedness": "Expected Relatedness"}))
  pl_somalier_mismatch = pl_somalier_pairs.filter(pl.col("Type") != "OK")
  mismatch_table = nw.from_native(pl_somalier_mismatch)

  print(xlsx_output_file_name)
  write_xlsx(xlsx_output_file_name, metrics_table, summary_table, concordance_table, mismatch_table)
  create_html_report(summary_table.to_pandas(), tsv_somalier_pairs, project_code, html_output_file_name)


if __name__ == "__main__":
  main()
