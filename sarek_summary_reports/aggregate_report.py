#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the aggregate report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
from glob import glob
from sarek_summary_reports.lib.common import append_common_args, get_output_file_name, DEFAULT_HTML_OUTPUT
from lib.aggregate.create_summary_table import create_summary_table
from lib.aggregate.create_html_report import create_html_report

DEFAULT_AGG_REPORT_FILE_ENDING = "aggregation_metrics_report.html"


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  append_common_args(parser, default_output_file=DEFAULT_HTML_OUTPUT)
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


def gather_input_files(args: argparse.Namespace) -> tuple[list[str], list[str], list[str], list[str], list[str]]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns tuple[5 lists]: 5 lists of file paths (mos_files, cnr_files, cns_files, gene_files, vcf_files) 
  """
  cumcov_dist_file = ""
  duplicate_marking_files = []
  mos_files = []
  samtool_files = []

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

  return (cumcov_dist_file, duplicate_marking_files, mos_files, samtool_files, tsv_somalier_pairs)


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (project_code, output_file_name) = get_output_file_name(args, DEFAULT_AGG_REPORT_FILE_ENDING, DEFAULT_HTML_OUTPUT)
  (cumcov_dist_file, duplicate_marking_files, mos_files, samtool_files, tsv_somalier_pairs) = gather_input_files(args)
  summary_table = create_summary_table(cumcov_dist_file, duplicate_marking_files, mos_files, samtool_files)
  create_html_report(summary_table, tsv_somalier_pairs, project_code, output_file_name)


if __name__ == "__main__":
  main()
