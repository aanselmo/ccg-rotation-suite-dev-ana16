#!/data/python-venv/for_ccg_rotation_suite/bin/python

"""
Use this script to generate the aggregate report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie, Neil A. Patel"]

import argparse
import re
import narwhals as nw
from glob import glob
from os.path import basename
from lib.common import (append_common_args, append_bucket_arg, get_output_file_name,
                        DEFAULT_HTML_OUTPUT, DEFAULT_EXCEL_OUTPUT, BASE_CLOUD_DIR)
from lib.aggregate.create_metrics_table import create_metrics_table
from lib.aggregate.create_concordance_matrix import *
from lib.aggregate.create_html_report import create_html_report
from lib.aggregate.write_aggregate_xlsx import write_xlsx

from lib.setup.google import get_storage_client, get_files_in_folder
from lib.aggregate.get_google_sarek import get_google_picard, get_google_somalier, get_google_duplicates
from lib.setup.lims_df import get_all_lims_data

DEFAULT_AGG_HTML_REPORT_FILE_ENDING = "aggregation_metrics_report.html"
DEFAULT_AGG_EXCEL_REPORT_FILE_ENDING = "aggregation_metrics_report.xlsx"


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  append_common_args(parser, default_output_file=DEFAULT_HTML_OUTPUT)
  append_bucket_arg(parser)
  parser.add_argument("-d", "--duplicate-marking-uri", dest="duplicate_marking_uri", default=[], action="extend",
                      nargs="+", help="A list of gs:// paths to .md.cram.metrics files PERCENT_DUPLICATION in them.")
  parser.add_argument("-t", "--somalier_pairs_uri", dest="somalier_pairs_uri", default="", type=str,
                      help="One gs:// uri to a somalier.pairs.tsv file.")
  parser.add_argument("-i", "--picard_uri", dest="picard_uri", default="", type=str,
                      help="A list of gs:// paths to /reports/picard/*hsmetrics*.coverage_metrics files.")

  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if not args.project_code:
    if not (args.duplicate_marking_uri or args.somalier_pairs_uri or args.picard_uri):
      parser.error("Please specify a --project or individual files (--duplicate-marking-uri, --somalier_pairs_uri, "
                   + "--picard_uri.")


def gather_input_files(args: argparse.Namespace, storage_client) -> tuple[list[str], list[str], list[str]]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns tuple[3 lists]: 3 lists of google cloud uris
  """
  duplicate_marking_uri = []
  somalier_pairs_uri = ""
  picard_uri = []

  project_uri = f"{BASE_CLOUD_DIR}/{args.project_code}/"
  google_files = get_files_in_folder(storage_client, project_uri)

  if not args.duplicate_marking_uri:
    r_duplicate = re.compile(".*reports/markduplicates/.*/.*.md.cram.metrics$")
    duplicate_marking_uri = list(filter(r_duplicate.match, google_files))
  else:
    duplicate_marking_files = args.duplicate_marking_files

  if not args.somalier_pairs_uri:
    r_somalier = re.compile(".*reports/somalier/somalier.pairs.tsv$")
    somalier_pairs_uri = list(filter(r_somalier.match, google_files))[0]
  else:
    somalier_pairs_uri = args.somalier_pairs_uri

  if not args.picard_uri:
    r_picard = re.compile(".*reports/picard/.*hsmetrics.*.coverage_metrics$")
    picard_uri = list(filter(r_picard.match, google_files))
  else:
    picard_uri = args.picard_uri

  return (duplicate_marking_uri, somalier_pairs_uri, picard_uri)


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (project_code, html_output_file_name) = get_output_file_name(
      args, DEFAULT_AGG_HTML_REPORT_FILE_ENDING, DEFAULT_HTML_OUTPUT)

  args_xlsx = args
  args_xlsx.default_output_file = 'output.xlsx'
  args_xlsx.output_file = 'output.xlsx'
  args_xlsx.output_file_name = 'output.xlsx'
  (project_code, xlsx_output_file_name) = get_output_file_name(
      args, DEFAULT_AGG_EXCEL_REPORT_FILE_ENDING, DEFAULT_EXCEL_OUTPUT)

  storage_client = get_storage_client()
  (duplicate_marking_uri, somalier_pairs_uri, picard_uri) = gather_input_files(args, storage_client)

  data_picard = get_google_picard(storage_client, picard_uri)
  data_lims = nw.from_native(get_all_lims_data(project_code))
  data_duplicates = get_google_duplicates(storage_client, duplicate_marking_uri)
  (somalier_pairs, somalier_mismatch) = get_google_somalier(storage_client, somalier_pairs_uri)

  (metrics_table, summary_table) = create_metrics_table(data_lims, data_picard, data_duplicates)
  concordance_table = create_concordance_matrix(somalier_pairs)

  print(xlsx_output_file_name)
  write_xlsx(xlsx_output_file_name, metrics_table, summary_table, concordance_table, somalier_mismatch)
  # create_html_report(summary_table.to_pandas(), tsv_somalier_pairs, project_code, html_output_file_name)


if __name__ == "__main__":
  main()
