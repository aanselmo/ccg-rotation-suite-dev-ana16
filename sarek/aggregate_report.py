#!/data/python-venv/for_ccg_rotation_suite/bin/python

"""
Use this script to generate the aggregate report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie, Neil A. Patel"]

import argparse
import numpy as np
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
from lib.aggregate.somalier_processing import append_somalier_sample_info_to_mismatch

from lib.setup.google import get_storage_client, get_files_in_folder
from lib.aggregate.get_google_sarek import (get_google_picard, get_google_somalier_pairs, get_google_somalier_samples,
                                            get_google_duplicates, get_google_lanes)
from lib.setup.lims_df import get_all_lims_data, get_illumina_run

DEFAULT_AGG_HTML_REPORT_FILE_ENDING = "aggregation_metrics_report.html"
DEFAULT_AGG_EXCEL_REPORT_FILE_ENDING = "aggregation_metrics_report.xlsx"


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  append_common_args(parser, default_output_file=DEFAULT_HTML_OUTPUT)
  append_bucket_arg(parser)
  parser.add_argument("-d", "--duplicate-marking-uri", dest="duplicate_marking_uri", default=[], action="extend",
                      help="A list of gs:// paths to .md.cram.metrics files PERCENT_DUPLICATION in them.")
  parser.add_argument("-t", "--somalier-pairs-uri", dest="somalier_pairs_uri", default="", type=str,
                      help="One gs:// uri to a somalier.pairs.tsv file.")
  parser.add_argument("-v", "--somalier-samples-uri", dest="somalier_samples_uri", default="", type=str,
                      help="One gs:// uri to a somalier.samples.tsv file.")
  parser.add_argument("-i", "--picard-uri", dest="picard_uri", default=[], type=str,
                      help="A list of gs:// paths to /reports/picard/*hsmetrics*.coverage_metrics files.")
  parser.add_argument("-s", "--samplesheet-uri", dest="samplesheet_uri", default="", type=str,
                      help="One gs:// path to a samplesheet.csv file.")

  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if not args.project_code:
    if not (args.duplicate_marking_uri or args.somalier_pairs_uri or args.somalier_samplse_uri or args.picard_uri or 
            args.samplesheet_uri):
      parser.error("Please specify a --project or individual files (--duplicate-marking-uri, --somalier_pairs_uri, "
                   + "--somalier_samples_uri, --picard_uri, --samplesheet-uri.")


def gather_input_files(args: argparse.Namespace, storage_client) -> tuple[list[str], list[str], list[str], list[str], list[str]]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns tuple[4 lists]: 4 lists of google cloud uris
  """
  duplicate_marking_uri = []
  somalier_pairs_uri = ""
  somalier_samples_uri = ""
  picard_uri = []

  project_uri = f"{BASE_CLOUD_DIR}/{args.project_code}/"
  google_files = get_files_in_folder(storage_client, project_uri)

  if not args.duplicate_marking_uri:
    r_duplicate = re.compile(".*reports/markduplicates/.*/.*.md.cram.metrics$")
    duplicate_marking_uri = list(filter(r_duplicate.match, google_files))
  else:
    duplicate_marking_files = args.duplicate_marking_files

  if not args.somalier_pairs_uri:
    r_somalier_pairs = re.compile(".*reports/somalier/somalier.pairs.tsv$")
    somalier_pairs_uri = list(filter(r_somalier_pairs.match, google_files))[0]
  else:
    somalier_pairs_uri = args.somalier_pairs_uri

  if not args.somalier_samples_uri:
    r_somalier_samples = re.compile(".*reports/somalier/somalier.samples.tsv$")
    somalier_samples_uri = list(filter(r_somalier_samples.match, google_files))[0]
  else: 
    somalier_samples_uri = args.somalier_samples_uri

  if not args.picard_uri:
    r_picard = re.compile(".*reports/picard/.*hsmetrics.*.coverage_metrics$")
    picard_uri = list(filter(r_picard.match, google_files))
  else:
    picard_uri = args.picard_uri

  if not args.samplesheet_uri: 
    r_samplesheet = re.compile(".*/csv/.*samplesheet.csv$")
    samplesheet_uri = list(filter(r_samplesheet.match, google_files))[0]


  return (duplicate_marking_uri, somalier_pairs_uri, somalier_samples_uri, picard_uri, samplesheet_uri)


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (project_code, html_output_file_name) = get_output_file_name(
      args, DEFAULT_AGG_HTML_REPORT_FILE_ENDING, DEFAULT_HTML_OUTPUT)

  run_name = get_illumina_run(project_code)

  if run_name:
    clean_run_name = re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F]", "-", run_name)
    excel_file_ending = f"{clean_run_name}_{DEFAULT_AGG_EXCEL_REPORT_FILE_ENDING}"
    data_lims = nw.from_native(get_all_lims_data(project_code))
  else: 
    excel_file_ending = DEFAULT_AGG_EXCEL_REPORT_FILE_ENDING
    data_lims = []

  args_xlsx = args
  args_xlsx.default_output_file = 'output.xlsx'
  args_xlsx.output_file = 'output.xlsx'
  args_xlsx.output_file_name = 'output.xlsx'
  (project_code, xlsx_output_file_name) = get_output_file_name(
      args, excel_file_ending, DEFAULT_EXCEL_OUTPUT)

  storage_client = get_storage_client()
  (duplicate_marking_uri, somalier_pairs_uri, somalier_samples_uri, picard_uri, samplesheet_uri) = gather_input_files(args, storage_client)

  data_picard = get_google_picard(storage_client, picard_uri)
  
  data_duplicates = get_google_duplicates(storage_client, duplicate_marking_uri)
  (somalier_pairs, somalier_mismatch) = get_google_somalier_pairs(storage_client, somalier_pairs_uri, data_lims)
  somalier_samples =  get_google_somalier_samples(storage_client, somalier_samples_uri, data_lims)
  somalier_mismatch = append_somalier_sample_info_to_mismatch(somalier_mismatch, somalier_samples)
  data_num_lanes = get_google_lanes(storage_client, samplesheet_uri)
  (metrics_table, summary_table) = create_metrics_table(data_lims, data_picard, data_duplicates, data_num_lanes)
  concordance_table = create_concordance_matrix(somalier_pairs)

  num_somalier_sites = np.unique(somalier_samples['n_possible_sites'].to_list()).tolist()
  if len(num_somalier_sites) == 1:
    num_somalier_sites = num_somalier_sites[0]

  write_xlsx(xlsx_output_file_name, metrics_table, summary_table, concordance_table, somalier_mismatch, num_somalier_sites)
  # create_html_report(summary_table.to_pandas(), tsv_somalier_pairs, project_code, html_output_file_name)


if __name__ == "__main__":
  main()
