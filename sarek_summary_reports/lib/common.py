#!/usr/local/miniforge3/bin/python

"""
This script contains common arguments used by most scripts in the sarek summary reports.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from argparse import ArgumentParser
from os import makedirs
from datetime import date
from os.path import join

AVAILABLE_GENOME_VERSIONS = ["GRCh37", "GRCh38"]
DEFAULT_HTML_OUTPUT = "output.html"
DEFAULT_EXCEL_OUTPUT = "output.xlsx"


def append_common_args(parser: ArgumentParser, default_output_file=DEFAULT_EXCEL_OUTPUT):
  parser.add_argument("-o", "--output-file", type=str, default=default_output_file,
                      dest="output_file_name", help="The name of the output file. Overrides project_code argument.")
  parser.add_argument("-p", "--project", dest="project_code", type=str, default="",
                      help="The project code for the report. Used to generate filename with the project_code in it.")
  parser.add_argument("-b", "--bucket-folder", dest="bucket_folder", default="", type=str,
                      help=("The sarek base folder with all of the other files within it. You can either enter this"
                            + "or pass all the files individually."))


def append_reference_arg(parser: ArgumentParser):
  parser.add_argument("-r", "--reference", dest="reference_genome", type=str, default="GRCh37",
                      choices=AVAILABLE_GENOME_VERSIONS,
                      help="The Genome Version. Options are 'GRCh37' or 'GRCh38'.")


def get_output_file_name(args, file_ending: str, default_file_name: str) -> tuple[str, str]:
  project_report_dir = join("reports", "ccg_reports")
  output_file_name = args.output_file_name
  today = date.today().strftime('%Y%m%d')
  project_code = "test"
  if output_file_name == default_file_name and args.project_code:
    project_code = args.project_code
    project_report_dir = join("reports", f"ccg_{project_code.lower()}_reports")
    makedirs(project_report_dir, exist_ok=True)
    output_file_name = join(project_report_dir, f"{today}_{project_code.lower()}_{file_ending}")
  return (project_code, output_file_name)
