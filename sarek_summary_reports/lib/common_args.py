#!/usr/local/miniforge3/bin/python

"""
This script contains common arguments used by most scripts in the sarek summary reports.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from argparse import ArgumentParser

AVAILABLE_GENOME_VERSIONS = ["GRCh37", "GRCh38"]


def append_common_args(parser: ArgumentParser):
  parser.add_argument("-o", "--output-file", type=str, default="output.xlsx",
                      dest="output_file_name", help="The name of the output file. Overrides project_code argument.")
  parser.add_argument("-p", "--project", dest="project_code", type=str, default="",
                      help="The project code for the report. Used to generate filename with the project_code in it.")
  parser.add_argument("-b", "--bucket-folder", dest="bucket_folder", default="", type=str,
                      help=("The sarek base folder with all of the other files within it. You can either enter this"
                            + "or pass all the files individually."))


def append_reference_arg(parser: ArgumentParser):
  parser.add_argument("-r", "--reference", dest="genome_reference", type=str, default="GRCh37",
                      choices=AVAILABLE_GENOME_VERSIONS,
                      help="The Genome Version. Options are 'GRCh37' or 'GRCh38'.")
