#!/data/python-venv/for_ccg_rotation_suite/bin/python

"""
Use this script to generate a samplesheet for Sarek.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import json
import csv
import requests
import argparse
import os
import sys
import tempfile
from lib.common import (
    append_common_args, append_reference_arg, get_output_file_name, DEFAULT_EXCEL_OUTPUT
)


DEFAULT_SAREK_SAMPLESHEET_FILE_ENDING = "sarek_samplesheet.csv"
SECRET_KEY = ""
DOMAIN = "https://ccgd2.dfci.harvard.edu/minilims/"
FLOWCELL_API_URL = DOMAIN + "api/get/10x-samplesheet/from/flowcell.php"
PROJECT_API_URL = DOMAIN + "api/get/10x-samplesheet/from/project.php"


def setup_arguments():
  parser = argparse.ArgumentParser(description="")
  project_arg = append_common_args(parser)
  project_arg.required = True
  return parser.parse_args()

def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (project_code, output_file_name) = get_output_file_name(
      args,
      DEFAULT_SAREK_SAMPLESHEET_FILE_ENDING,
      DEFAULT_EXCEL_OUTPUT
  )


if __name__ == "__main__":
  main()
