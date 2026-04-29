#!/usr/local/miniforge3/bin/python

"""
Use this script to rename the fastqs sent by MBCF for Sarek.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import argparse
import os
import re
import gzip
import warnings
import sys
from glob import glob
from os.path import join, basename, abspath, dirname
from concurrent.futures import ProcessPoolExecutor
from google.cloud import storage
from lib.common import (
    append_common_args, append_reference_arg, get_output_file_name, DEFAULT_EXCEL_OUTPUT
)


git_repo_dir = abspath(dirname(dirname(__file__)))
SERVICE_ACCOUNT = join(git_repo_dir, "secrets", "rnaseq_ccgd-174313-98def6888a16.json")
# The first read type should be one that all samples have. It should always include the indices in the read names
# of the fastq
ALL_READ_TYPES = ["R1", "R2", "R3", "I1", "I2"]
BEFORE_LANE_START = -21
READ_TYPE_START = -15
READ_TYPE_END = -13
EXTERNAL_NAME_REGEX = re.compile("S[0-9]+")
NUMBER_OF_ROWS = 0
NUMBER_OF_COLUMNS = 1


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
  # 1. Get the folder of MBCF Fastqs
  # 2. Bin the fastqs.
  # 3. Copy the fastqs to the correct directory.
  # 4. create a table of how it did the matching. Output the .csv


if __name__ == "__main__":
  main()