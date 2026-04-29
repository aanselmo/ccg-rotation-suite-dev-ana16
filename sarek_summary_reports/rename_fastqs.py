#!/usr/local/miniforge3/bin/python

"""
Use this script to rename the fastqs sent by MBCF for Sarek. The fastqs received from MBCF a lot of times have random
names. This makes sure that the collaborator sample name is in the file name for the fastq.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import argparse
import re
import gzip
import warnings
import sys
from glob import glob
from datetime import date
from os import makedirs
from os.path import join, basename, abspath, dirname
from concurrent.futures import ProcessPoolExecutor
from google.cloud import storage
from lib.common import (
    append_common_args, append_reference_arg, DEFAULT_CSV_OUTPUT
)


git_repo_dir = abspath(dirname(dirname(__file__)))
SERVICE_ACCOUNT = join(git_repo_dir, "secrets", "rnaseq_ccgd-174313-98def6888a16.json")
PROJECT_FASTQS_DIR = join("preprocessing", "fastq")
DEFAULT_FASTQ_MAP_CSV = "mbcf_fastq_map.csv"
# The first read type should be one that all samples have. It should always include the indices in the read names
# of the fastq
ALL_READ_TYPES = ["R1", "R2", "R3", "I1", "I2"]
BEFORE_LANE_START = -21
READ_TYPE_START = -15
READ_TYPE_END = -13
EXTERNAL_NAME_REGEX = re.compile("S[0-9]+")
NUMBER_OF_ROWS = 0
NUMBER_OF_COLUMNS = 1
DOMAIN = "https://ccgd2.dfci.harvard.edu/minilims/"
# project.php is a POST request with a secretKey and project in the body of the request (as multipart form data).
PROJECT_API_URL = DOMAIN + "api/get/10x-samplesheet/from/project.php"
# getIndicesFromSamples.php is a GET request with sample query params in the url
# Example: https://ccgd2.dfci.harvard.edu/minilims/api/getIndicesFromSamples.php?samples=["S_035150"]
INDICES_API_URL = DOMAIN + "api/getIndicesFromSamples.php"


def setup_arguments():
  parser = argparse.ArgumentParser(description=(
      "Use this script to rename the fastqs sent by MBCF for Sarek. The fastqs received from MBCF a lot of times have"
      + "random names. This makes sure that the collaborator sample name is in the file name for the fastq."
  ))
  project_arg = append_common_args(parser)
  project_arg.required = True
  parser.add_argument("-m", "--mbcf-dir", dest="mbcf_dir", type=str, default="", required=True,
                      help="The MBCF Bucket Directory with all the fastqs in them.")
  return parser.parse_args()


def get_output_file_name(args, file_ending: str, default_file_name: str) -> tuple[str, str]:
  output_file_name = args.output_file_name
  today = date.today().strftime('%Y%m%d')
  project_code = "test"
  if output_file_name == default_file_name and args.project_code:
    project_code = args.project_code
    makedirs(PROJECT_FASTQS_DIR, exist_ok=True)
    output_file_name = join(PROJECT_FASTQS_DIR, f"{today}_{project_code.lower()}_{file_ending}")
  return (project_code, output_file_name)



def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (project_code, output_file_name) = get_output_file_name(
      args,
      DEFAULT_FASTQ_MAP_CSV,
      DEFAULT_CSV_OUTPUT
  )
  # 1. Get the folder of MBCF Fastqs
  # 2. Bin the fastqs.
  # 3. Copy the fastqs to the correct directory.
  # 4. create a table of how it did the matching. Output the .csv


if __name__ == "__main__":
  main()