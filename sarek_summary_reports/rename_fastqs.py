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
import io
import warnings
import sys
from glob import glob
from datetime import date
from os import makedirs
from os.path import join, basename, abspath, dirname
from concurrent.futures import ProcessPoolExecutor
from google.cloud import storage
import polars as pl
import narwhals.stable.v2 as nw
from lib.common import append_common_args, DEFAULT_CSV_OUTPUT
from lib.lims import get_one_project_metadata, get_indices_of_samples
from lib.google import get_storage_client, get_head_of_gcp_gzip_file, get_files_in_folder


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
R1_FASTQ_COL = "r1_fastqs"
FIRST_BARCODE_COL = "first_barcode"
SECOND_BARCODE_COL = "second_barcode"


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


def determine_file_indices(gcp_storage_client, fastq_file_path: str) -> list[str]:
  """
  Find the molecular barcodes used in the fastq.
  :return list[str]: Can be a list of 1 if there is 1 barcode. Usually, however, in CCG we use dual-indices.
      e.g. [TGAGCTAG, AACCGTTC]
  """
  indices_in_file: list[str] = []
  with gzip.open(fastq_file_path, 'rt') as fastq_file:
    readname_line = fastq_file.readline()
    position_result = re.search(r":[A-Z+]+$", readname_line)  # last section of read name with index/indices in them.
    indices_in_file = readname_line[position_result.start() + 1:].strip()
    indices_in_file = indices_in_file.split("+")
    indices_in_file = [index.replace("N", "[A-Z]") for index in indices_in_file]
  return indices_in_file


def get_all_fastq_indices(gcp_storage_client, mbcf_directory: str) -> nw.DataFrame:
  all_fastqs = get_files_in_folder(gcp_storage_client, mbcf_directory)
  r1_read_regex = re.compile(f".*{ALL_READ_TYPES[0]}_001.fastq.gz")
  r1s_in_dir = list(filter(r1_read_regex.match, all_fastqs))
  first_barcodes = []
  second_barcodes = []
  for r1_file in r1s_in_dir:
    indices_in_file = determine_file_indices(gcp_storage_client, r1_file)
    first_barcodes.append(indices_in_file[0])
    if len(indices_in_file) > 0:
      second_barcodes.append(indices_in_file[1])
    else:
      second_barcodes.append("")
  r1_fastqs_with_indices = nw.DataFrame.from_dict({
    R1_FASTQ_COL: r1s_in_dir,
    FIRST_BARCODE_COL: first_barcodes,
    SECOND_BARCODE_COL: second_barcodes
  }, backend="polars")
  return r1_fastqs_with_indices


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
  gcp_storage_client = get_storage_client()
  # 1. Get the folder of MBCF Fastqs
  mbcf_dir = args.mbcf_dir
  # 2a. Get indices in fastqs.
  fastq_indices = get_all_fastq_indices(gcp_storage_client, mbcf_dir)
  # 2b. Get indices in LIMS.
  # 2c. Outer join in polars?
  # 3. Copy the fastqs to the correct directory.
  # 4. create a table of how it did the matching. Output the .csv


if __name__ == "__main__":
  main()