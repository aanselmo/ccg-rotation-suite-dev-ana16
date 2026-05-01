#!/usr/local/miniforge3/bin/python

"""
Use this script to rename the fastqs sent by MBCF for Sarek. The fastqs received from MBCF a lot of times have random
names. This makes sure that the collaborator sample name is in the file name for the fastq.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import argparse
import re
import io
import warnings
from tempfile import TemporaryFile
from datetime import date
from os import makedirs
from os.path import join, basename, abspath, dirname
from concurrent.futures import ProcessPoolExecutor
import narwhals.stable.v2 as nw
from lib.common import append_common_args, DEFAULT_CSV_OUTPUT, BASE_CLOUD_DIR
from lib.lims import get_one_project_metadata, get_indices_of_samples
from lib.google import get_storage_client, get_head_of_gcp_gzip_file, get_files_in_folder, get_blob_from_gs_uri


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
LIMS_SAMPLE_NAMES_COL = "lims_sample_names"
COLLABORATOR_SAMPLE_NAMES_COL = "collaborator_sample_names"


def setup_arguments():
  parser = argparse.ArgumentParser(description=(
      "Use this script to rename the fastqs sent by MBCF for Sarek. The fastqs received from MBCF a lot of times have"
      + "random names. This makes sure that the collaborator sample name is in the file name for the fastq."
  ))
  project_arg = append_common_args(parser)
  project_arg.required = True
  parser.add_argument("-m", "--mbcf-dir", dest="mbcf_dir", type=str, default="", required=True,
                      help="The MBCF Bucket Directory with all the fastqs in them (e.g. gs://mbcgs/260424_AB14103_fastq).")
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


def determine_file_indices(gcp_storage_client, fastq_blob) -> list[str]:
  """
  Find the molecular barcodes used in the fastq.
  :return list[str]: Can be a list of 1 if there is 1 barcode. Usually, however, in CCG we use dual-indices.
      e.g. [TGAGCTAG, AACCGTTC]
  """
  indices_in_file: list[str] = []
  readname_line = get_head_of_gcp_gzip_file(fastq_blob, num_lines=1)
  position_result = re.search(r":[A-Z+]+$", readname_line)  # last section of read name with index/indices in them.
  indices_in_file = readname_line[position_result.start() + 1:].strip()
  indices_in_file = indices_in_file.split("+")
  indices_in_file = [index.replace("N", "[A-Z]") for index in indices_in_file]
  return indices_in_file


def get_all_fastq_indices(gcp_storage_client, mbcf_directory: str) -> nw.DataFrame:
  all_fastqs = get_files_in_folder(gcp_storage_client, mbcf_directory)
  r1_read_regex = re.compile(f".*{ALL_READ_TYPES[0]}_001.fastq.gz")
  # These are parallel arrays that will later be put into a DataFrame
  r1s_in_dir = []
  for fastq in all_fastqs:
    if r1_read_regex.match(fastq.name):
      r1s_in_dir.append(fastq)
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


def get_all_lims_indices(project_code: str) -> nw.DataFrame:
  # These are parallel arrays that will later be put into a DataFrame
  lims_sample_names = []
  collaborator_sample_names = []
  first_barcodes = []
  second_barcodes = []
  project_metadata = get_one_project_metadata(project_code)
  all_samples = project_metadata["samples"]
  for one_sample in all_samples:
    lims_sample_names.append(one_sample["name"])
    collaborator_sample_names.append(one_sample["externalId"])
  first_barcodes = [""] * len(lims_sample_names)
  second_barcodes = [""] * len(lims_sample_names)
  indices_response = get_indices_of_samples(lims_sample_names)
  for one_sample in indices_response:
    one_lims_sample_name = one_sample["name"]
    one_first_sequence = one_sample["libraries"]["index"]["firstSequence"]
    one_second_sequence = one_sample["libraries"]["index"]["secondSequence"]
    row_index = lims_sample_names.index(one_lims_sample_name)
    first_barcodes[row_index] = one_first_sequence
    second_barcodes[row_index] = one_second_sequence
  lims_indices = nw.DataFrame.from_dict({
    LIMS_SAMPLE_NAMES_COL: lims_sample_names,
    COLLABORATOR_SAMPLE_NAMES_COL: collaborator_sample_names,
    FIRST_BARCODE_COL: first_barcodes,
    SECOND_BARCODE_COL: second_barcodes
  }, backend="polars")
  return lims_indices


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
  lims_indices = get_all_lims_indices(project_code)
  # 2c. Outer join in polars?
  mapped_indices = fastq_indices.join(lims_indices, on=[FIRST_BARCODE_COL, SECOND_BARCODE_COL], how="full")
  # 3. create a table of how it did the matching. Output the .csv
  with TemporaryFile(mode='w+') as f:
    mapped_indices.write_csv(f)
    upload_location = get_blob_from_gs_uri(join(BASE_CLOUD_DIR, project_code, output_file_name))
    upload_location.upload_from_filename(f)
  # 4. Copy the fastqs to the correct directory.


if __name__ == "__main__":
  main()