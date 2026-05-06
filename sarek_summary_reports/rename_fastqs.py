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
from tempfile import NamedTemporaryFile
from itertools import cycle
from datetime import date
from os import makedirs
from os.path import join, basename, abspath, dirname
from concurrent.futures import ProcessPoolExecutor
import polars as pl
from lib.common import append_common_args, append_processes_arg, DEFAULT_CSV_OUTPUT, BASE_CLOUD_DIR
from lib.lims import get_one_project_metadata, get_indices_of_samples
from lib.google import (
    get_storage_client, get_head_of_gcp_gzip_file, get_files_in_folder, get_blob_from_gs_uri, copy_gcp_file
)


PROJECT_FASTQS_DIR = join("preprocessing", "fastq")
DEFAULT_FASTQ_MAP_CSV = "mbcf_fastq_map.csv"
# The first read type should be one that all samples have. It should always include the indices in the read names
# of the fastq
ALL_READ_TYPES = ["R1", "R2", "R3", "I1", "I2"]
BEFORE_LANE_START = -21
READ_TYPE_START = -15
READ_TYPE_END = -13
EXTERNAL_NAME_REGEX = r"S[0-9]+"
NUMBER_OF_ROWS = 0
NUMBER_OF_COLUMNS = 1
R1_FASTQ_COL = "r1_fastqs"
FIRST_BARCODE_LIMS_COL = "first_barcode_lims"
SECOND_BARCODE_LIMS_COL = "second_barcode_lims"
FIRST_BARCODE_FASTQ_COL = "first_barcode_fastq"
SECOND_BARCODE_FASTQ_COL = "second_barcode_fastq"
SECOND_BARCODE_REVCOMP_COL = "second_barcode_lims_revcomp"
LIMS_SAMPLE_NAMES_COL = "lims_sample_names"
COLLABORATOR_SAMPLE_NAMES_COL = "collaborator_sample_names"


def setup_arguments():
  parser = argparse.ArgumentParser(description=(
      "Use this script to rename the fastqs sent by MBCF for Sarek. The fastqs received from MBCF a lot of times have"
      + "random names. This makes sure that the collaborator sample name is in the file name for the fastq."
  ))
  project_arg = append_common_args(parser)
  append_processes_arg(parser)
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


def determine_file_indices(fastq_uri: str) -> list[str]:
  """
  Find the molecular barcodes used in the fastq.
  :return list[str]: Can be a list of 1 if there is 1 barcode. Usually, however, in CCG we use dual-indices.
      e.g. [TGAGCTAG, AACCGTTC]
  """
  indices_in_file: list[str] = ["", ""]
  storage_client = get_storage_client()
  readname_line = get_head_of_gcp_gzip_file(storage_client, fastq_uri, num_lines=1)
  position_result = re.search(r":[A-Z+]+$", readname_line)  # last section of read name with index/indices in them.
  if position_result:
    indices_in_file = readname_line[position_result.start() + 1:].strip()
    indices_in_file = indices_in_file.split("+")
    indices_in_file = [index.replace("N", "[A-Z]") for index in indices_in_file]
  return indices_in_file


def get_all_fastq_indices(gcp_storage_client, mbcf_directory: str, max_workers=1) -> pl.DataFrame:
  all_fastqs = get_files_in_folder(gcp_storage_client, mbcf_directory)
  r1_read_regex = re.compile(f".*{ALL_READ_TYPES[0]}_001.fastq.gz")
  r1s_in_dir = list(filter(r1_read_regex.match, all_fastqs))
  # These are parallel arrays that will later be put into a DataFrame
  first_barcodes = []
  second_barcodes = []
  all_indices = []
  with ProcessPoolExecutor(max_workers=max_workers) as determine_indices_pool:
    all_indices = list(determine_indices_pool.map(determine_file_indices, r1s_in_dir))
  if all_indices:
    first_barcodes, second_barcodes = zip(*all_indices)
  r1_fastqs_with_indices = pl.from_dict({
    R1_FASTQ_COL: r1s_in_dir,
    FIRST_BARCODE_FASTQ_COL: first_barcodes,
    SECOND_BARCODE_FASTQ_COL: second_barcodes
  })
  return r1_fastqs_with_indices


def get_all_lims_indices(project_code: str) -> pl.DataFrame:
  # These are parallel arrays that will later be put into a DataFrame
  lims_sample_names = []
  collaborator_sample_names = []
  first_barcodes = []
  second_barcodes = []
  FIRST_ELEMENT = 0
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
    one_first_sequence = one_sample["libraries"][FIRST_ELEMENT]["index"]["firstSequence"][FIRST_ELEMENT]
    one_second_sequence = one_sample["libraries"][FIRST_ELEMENT]["index"]["secondSequence"][FIRST_ELEMENT]
    row_index = lims_sample_names.index(one_lims_sample_name)
    first_barcodes[row_index] = one_first_sequence
    second_barcodes[row_index] = one_second_sequence
  second_revcomp = map(reverse_complement, second_barcodes)
  lims_indices = pl.from_dict({
    LIMS_SAMPLE_NAMES_COL: lims_sample_names,
    COLLABORATOR_SAMPLE_NAMES_COL: collaborator_sample_names,
    FIRST_BARCODE_LIMS_COL: first_barcodes,
    SECOND_BARCODE_LIMS_COL: second_barcodes,
    SECOND_BARCODE_REVCOMP_COL: second_revcomp
  })
  return lims_indices


def copy_one_fastq(all_args: tuple[str, str, str]):
  (fastq, sample_name, destination_folder) = all_args
  gcp_storage_client = get_storage_client()
  external_name_matches = [match for match in re.finditer(EXTERNAL_NAME_REGEX, fastq)]
  tail_filename = fastq[external_name_matches[-1].start():]
  dest_file_name = sample_name + "_" + tail_filename
  destination_uri = join(destination_folder, dest_file_name)
  copy_gcp_file(gcp_storage_client, fastq, destination_uri)


def copy_all_fastqs(mapped_indices: pl.DataFrame, destination_folder: str, max_workers=1):
  mapped_indices = mapped_indices.drop_nulls([R1_FASTQ_COL, COLLABORATOR_SAMPLE_NAMES_COL])
  r1_fastqs = mapped_indices.get_column(R1_FASTQ_COL).to_list()
  sample_names = mapped_indices.get_column(COLLABORATOR_SAMPLE_NAMES_COL).to_list()
  all_copy_args = zip(r1_fastqs, sample_names, cycle([destination_folder]))
  with ProcessPoolExecutor(max_workers=max_workers) as copy_pool:
    copy_pool.map(copy_one_fastq, all_copy_args)


def reverse_complement(barcode: str) -> str:
  """
  A fast way to find the reverse complement of a sequence. See:
  https://stackoverflow.com/questions/8776924/finding-the-complement-of-a-dna-sequence
  """
  REVCOMP_TABLE = str.maketrans("ACGT", "TGCA")
  return barcode.translate(REVCOMP_TABLE)[::-1]


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
  mbcf_dir = args.mbcf_dir
  print("Get indices and collaborator names from LIMS.")
  lims_indices = get_all_lims_indices(project_code)
  print(lims_indices.head())
  print("Get indices in fastqs.")
  fastq_indices = get_all_fastq_indices(gcp_storage_client, mbcf_dir, args.number_processes)
  fastq_indices = fastq_indices.filter(pl.col(FIRST_BARCODE_FASTQ_COL) != "")
  fastq_indices = fastq_indices.filter(pl.col(SECOND_BARCODE_FASTQ_COL) != "")
  print(fastq_indices.head())
  print("Mapping collaborator names in LIMS to fastq indices.")
  mapped_indices = fastq_indices.join_where(
      lims_indices,
      (
          pl.col(FIRST_BARCODE_LIMS_COL).str.contains(pattern=pl.col(FIRST_BARCODE_FASTQ_COL))
          & (pl.col(SECOND_BARCODE_LIMS_COL).str.contains(pattern=pl.col(SECOND_BARCODE_FASTQ_COL))
             | pl.col(SECOND_BARCODE_REVCOMP_COL).str.contains(pattern=pl.col(SECOND_BARCODE_FASTQ_COL)))
      )
  )
  print(mapped_indices.head())
  print("Writing mapping to output file.")
  with NamedTemporaryFile(mode='w+') as map_temp_file:
    mapped_indices.write_csv(map_temp_file.name)
    upload_location = get_blob_from_gs_uri(gcp_storage_client, join(BASE_CLOUD_DIR, project_code, output_file_name))
    upload_location.upload_from_filename(map_temp_file.name)
    breakpoint()
  print("Copy the fastqs to the correct directory.")
  detination_folder = join(BASE_CLOUD_DIR, project_code, PROJECT_FASTQS_DIR)
  # copy_all_fastqs(mapped_indices, detination_folder, args.number_processes)


if __name__ == "__main__":
  main()