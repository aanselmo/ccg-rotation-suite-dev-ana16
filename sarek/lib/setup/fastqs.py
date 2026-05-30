#!/usr/bin/env python

"""
Part of the setup script that interacts with fastqs.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import re
from os.path import join
from concurrent.futures import ProcessPoolExecutor, as_completed
import polars as pl
from enum import StrEnum
from pip._vendor.rich.progress import Progress
from .google import (
    get_storage_client, get_head_of_gcp_gzip_file, get_files_in_folder, get_blob_from_gs_uri, copy_gcp_file
)
from ..common import BASE_CLOUD_DIR
from .lims_df import COLLABORATOR_SAMPLE_NAMES_COL


PROJECT_FASTQS_DIR = join("preprocessing", "fastq")
class AllReadTypes(StrEnum):
  # R1 contains the indices and read names
  R1 = "R1"
  R2 = "R2"
  R3 = "R3"
  I1 = "I1"
  I2 = "I2"
START_LANE_NUMBER = -20
READ_TYPE_START = -15
READ_TYPE_END = -13
EXTERNAL_NAME_REGEX = r"S[0-9]+"
FASTQ_COL = "fastqs"
DESTINATION_URI_COL = "destination_uris"
FIRST_BARCODE_FASTQ_COL = "first_barcode_fastq"
SECOND_BARCODE_FASTQ_COL = "second_barcode_fastq"


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
  r1_read_regex = re.compile(f".*{AllReadTypes.R1.value}_001.fastq.gz")
  r1s_in_dir = list(filter(r1_read_regex.match, all_fastqs))
  others_in_dir = list(filter(lambda x: not r1_read_regex.match(x), all_fastqs))
  # These are parallel arrays that will later be put into a DataFrame
  first_barcodes = []
  second_barcodes = []
  all_indices = []
  with ProcessPoolExecutor(max_workers=max_workers) as determine_indices_pool:
    all_indices = list(determine_indices_pool.map(determine_file_indices, r1s_in_dir))
  if all_indices:
    first_barcodes, second_barcodes = zip(*all_indices)
  r1_fastqs_with_indices = pl.from_dict({
      FASTQ_COL: r1s_in_dir,
      FIRST_BARCODE_FASTQ_COL: first_barcodes,
      SECOND_BARCODE_FASTQ_COL: second_barcodes
  }, schema={
      FASTQ_COL: pl.String,
      FIRST_BARCODE_FASTQ_COL: pl.String,
      SECOND_BARCODE_FASTQ_COL: pl.String   
  })
  return (r1_fastqs_with_indices, others_in_dir)


def copy_one_fastq(all_args: tuple[str, str]):
  (fastq, destination_uri) = all_args
  gcp_storage_client = get_storage_client()
  copy_gcp_file(gcp_storage_client, fastq, destination_uri)


def copy_all_fastqs(mapped_indices: pl.DataFrame, max_workers=1):
  fastqs = mapped_indices.get_column(FASTQ_COL).to_list()
  all_destination_uris = mapped_indices.get_column(DESTINATION_URI_COL).to_list()
  all_copy_args = zip(fastqs, all_destination_uris)
  with ProcessPoolExecutor(max_workers=max_workers) as copy_pool:
    futures = [
        copy_pool.submit(copy_one_fastq, one_arg_set)
        for one_arg_set in all_copy_args
    ]
    with Progress() as progress:
      copy_progress_bar = progress.add_task("Copying Fastq Files...", total=len(futures))
      for future in as_completed(futures):
        future.result()  # raises exceptions if any
        progress.update(copy_progress_bar, advance=1)


def get_tail_filename(fastq):
  """
  :return str: returns everything after the last S1, S2, S3 number in the filename.
  """
  external_name_matches = [match for match in re.finditer(EXTERNAL_NAME_REGEX, fastq)]
  tail_filename = fastq[external_name_matches[-1].start():]
  return tail_filename


def get_all_destination_uris(mapped_indices: pl.DataFrame, project_code: str, other_fastqs: list[str]) -> pl.DataFrame:
  DESTINATION_FOLDER = join(BASE_CLOUD_DIR, project_code, PROJECT_FASTQS_DIR)
  r1_fastq_col = mapped_indices.get_column(FASTQ_COL).to_list()
  collaborator_name_col = mapped_indices.get_column(COLLABORATOR_SAMPLE_NAMES_COL).to_list()
  all_r1_dest = []
  all_other_mapped_fastqs = []
  all_other_dest = []
  all_other_collaborator_names = []
  for fastq_index, r1_fastq in enumerate(r1_fastq_col):
    tail_filename = get_tail_filename(r1_fastq)
    r1_dest_name = join(DESTINATION_FOLDER, collaborator_name_col[fastq_index] + "_" + tail_filename)
    all_r1_dest.append(r1_dest_name)
    not_r1_read = [read_type.value for read_type in AllReadTypes]
    not_r1_read.remove(AllReadTypes.R1.value)
    not_r1_read = "|".join(not_r1_read)
    not_r1_read = f"({not_r1_read})"
    fastq_search = re.compile(f"{r1_fastq[:READ_TYPE_START]}{not_r1_read}{r1_fastq[READ_TYPE_END:]}")
    other_fastqs_for_sample = list(filter(fastq_search.match, other_fastqs))
    for other in other_fastqs_for_sample:
      other_tail = get_tail_filename(other)
      other_dest_name = join(DESTINATION_FOLDER, collaborator_name_col[fastq_index] + "_" + other_tail)
      all_other_mapped_fastqs.append(other)
      all_other_dest.append(other_dest_name)
      all_other_collaborator_names.append(collaborator_name_col[fastq_index])
  mapped_indices = mapped_indices.with_columns(pl.Series(name=DESTINATION_URI_COL, values=all_r1_dest))
  other_mapped_fastqs = pl.from_dict({
      FASTQ_COL: all_other_mapped_fastqs,
      DESTINATION_URI_COL: all_other_dest,
      COLLABORATOR_SAMPLE_NAMES_COL: all_other_collaborator_names
  })
  all_destinations = pl.concat([mapped_indices, other_mapped_fastqs], how="diagonal_relaxed")
  return all_destinations
