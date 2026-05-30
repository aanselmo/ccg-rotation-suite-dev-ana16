#!/data/python-venv/for_ccg_rotation_suite/bin/python

"""
Use this script to rename the fastqs sent by MBCF for Sarek and then create a Sarek Samplesheet. The fastqs received
from MBCF a lot of times have random names. This makes sure that the collaborator sample name is in the file name
for the fastq.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import argparse
import re
import warnings
from tempfile import NamedTemporaryFile
from datetime import date
from os import makedirs
from os.path import join
import polars as pl
from lib.common import append_project_arg, append_bucket_arg, DEFAULT_CSV_OUTPUT, BASE_CLOUD_DIR
from lib.setup.google import get_storage_client, get_blob_from_gs_uri
from lib.setup.lims_df import (
    FIRST_BARCODE_LIMS_COL, SECOND_BARCODE_LIMS_COL, SECOND_BARCODE_REVCOMP_COL, get_all_lims_data,
    PATIENT_COL, STATUS_COL, COLLABORATOR_SAMPLE_NAMES_COL
)
from lib.setup.fastqs import (
    FIRST_BARCODE_FASTQ_COL, SECOND_BARCODE_FASTQ_COL, copy_all_fastqs, get_all_destination_uris, get_all_fastq_indices,
    PROJECT_FASTQS_DIR, DESTINATION_URI_COL, AllReadTypes, BEFORE_LANE_START
)


PROJECT_CSV_DIR = "csv"
LANE_COL = "lane"
FASTQ1_COL = "fastq_1"
FASTQ2_COL = "fastq_2"


def setup_arguments():
  parser = argparse.ArgumentParser(description=(
      "Use this script to rename the fastqs sent by MBCF for Sarek and then create a Sarek Samplesheet. The fastqs"
      + " received from MBCF a lot of times have random names. This makes sure that the collaborator sample name"
      + " is in the file name for the fastq."
  ))
  project_arg = append_project_arg(parser)
  append_bucket_arg(parser)
  parser.add_argument("-n", "--number-processes", dest="number_processes",
                      default=200, type=int,
                      help=(f"Number of processes available on the multiprocessing sections."
                            + " Default is 200. This is what I've found optimal from testing."
                            + " Be careful! If a CPU heavy task is introduced to this script, lower the default."))
  project_arg.required = True
  parser.add_argument("-m", "--mbcf-dir", dest="mbcf_dir", type=str, default="", required=True,
                      help="The MBCF Bucket Directory with all the fastqs in them (e.g. gs://mbcgs/260424_AB14103_fastq).")
  return parser.parse_args()


def get_all_output_file_names(project_code: str) -> tuple[str, str]:
  today = date.today().strftime('%Y%m%d')
  fastq_map_name = join(PROJECT_FASTQS_DIR, f"{today}_{project_code.lower()}_mbcf_fastq_map.csv")
  samplesheet_name = join(PROJECT_CSV_DIR, f"{today}_{project_code.lower()}_sarek_samplesheet.csv")
  return (fastq_map_name, samplesheet_name)


def map_lims_to_fastq(project_code, bucket_folder, r1_fastq_indices, lims_indices, other_fastqs) -> pl.DataFrame:
  """
  Join tables together based on the indices from LIMS and the indices from the fastqs.
  """
  mapped_indices = r1_fastq_indices.join_where(
      lims_indices,
      (
          pl.col(FIRST_BARCODE_LIMS_COL).str.contains(pattern=pl.col(FIRST_BARCODE_FASTQ_COL))
          & (pl.col(SECOND_BARCODE_LIMS_COL).str.contains(pattern=pl.col(SECOND_BARCODE_FASTQ_COL))
             | pl.col(SECOND_BARCODE_REVCOMP_COL).str.contains(pattern=pl.col(SECOND_BARCODE_FASTQ_COL)))
      )
  )
  if bucket_folder:
    mapped_indices = get_all_destination_uris(mapped_indices, bucket_folder, other_fastqs)
  else:
    mapped_indices = get_all_destination_uris(mapped_indices, project_code, other_fastqs)
  if mapped_indices.shape[0] < lims_indices.shape[0]:
    warnings.warn("Some samples' fastqs were not found.")
  print(mapped_indices.head())
  return mapped_indices


def create_sarek_samplesheet(mapped_indices: pl.DataFrame) -> pl.DataFrame:
  samplesheet = mapped_indices.select([PATIENT_COL, STATUS_COL, COLLABORATOR_SAMPLE_NAMES_COL, DESTINATION_URI_COL])
  # Drop all others besides R1 & R2 fastqs if exist
  r1_r2_read = f"({AllReadTypes.R1.value}|{AllReadTypes.R2.value})"
  r1_r2_read = f"{r1_r2_read}_001.fastq.gz$"
  r1_r2_regex = re.compile(r1_r2_read)
  r1_regex = re.compile(f"{AllReadTypes.R1.value}_001.fastq.gz$")
  r2_regex = re.compile(f"{AllReadTypes.R2.value}_001.fastq.gz$")
  samplesheet = samplesheet.with_columns(
      pl.col(DESTINATION_URI_COL).str.contains(r1_r2_regex).alias(DESTINATION_URI_COL)
  )
  # Group fastqs by unique sample name
  samplesheet = samplesheet.group_by(COLLABORATOR_SAMPLE_NAMES_COL).agg(
      pl.col(DESTINATION_URI_COL),
      pl.col(PATIENT_COL).unique(),
      pl.col(STATUS_COL).unique(),
  )
  # R1 fastqs are fastq_1 and R2 fastqs are fastq_2
  samplesheet = samplesheet.with_columns(
      pl.col(DESTINATION_URI_COL).list.eval(pl.element().filter(
          pl.element().str.contains(r1_regex)
      )).list.first().alias(FASTQ1_COL),
      pl.col(DESTINATION_URI_COL).list.eval(pl.element().filter(
          pl.element().str.contains(r2_regex)
      )).list.first().alias(FASTQ2_COL),
  )
  samplesheet = samplesheet.with_columns(pl.col(FASTQ1_COL).str.slice(BEFORE_LANE_START, 3).alias(LANE_COL))
  samplesheet = samplesheet.select([
      PATIENT_COL, STATUS_COL, COLLABORATOR_SAMPLE_NAMES_COL, LANE_COL, FASTQ1_COL, FASTQ2_COL
  ])
  return samplesheet


def write_output_files(args, mapped_indices: pl.DataFrame, gcp_storage_client, fastq_map_name, samplesheet_name):
  print("Writing mapping to output file.")
  with NamedTemporaryFile(mode='w+') as map_temp_file:
    mapped_indices.write_csv(map_temp_file.name)
    if args.bucket_folder:
      upload_location = get_blob_from_gs_uri(gcp_storage_client,
                                             join(BASE_CLOUD_DIR, args.bucket_folder, fastq_map_name))
    else:
      upload_location = get_blob_from_gs_uri(gcp_storage_client,
                                             join(BASE_CLOUD_DIR, args.project_code, fastq_map_name))
    upload_location.upload_from_filename(map_temp_file.name)
  print("Writing sarek samplesheet to output file.")
  with NamedTemporaryFile(mode='w+') as map_temp_file:
    mapped_indices.write_csv(map_temp_file.name)
    if args.bucket_folder:
      upload_location = get_blob_from_gs_uri(gcp_storage_client,
                                             join(BASE_CLOUD_DIR, args.bucket_folder, samplesheet_name))
    else:
      upload_location = get_blob_from_gs_uri(gcp_storage_client,
                                             join(BASE_CLOUD_DIR, args.project_code, samplesheet_name))
    upload_location.upload_from_filename(samplesheet_name.name)


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (fastq_map_name, samplesheet_name) = get_all_output_file_names(args.project_code)
  gcp_storage_client = get_storage_client()
  mbcf_dir = args.mbcf_dir
  print("Get indices and collaborator names from LIMS.")
  lims_indices = get_all_lims_data(project_code)
  print(lims_indices.head())
  print("Get indices in fastqs.")
  (r1_fastq_indices, other_fastqs) = get_all_fastq_indices(gcp_storage_client, mbcf_dir, args.number_processes)
  r1_fastq_indices = r1_fastq_indices.filter(pl.col(FIRST_BARCODE_FASTQ_COL) != "")
  r1_fastq_indices = r1_fastq_indices.filter(pl.col(SECOND_BARCODE_FASTQ_COL) != "")
  print(r1_fastq_indices.head())
  print("Mapping collaborator names in LIMS to fastq indices.")
  mapped_indices = map_lims_to_fastq(
      args.project_code, args.bucket_folder, r1_fastq_indices, lims_indices, other_fastqs
  )
  write_output_files(args, mapped_indices, gcp_storage_client, fastq_map_name, samplesheet_name)
  print("Copy the fastqs to the project directory.")
  copy_all_fastqs(mapped_indices, args.number_processes)
  print("Done.")

if __name__ == "__main__":
  main()
