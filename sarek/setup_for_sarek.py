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
from operator import ne as not_equal_operator
from tempfile import NamedTemporaryFile
from datetime import date
from os.path import join
import polars as pl
from lib.common import (
    append_project_arg, append_bucket_arg, BASE_CLOUD_DIR, CONSOLE_RED, CONSOLE_BOLD, CONSOLE_FONT_RESET
)
from lib.setup.google import get_storage_client, get_blob_from_gs_uri
from lib.setup.lims_df import (
    FIRST_BARCODE_LIMS_COL, SECOND_BARCODE_LIMS_COL, SECOND_BARCODE_REVCOMP_COL, get_all_lims_data,
    PATIENT_COL, STATUS_COL, COLLABORATOR_SAMPLE_NAMES_COL, get_reverse_complement
)
from lib.setup.fastqs import (
    FIRST_BARCODE_FASTQ_COL, SECOND_BARCODE_FASTQ_COL, copy_all_fastqs, get_all_destination_uris, get_all_fastq_indices,
    PROJECT_FASTQS_DIR, DESTINATION_URI_COL, AllReadTypes, START_LANE_NUMBER
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
                      help=("The MBCF Bucket Directory with all the fastqs in them"
                            + " (e.g. gs://mbcgs/260424_AB14103_fastq)."))
  parser.add_argument("-l", "--lims-override", dest="lims_override", type=argparse.FileType('r', encoding='UTF-8'),
                      help=("If you're working on an R&D project, it won't be in LIMS. Pass a .csv with columns:"
                            + " patient,status,sample,lane,first_barcode,second_barcode"))
  parser.add_argument("-d", "--max-hamming-distance", dest="max_hamming_distance", type=int, default=1,
                      help=("The maximum allowed hamming distance"
                            + " between the indices in the mbcf_dir and indices in lims."))
  return parser.parse_args()


def get_all_output_file_names(project_code: str) -> tuple[str, str]:
  today = date.today().strftime('%Y%m%d')
  fastq_map_name = join(PROJECT_FASTQS_DIR, f"{today}_{project_code.lower()}_mbcf_fastq_map.csv")
  samplesheet_name = join(PROJECT_CSV_DIR, f"{today}_{project_code.lower()}_sarek_samplesheet.csv")
  return (fastq_map_name, samplesheet_name)


def get_hamming(row: tuple):
  row_hamming = sum(map(not_equal_operator, (row[0] + row[1]), (row[2] + row[3])))
  rev_comp_hamming = sum(map(not_equal_operator, (row[0] + row[1]), (row[2] + row[4])))
  row_hamming = rev_comp_hamming if rev_comp_hamming < row_hamming else row_hamming
  return row_hamming


def get_hamming_column(all_potential_matches: pl.DataFrame) -> pl.DataFrame:
  all_potential_matches = all_potential_matches.select([
      FIRST_BARCODE_FASTQ_COL, SECOND_BARCODE_FASTQ_COL,
      FIRST_BARCODE_LIMS_COL, SECOND_BARCODE_LIMS_COL, SECOND_BARCODE_REVCOMP_COL
  ])
  hamming_results = all_potential_matches.map_rows(get_hamming)
  return hamming_results


def map_lims_to_fastq(
    project_code: str, bucket_folder: str, r1_fastq_indices: pl.DataFrame, lims_indices: pl.DataFrame,
    other_fastqs: list[str], max_hamming_distance: int
) -> pl.DataFrame:
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
  remaining_fastq_indices = r1_fastq_indices.join(mapped_indices, how="anti",
                                            on=[FIRST_BARCODE_FASTQ_COL, SECOND_BARCODE_FASTQ_COL])
  remaining_fastq_indices = remaining_fastq_indices.join(lims_indices, how="cross")
  hamming_results = get_hamming_column(remaining_fastq_indices)
  remaining_fastq_indices = pl.concat([remaining_fastq_indices, hamming_results], how="horizontal")
  remaining_fastq_indices = remaining_fastq_indices.filter(pl.col("map") <= max_hamming_distance)
  remaining_fastq_indices = remaining_fastq_indices.drop("map")
  mapped_indices = pl.concat([mapped_indices, remaining_fastq_indices], how="vertical")
  if bucket_folder:
    mapped_indices = get_all_destination_uris(mapped_indices, bucket_folder, other_fastqs)
  else:
    mapped_indices = get_all_destination_uris(mapped_indices, project_code, other_fastqs)
  remainder_fastqs = mapped_indices.shape[0] % lims_indices.shape[0]
  if remainder_fastqs != 0:
    num_missing_fastqs = lims_indices.shape[0] - remainder_fastqs
    warnings.warn(
        f"Some samples' fastqs were not found."
        f" {CONSOLE_RED}{CONSOLE_BOLD}There are {num_missing_fastqs} missing fastqs{CONSOLE_FONT_RESET}"
        + " if you assume there are an equal number of fastqs per sample.",
    )
  return mapped_indices


def create_sarek_samplesheet(mapped_indices: pl.DataFrame) -> pl.DataFrame:
  samplesheet = mapped_indices.select([PATIENT_COL, STATUS_COL, COLLABORATOR_SAMPLE_NAMES_COL, DESTINATION_URI_COL])
  # Drop all others besides R1 & R2 fastqs if exist
  r1_r2_regex = f"({AllReadTypes.R1.value}|{AllReadTypes.R2.value})"
  r1_r2_regex = f"{r1_r2_regex}_001.fastq.gz$"
  r1_regex = f"{AllReadTypes.R1.value}_001.fastq.gz$"
  r2_regex = f"{AllReadTypes.R2.value}_001.fastq.gz$"
  samplesheet = samplesheet.filter(
      pl.col(DESTINATION_URI_COL).str.contains(r1_r2_regex)
  )
  # Group fastqs by unique sample name
  samplesheet = samplesheet.group_by(COLLABORATOR_SAMPLE_NAMES_COL).agg(
      pl.col(DESTINATION_URI_COL),
      pl.col(PATIENT_COL).drop_nulls(),
      pl.col(STATUS_COL).drop_nulls(),
  )
  samplesheet = samplesheet.with_columns(
      pl.col(PATIENT_COL).list.first(),
      pl.col(STATUS_COL).list.first()
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
  samplesheet = samplesheet.with_columns(pl.col(FASTQ1_COL).str.slice(START_LANE_NUMBER, 4).alias(LANE_COL))
  samplesheet = samplesheet.select([
      PATIENT_COL, STATUS_COL, COLLABORATOR_SAMPLE_NAMES_COL, LANE_COL, FASTQ1_COL, FASTQ2_COL
  ])
  return samplesheet


def write_output_files(args, gcp_storage_client,
                       mapped_indices: pl.DataFrame, fastq_map_name,
                       sarek_samplesheet: pl.DataFrame, samplesheet_name):
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
  with NamedTemporaryFile(mode='w+') as samplesheet_temp_file:
    sarek_samplesheet.write_csv(samplesheet_temp_file.name, quote_style="non_numeric")
    if args.bucket_folder:
      upload_location = get_blob_from_gs_uri(gcp_storage_client,
                                             join(BASE_CLOUD_DIR, args.bucket_folder, samplesheet_name))
    else:
      upload_location = get_blob_from_gs_uri(gcp_storage_client,
                                             join(BASE_CLOUD_DIR, args.project_code, samplesheet_name))
    upload_location.upload_from_filename(samplesheet_temp_file.name)


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (fastq_map_name, samplesheet_name) = get_all_output_file_names(args.project_code)
  gcp_storage_client = get_storage_client()
  mbcf_dir = args.mbcf_dir
  print("Get indices and collaborator names from LIMS.")
  lims_indices = pl.DataFrame()
  if args.lims_override:
    lims_indices = pl.read_csv(args.lims_override)
    rev_comp = map(get_reverse_complement, lims_indices.get_column(SECOND_BARCODE_LIMS_COL).to_list())
    rev_comp = pl.Series(name=SECOND_BARCODE_REVCOMP_COL, values=rev_comp)
    lims_indices.insert_column(-1, rev_comp)
  else:
    lims_indices = get_all_lims_data(args.project_code)
  print("Get indices in fastqs.")
  (r1_fastq_indices, other_fastqs) = get_all_fastq_indices(gcp_storage_client, mbcf_dir, args.number_processes)
  r1_fastq_indices = r1_fastq_indices.filter(pl.col(FIRST_BARCODE_FASTQ_COL) != "")
  r1_fastq_indices = r1_fastq_indices.filter(pl.col(SECOND_BARCODE_FASTQ_COL) != "")
  print("Mapping collaborator names in LIMS to fastq indices.")
  mapped_indices = map_lims_to_fastq(
      args.project_code, args.bucket_folder, r1_fastq_indices, lims_indices, other_fastqs, args.max_hamming_distance
  )
  sarek_samplesheet = create_sarek_samplesheet(mapped_indices)
  write_output_files(args, gcp_storage_client, mapped_indices, fastq_map_name, sarek_samplesheet, samplesheet_name)
  print("Copy the fastqs to the project directory.")
  copy_all_fastqs(mapped_indices, args.number_processes)
  print("Done.")


if __name__ == "__main__":
  main()
