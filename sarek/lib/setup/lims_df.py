#!/usr/bin/env python

"""
Converts LIMS JSON into a DataFrame
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import polars as pl
import re
from .lims import get_one_project_metadata, get_indices_of_samples, get_capture_names, get_run_names, get_illumina_names
from ..common import debug_df

LIMS_SAMPLE_NAMES_COL = "lims_sample_names"
COLLABORATOR_SAMPLE_NAMES_COL = "sample"
LIBRARY_NAMES_COL = "library"
BARCODE_SET_COL = "barcode_set"
BARCODE_MEMBER_COL = "barcode_member"
FIRST_BARCODE_LIMS_COL = "first_barcode"
SECOND_BARCODE_LIMS_COL = "second_barcode"
SECOND_BARCODE_REVCOMP_COL = "second_barcode_revcomp"
CAPTURE_NAME_COL = "capture_name"
LIMS_SUBJECT_COL = "subject_id"
PATIENT_COL = "patient"
ORIGINAL_PATIENT_COL = "patient_original"
STATUS_COL = "status"
REVCOMP_TABLE = str.maketrans("ACGT", "TGCA")

FIRST_ELEMENT = 0

def get_reverse_complement(barcode: str) -> str:
  """
  A fast way to find the reverse complement of a sequence. See:
  https://stackoverflow.com/questions/8776924/finding-the-complement-of-a-dna-sequence
  """
  return barcode.translate(REVCOMP_TABLE)[::-1]


def get_all_lims_data(project_code: str) -> pl.DataFrame:
  # These are parallel arrays that will later be put into a DataFrame

  lims_sample_names = []
  collaborator_sample_names = []
  lims_subject_ids = []
  collaborator_subject_ids = []
  original_collaborator_subject_ids = []
  sample_classes = []

  project_metadata = get_one_project_metadata(project_code)
  all_samples = project_metadata["samples"]
  for one_sample in all_samples:
    lims_sample_names.append(one_sample["name"])
    one_external_id = one_sample["externalId"]
    if type(one_external_id) is not str:
      collaborator_sample_names.append(f"S{one_external_id}")
    else:
      collaborator_sample_names.append(one_external_id)
    lims_subject_ids.append(one_sample["subjectId"])
    one_external_subject_id = one_sample["externalSubjectId"]
    original_collaborator_subject_ids.append(f"{one_external_subject_id}")
    if type(one_external_subject_id) is not str:
      collaborator_subject_ids.append(f"P{one_external_subject_id}")
    else:
      collaborator_subject_ids.append(one_external_subject_id)
    if one_sample["sampleClass"] == "Normal":
      sample_classes.append(0)
    else:
      sample_classes.append(1)

  library_names = [""] * len(lims_sample_names)
  barcode_member_names = [""] * len(lims_sample_names)
  barcode_set_names = [""] * len(lims_sample_names)
  first_barcodes = [""] * len(lims_sample_names)
  second_barcodes = [""] * len(lims_sample_names)
  indices_response = get_indices_of_samples(lims_sample_names)
  for one_sample in indices_response:
    one_lims_sample_name = one_sample["name"]
    one_library_name = one_sample["libraries"][FIRST_ELEMENT]["name"]
    one_barcode_member_name = one_sample["libraries"][FIRST_ELEMENT]["index"]["name"]
    one_first_sequence = one_sample["libraries"][FIRST_ELEMENT]["index"]["firstSequence"][FIRST_ELEMENT]
    one_second_sequence = one_sample["libraries"][FIRST_ELEMENT]["index"]["secondSequence"][FIRST_ELEMENT]
    row_index = lims_sample_names.index(one_lims_sample_name)
    library_names[row_index] = one_library_name
    barcode_member_names[row_index] = one_barcode_member_name
    barcode_set_names[row_index] = re.sub("_[A-Z][0-9]+$", "", one_barcode_member_name)
    first_barcodes[row_index] = one_first_sequence
    second_barcodes[row_index] = one_second_sequence

  second_revcomp = map(get_reverse_complement, second_barcodes)
  capture_names = [""] * len(lims_sample_names)

  capture_response = get_capture_names(library_names)

  for one_library in capture_response:
    if ("libraryPools" in one_library):
      one_capture_name = one_library["libraryPools"][1]["hybridCaptures"][FIRST_ELEMENT]["name"]
    else:
      one_capture_name = ""  # set a placeholder default if no capture name is present
    row_index = library_names.index(one_library["name"])
    capture_names[row_index] = one_capture_name

  lims_indices = pl.from_dict({
      LIMS_SAMPLE_NAMES_COL: lims_sample_names,
      COLLABORATOR_SAMPLE_NAMES_COL: collaborator_sample_names,
      LIBRARY_NAMES_COL: library_names,
      BARCODE_SET_COL: barcode_set_names,
      BARCODE_MEMBER_COL: barcode_member_names,
      FIRST_BARCODE_LIMS_COL: first_barcodes,
      SECOND_BARCODE_LIMS_COL: second_barcodes,
      SECOND_BARCODE_REVCOMP_COL: second_revcomp,
      CAPTURE_NAME_COL: capture_names,
      LIMS_SUBJECT_COL: lims_subject_ids,
      PATIENT_COL: collaborator_subject_ids,
      ORIGINAL_PATIENT_COL: original_collaborator_subject_ids,
      STATUS_COL: sample_classes,
  }, schema={
      LIMS_SAMPLE_NAMES_COL: pl.String,
      COLLABORATOR_SAMPLE_NAMES_COL: pl.String,
      LIBRARY_NAMES_COL: pl.String,
      BARCODE_SET_COL: pl.String,
      BARCODE_MEMBER_COL: pl.String,
      FIRST_BARCODE_LIMS_COL: pl.String,
      SECOND_BARCODE_LIMS_COL: pl.String,
      SECOND_BARCODE_REVCOMP_COL: pl.String,
      CAPTURE_NAME_COL: pl.String,
      LIMS_SUBJECT_COL: pl.String,
      PATIENT_COL: pl.String,
      ORIGINAL_PATIENT_COL: pl.String,
      STATUS_COL: pl.Int8,
  })
  debug_df("lims_indices", lims_indices)
  return lims_indices

def get_illumina_run(project: str):
  project_run_name = get_run_names([project])
  if 'name' in project_run_name[FIRST_ELEMENT]:
    run_name = project_run_name[FIRST_ELEMENT]['illuminaRuns'][FIRST_ELEMENT]
    return get_illumina_names([run_name])
  else:
    print(project_run_name)
    print(f"Project ${project} not found in miniLIMS")
    return ""
  

