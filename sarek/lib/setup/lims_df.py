#!/usr/bin/env python

"""
Converts LIMS JSON into a DataFrame
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import polars as pl
from .lims import get_one_project_metadata, get_indices_of_samples


FIRST_BARCODE_LIMS_COL = "first_barcode"
SECOND_BARCODE_LIMS_COL = "second_barcode"
SECOND_BARCODE_REVCOMP_COL = "second_barcode_revcomp"
LIMS_SAMPLE_NAMES_COL = "lims_sample_names"
COLLABORATOR_SAMPLE_NAMES_COL = "sample"
LIMS_SUBJECT_COL = "subject_id"
PATIENT_COL = "patient"
STATUS_COL = "status"
REVCOMP_TABLE = str.maketrans("ACGT", "TGCA")


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
  sample_classes = []
  first_barcodes = []
  second_barcodes = []
  FIRST_ELEMENT = 0
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
    if type(one_external_subject_id) is not str:
      collaborator_subject_ids.append(f"P{one_external_subject_id}")
    else:
      collaborator_subject_ids.append(one_external_subject_id)
    if one_sample["sampleClass"] == "Normal":
      sample_classes.append(0)
    else:
      sample_classes.append(1)
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
  second_revcomp = map(get_reverse_complement, second_barcodes)
  lims_indices = pl.from_dict({
      LIMS_SAMPLE_NAMES_COL: lims_sample_names,
      COLLABORATOR_SAMPLE_NAMES_COL: collaborator_sample_names,
      FIRST_BARCODE_LIMS_COL: first_barcodes,
      SECOND_BARCODE_LIMS_COL: second_barcodes,
      SECOND_BARCODE_REVCOMP_COL: second_revcomp,
      LIMS_SUBJECT_COL: lims_subject_ids,
      PATIENT_COL: collaborator_subject_ids,
      STATUS_COL: sample_classes,
  }, schema={
      LIMS_SAMPLE_NAMES_COL: pl.String,
      COLLABORATOR_SAMPLE_NAMES_COL: pl.String,
      FIRST_BARCODE_LIMS_COL: pl.String,
      SECOND_BARCODE_LIMS_COL: pl.String,
      SECOND_BARCODE_REVCOMP_COL: pl.String,
      LIMS_SUBJECT_COL: pl.String,
      PATIENT_COL: pl.String,
      STATUS_COL: pl.Int8,  
  })
  return lims_indices
