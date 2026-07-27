#!/usr/local/miniforge3/bin/python

"""
A collection of functions to interact with MiniLIMS.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel, Brandon Sie"]

import json
import requests


# DOMAIN = "https://ccgd-build.dfci.harvard.edu/minilims_nap24/"
DOMAIN = "https://ccgd2.dfci.harvard.edu/minilims/"

# getIndicesFromSamples.php is a GET request with sample query params in the url
# Example: https://ccgd2.dfci.harvard.edu/minilims/api/getIndicesFromSamples.php?samples=["S_035150"]
INDICES_API_URL = DOMAIN + "api/getIndicesFromSamples.php"
# This API is used by the lab to generate a Project --> Sample (External ID, Sample Class)-->Library (Yield & Index)
PROJECT_METADATA_URL = DOMAIN + "api/getProjectMetadata.php"
# getCapturesFromLibraries takes a library ID (e.g. L_179592) and returns capture name and ID
#   e.g 260421_SCB1108l_TWHC_Pool3_AB
# Example: https://ccgd-build.dfci.harvard.edu/minilims_nap24/api/getCapturesFromLibraries.php?libraries=["L_179592"]
LIBRARY_CAPTURES_URL = DOMAIN + "/api/getCapturesFromLibraries.php"
# Get Illumina Run ID from Project code. e.g. SCB1108m -> 20260618_LH00547_0000_N\/A260622_AB14251
# e.g. https://ccgd2.dfci.harvard.edu/minilims/api/getIlluminaRunsFromProjects.php?projects=[%22SCB1108m%22]
ILLUMINA_RUNS_URL = DOMAIN + "api/getIlluminaRunsFromProjects.php"
ILLUMINA_RUNS_SHORT_URL = DOMAIN + "api/getIlluminaRunsFromRunNames.php"

FIRST_ELEMENT = 0

def get_all_projects_metadata(all_projects: list[str]) -> list:
  """
  Project --> Sample (External ID, Sample Class) --> Library (Yield & Index). Note that the Index is just the name of
  of the Index - not the actual barcode sequence with ATGCs.
  :param all_projects list[str]: A list of project codes (parent + subproject e.g. SCB1108l).
  :return list: A list of Samples
  """
  projects_query_string = json.dumps(all_projects).replace(" ", "")
  response = requests.get(PROJECT_METADATA_URL, params={"projects": projects_query_string})
  all_projects_metadata = response.json()["data"]
  return all_projects_metadata


def get_one_project_metadata(one_project: str) -> list:
  """
  Get's metadata for 1 project only.
  """
  FIRST_PROJECT = 0
  project_metadata = get_all_projects_metadata([one_project])
  return project_metadata[FIRST_PROJECT]


def get_indices_of_samples(all_samples: list[str]) -> list:
  """
  :param all_samples list[str]: A list of samples in MiniLIMS (e.g. S_035162).
  :return list: A list of barcodes with their names and sequences (ATGCs).
  """
  samples_query_string = json.dumps(all_samples).replace(" ", "")
  response = requests.get(INDICES_API_URL, params={"samples": samples_query_string})
  all_indices_of_samples = response.json()["data"]
  return all_indices_of_samples


def get_capture_names(library: list[str]) -> list:
  """ 
  :param library list[str]: A list of library IDs in miniLIMS (e.g. L_179592)
  :return list: A list capture names (e.g HC_TW_250620_SCB1108j_CapPool1)
  """
  capture_query_string = json.dumps(library).replace(" ", "")
  response = requests.get(LIBRARY_CAPTURES_URL, params={"libraries": capture_query_string})
  all_capture_names = response.json()["data"]
  return all_capture_names


def get_run_names(project: list[str]) -> str:
  """
  :param project list[str]: A list of project codes (e.g. ["SCB1108m", "SCB1108k"])
  :return list[str]: Illumina run IDs from field "illuminaRuns" (e.g. 20260618_LH00547_0000_NA260622_AB14251)
  """
  capture_query_string = json.dumps(project).replace(" ", "")
  response = requests.get(ILLUMINA_RUNS_URL, params={"projects": capture_query_string})
  all_run_names = response.json()["data"]
  return all_run_names


def get_illumina_names(run_name: list[str]) -> str:
  capture_query_string = json.dumps(run_name).replace(" ", "")
  response = requests.get(ILLUMINA_RUNS_SHORT_URL, params={"names": capture_query_string})
  all_illumina_names = response.json()["data"][FIRST_ELEMENT]['illuminaRun'][FIRST_ELEMENT]
  return all_illumina_names


def get_lanes_from_illumina_names(illumina_name: list[str]): 
  capture_query_string = json.dumps(illumina_name).replace(" ", "")
  response = requests.get(LANES_URL, params={"runs": capture_query_string})
  response_text = response.text
  if response_text.startswith("Array"):
    response_text = response_text[len("Array"):]
    all_lanes = json.loads(response_text)
  else: 
    all_lanes = response_text.json()
  return all_lanes