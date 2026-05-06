#!/usr/local/miniforge3/bin/python

"""
A collection of functions to interact with MiniLIMS.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import json
import requests


DOMAIN = "https://ccgd2.dfci.harvard.edu/minilims/"
# project.php is a POST request with a secretKey and project in the body of the request (as multipart form data).
PROJECT_API_URL = DOMAIN + "api/get/10x-samplesheet/from/project.php"
# getIndicesFromSamples.php is a GET request with sample query params in the url
# Example: https://ccgd2.dfci.harvard.edu/minilims/api/getIndicesFromSamples.php?samples=["S_035150"]
INDICES_API_URL = DOMAIN + "api/getIndicesFromSamples.php"
# This API is used by the lab to generate a Project --> Sample (External ID, Sample Class)-->Library (Yield & Index)
PROJECT_METADATA_URL = DOMAIN + "api/getProjectMetadata.php"


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