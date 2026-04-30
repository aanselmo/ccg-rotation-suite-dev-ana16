#!/usr/local/miniforge3/bin/python

"""
Use this module for interacting with GCP.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


from os.path import join, abspath, dirname
import gzip
import io
from google.oauth2 import service_account
from google.cloud import storage


GCP_PROJECT = "ccgd-174313"
git_repo_dir = abspath(dirname(dirname(dirname(__file__))))
JSON_KEY = join(git_repo_dir, "secrets", "rnaseq_ccgd-174313-98def6888a16.json")


def get_storage_client():
  """
  Get's a GCP storage client connection so you can do the transfer.
  """
  the_credentials = service_account.Credentials.from_service_account_file(JSON_KEY)
  storage_client = storage.Client(project=GCP_PROJECT, credentials=the_credentials)
  return storage_client


def split_uri(gcp_uri: str) -> tuple[str, str]:
  """
  Split a GCP uri into it's bucket and the rest of the path.
  gs://mbcgs/x_folder/y_file.gz -> (mbcgs, x_folder/y_file.gz)
  """
  split_path = gcp_uri.split("/")
  bucket_name = split_path[2]
  rest_path = "/".join(split_path[3:])
  return (bucket_name, rest_path)


def get_blob_from_gs_uri(gcp_path: str, storage_client: storage.Client) -> storage.Blob:
  (bucket_name, rest_path) = split_uri(gcp_path)
  bucket = storage_client.bucket(bucket_name)
  blob_of_path = storage.Blob(bucket=bucket, name=rest_path)
  return blob_of_path


def get_gs_uri_from_blob(blob: storage.Blob) -> str:
  return f"gs://{blob.bucket.name}/{blob.name}"


def get_files_in_folder(client: storage.Client, gcp_folder_uri: str, recursive=True):
  (bucket, folder_uri) = split_uri(gcp_folder_uri)
  all_files = []
  if recursive:
    all_files = client.list_blobs(bucket, prefix=folder_uri)
  else:
    all_files = client.list_blobs(bucket, prefix=folder_uri, delimiter="/")
  return all_files


def get_head_of_gcp_gzip_file(client: storage.Client, gcp_file_uri: str, size_bytes=200):
  """
  Get the head of a .gz file that is on a Google bucket.
  :param file_path str: should follow GCP conventions e.g. gs://mbcgs/x_folder/y_file.gz
  :param size_bytes int: size in Bytes to read
  """
  blob = get_blob_from_gs_uri(gcp_file_uri, client)
  partial_bytes = blob.download_as_bytes(start=0, end=size_bytes)
  head_of_file = ""
  with gzip.GzipFile(fileobj=io.BytesIO(partial_bytes)) as f:
      head_of_file = f.read().decode("utf-8")
  return head_of_file
