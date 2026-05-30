#!/usr/local/miniforge3/bin/python

"""
This script contains common arguments used by most scripts in the sarek summary reports.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from argparse import ArgumentParser
from os import makedirs, cpu_count
from datetime import date
from os.path import join
from polars import Enum as PolarsEnum

BASE_CLOUD_DIR = "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek"
AVAILABLE_GENOME_VERSIONS = ["GRCh37", "GRCh38"]
DEFAULT_HTML_OUTPUT = "output.html"
DEFAULT_EXCEL_OUTPUT = "output.xlsx"
DEFAULT_CSV_OUTPUT = "output.csv"
AVAILABLE_CHROMOSOMES = PolarsEnum([
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "X", "Y", "MT", "GL000207.1", "GL000226.1", "GL000229.1", "GL000231.1", "GL000210.1", "GL000239.1",
    "GL000235.1", "GL000201.1", "GL000247.1", "GL000245.1", "GL000197.1", "GL000203.1", "GL000246.1", "GL000249.1",
    "GL000196.1", "GL000248.1", "GL000244.1", "GL000238.1", "GL000202.1", "GL000234.1", "GL000232.1", "GL000206.1",
    "GL000240.1", "GL000236.1", "GL000241.1", "GL000243.1", "GL000242.1", "GL000230.1", "GL000237.1", "GL000233.1",
    "GL000204.1", "GL000198.1", "GL000208.1", "GL000191.1", "GL000227.1", "GL000228.1", "GL000214.1", "GL000221.1",
    "GL000209.1", "GL000218.1", "GL000220.1", "GL000213.1", "GL000211.1", "GL000199.1", "GL000217.1", "GL000216.1",
    "GL000215.1", "GL000205.1", "GL000219.1", "GL000224.1", "GL000223.1", "GL000195.1", "GL000212.1", "GL000222.1",
    "GL000200.1", "GL000193.1", "GL000194.1", "GL000225.1", "GL000192.1"
])



def append_project_arg(parser: ArgumentParser):
  project_arg = parser.add_argument("-p", "--project", dest="project_code", type=str, default="",
                    help="The project code for the report. Used to generate filename with the project_code in it.")
  return project_arg


def append_common_args(parser: ArgumentParser, default_output_file=DEFAULT_CSV_OUTPUT):
  parser.add_argument("-o", "--output-file", type=str, default=default_output_file,
                      dest="output_file_name", help="The name of the output file. Overrides project_code argument.")
  project_arg = append_project_arg(parser)
  return project_arg


def append_bucket_arg(parser: ArgumentParser):
  parser.add_argument("-b", "--bucket-folder", dest="bucket_folder", default="", type=str,
                      help=("The sarek base folder with all of the other files within it. You can either enter this"
                            + "or pass all the files individually."))


def append_reference_arg(parser: ArgumentParser):
  parser.add_argument("-r", "--reference", dest="reference_genome", type=str, default="GRCh37",
                      choices=AVAILABLE_GENOME_VERSIONS,
                      help="The Genome Version. Options are 'GRCh37' or 'GRCh38'.")


def append_processes_arg(parser: ArgumentParser):
  this_pc_cpus = cpu_count() - 1 if cpu_count() else 1  # subtract 1 for this script.
  print(f"Process (CPU) Count Available: {this_pc_cpus}")
  parser.add_argument("-n", "--number-processes", dest="number_processes",
                      default=this_pc_cpus, type=int,
                      help=(f"Number of processes available on the multiprocessing sections."
                          + " Default is number of processes minus 1 (for this process)."
                          + " I haven't added threads because"
                          + " the documentation for google.cloud.storage says it's unlikely to change anything."))



def get_output_file_name(args, file_ending: str, default_file_name: str) -> tuple[str, str]:
  project_report_dir = join("reports", "ccg_reports")
  output_file_name = args.output_file_name
  today = date.today().strftime('%Y%m%d')
  project_code = "test"
  if output_file_name == default_file_name and args.project_code:
    project_code = args.project_code
    project_report_dir = join("reports", f"ccg_{project_code.lower()}_reports")
    makedirs(project_report_dir, exist_ok=True)
    output_file_name = join(project_report_dir, f"{today}_{project_code.lower()}_{file_ending}")
  return (project_code, output_file_name)
