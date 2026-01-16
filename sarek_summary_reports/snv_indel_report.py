#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the SNV/Indel report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import io
from os import makedirs, cpu_count
from os.path import join
from glob import glob
from datetime import date
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from lib.common_args import append_common_args, append_reference_arg
from lib.snv_indel.lego_plot_sheet import graph_lego_plot, UNFILTERED_VCF_FILE_ENDING
from lib.snv_indel.snv_indel_sheet import get_variants, VEP_FILE_ENDING
from lib.snv_indel.variant_summary_sheet import get_variant_summary
from lib.snv_indel.analysis_legend_sheet import get_legend
from lib.snv_indel.write_snv_xlsx import write_xlsx


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  this_pc_cpus = cpu_count()-1 if cpu_count() else 1  # subtract 1 for this script.
  print(f"Process (CPU) Count Available: {this_pc_cpus}")
  append_common_args(parser)
  append_reference_arg(parser)
  parser.add_argument("-a", "--annotated-vcf-files", dest="annotated_vcf_files", default=[], action="extend", nargs="+",
                      help=f"A list of paths to {VEP_FILE_ENDING} files with the actual variant calls in them.")
  parser.add_argument("-n", "--number-processes", dest="number_processes",
                      default=this_pc_cpus, type=int,
                      help=f"A list of paths to {VEP_FILE_ENDING} files with the actual variant calls in them.")
  parser.add_argument("-u", "--unfiltered-vcf-files", dest="unfiltered_vcf_files", default=[], action="extend",
                      nargs="+",
                      help=f"A list of paths to {UNFILTERED_VCF_FILE_ENDING} files with the unfiltered calls in them.")
  parser.add_argument("-e", "--bed-file", dest="bed_file", type=str, required=True,
                      help="The targeted bed file for the bait set. If it's Whole genome sequencing, pass 'WGS'.")
  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if args.bucket_folder:
    if (args.annotated_vcf_files or args.unfiltered_vcf_files):
      parser.error("No files should be specified if --bucket-folder is specified.")
  else:
    if not (args.annotated_vcf_files or args.unfiltered_vcf_files):
      parser.error("Please specify a --bucket-folder or individual files (--vcf-files).")


def gather_input_files(args: argparse.Namespace) -> list[str]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns list[str]: A list of strings for vcf_files paths
  """
  annotated_vcf_files = []
  unfiltered_vcf_files = []

  if not args.annotated_vcf_files:
    annotated_vcf_files = glob(f"{args.bucket_folder}/annotation/mutect2/*/*{VEP_FILE_ENDING}")
  else:
    annotated_vcf_files = args.annotated_vcf_files
  
  if not args.unfiltered_vcf_files:
    unfiltered_vcf_files = glob(f"{args.bucket_folder}/variant_calling/mutect2/*/*{UNFILTERED_VCF_FILE_ENDING}")
  else:
    unfiltered_vcf_files = args.unfiltered_vcf_files
  
  return (annotated_vcf_files, unfiltered_vcf_files)


def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  output_file_name = args.output_file_name
  today = date.today().strftime('%Y%m%d')
  project_code = "test"
  if output_file_name == "output.xlsx" and args.project_code:
    project_code = args.project_code.lower()
    project_report_dir = f"{project_code}_reports"
    makedirs(project_report_dir, exist_ok=True)
    output_file_name = join(project_report_dir, f"{today}_{project_code}_snv_indel_report.xlsx")
  (annotated_vcf_files, unfiltered_vcf_files) = gather_input_files(args)
  bed_file = None
  lego_plot = io.BytesIO()
  number_processes = args.number_processes
  if args.number_processes > 2:
    number_processes = args.number_processes - 1  # subtract 1 for graph_task
  if args.bed_file != "WGS":
    bed_file = args.bed_file
  executor = ProcessPoolExecutor(max_workers=1, mp_context=get_context("spawn"))
  graph_task = executor.submit(graph_lego_plot, unfiltered_vcf_files, project_code, args.reference_genome, bed_file)
  all_variants = get_variants(annotated_vcf_files, number_processes)
  variant_summary_table = get_variant_summary(all_variants)
  legend = get_legend(args.reference_genome)
  lego_plot = graph_task.result()
  write_xlsx(output_file_name, all_variants, variant_summary_table, legend, lego_plot)


if __name__ == "__main__":
  main()
