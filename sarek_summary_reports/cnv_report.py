#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the CNV report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import gzip
import sys
from os import makedirs
from os.path import dirname, abspath, basename, join
from glob import glob
from datetime import date
import asyncio
import xlsxwriter
import polars as pl
import narwhals.stable.v2 as nw
from lib.cnv.qc_sheet import get_qc_metrics, graph_qc_metrics
from lib.cnv.variant_sheet import get_cnv_variants, CHROMOSOME_COL, START_COL, END_COL
from lib.cnv.gene_sheet import summarize_gene_level
from lib.cnv.chromosome_sheet import summarize_chromosome_arm_level, ARM_COL


AVAILABLE_GENOME_VERSIONS = ["GRCh37", "GRCh38"]


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  parser.add_argument("-o", "--output-file", type=str, default="output.xlsx",
                      dest="output_file_name", help="The name of the output file. Overrides project_code argument.")
  parser.add_argument("-p", "--project", dest="project_code", type=str, default="",
                      help="The project code for the report. Used to generate filename with the project_code in it.")
  parser.add_argument("-r", "--reference", dest="genome_reference", type=str, default="GRCh37",
                      choices=AVAILABLE_GENOME_VERSIONS,
                      help="The Genome Version. Options are 'GRCh37' or 'GRCh38'.")
  parser.add_argument("-b" "--bucket-folder", dest="bucket_folder", default="", type=str,
                      help=("The sarek base folder with all of the other files within it. You can either enter this"
                            + "or pass all the files individually."))
  # These 3 types of files are for the QC metrics sheet
  parser.add_argument("-m" "--mos-files", dest="mos_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .recal.mosdepth.summary.txt files with mean coverage in them.")
  parser.add_argument("-r" "--cnr-files", dest="cnr_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .cnr files with normalized log2 ratios in them.")
  parser.add_argument("-s" "--cns-files", dest="cns_files", default=[], action="extend", nargs="+",
                      help=("A list of paths to .cns files with segment calls in them. Use the most unfiltered ones"
                            + " i.e. the .cns files not .call.cns or .tumor_only.call.cns"))
  parser.add_argument("-t" "--tumor-cns-files", dest="tumor_cns_files", default=[], action="extend", nargs="+",
                      help=("A list of paths to .cns files with segment calls in them. Use the .tumor_only.call.cns."
                            + "These files have a cn column (copy number calls)."))
  # Gene summary sheet
  parser.add_argument("-g" "--gene-files", dest="gene_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .genemetrics.tsv files with gene summaries in them.")
  # Variant sheet
  parser.add_argument("-v" "--vcf-files", dest="vcf_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .vcf files with the actual variant calls in them.")
  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if args.bucket_folder:
    if (args.mos_files or args.cnr_files or args.cns_files, args.gene_files, args.vcf_files):
      parser.error("No files should be specified if --bucket-folder is specified.")


def gather_input_files(args: argparse.Namespace) -> tuple[list[str], list[str], list[str], list[str], list[str]]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns tuple[5 lists]: 5 lists of file paths (mos_files, cnr_files, cns_files, gene_files, vcf_files) 
  """
  mos_files = []
  cnr_files = []
  cns_files = []
  tumor_cns_files = []
  gene_files = []
  vcf_files = []
  if not args.mos_files:
    mos_files = glob(f"{args.bucket_folder}/reports/mosdepth/*/*.recal.mosdepth.summary.txt")
  else:
    mos_files = args.mos_files

  if not args.cnr_files:
    cnr_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cnr")
  else:
    cnr_files = args.cnr_files

  if not args.cns_files:
    cns_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cns")
    # Remove call.cns files which are basically duplicates. In this case, we want the ones without call in filenames
    # because they have more segment calls which will lower the median absolute deviation (MAD) scores.
    for one_file in cns_files:
      if "call" in one_file:
        cns_files.remove(one_file)
        if "tumor_only" in one_file:
          tumor_cns_files.append(one_file)
  else:
    cns_files = args.cns_files

  if not args.gene_files:
    gene_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.genemetrics.tsv")
  else:
    gene_files = args.gene_files

  if not args.vcf_files:
    vcf_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cnvcall.vcf")
  else:
    vcf_files = args.vcf_files
  return (mos_files, cnr_files, cns_files, tumor_cns_files, gene_files, vcf_files)


async def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  output_file_name = args.output_file_name
  today = date.today().strftime('%Y%m%d')
  if output_file_name == "output.xlsx" and args.project_code:
    project_report_dir = f"{args.project_code}_Reports"
    makedirs(project_report_dir, exist_ok=True)
    output_file_name = join(project_report_dir, f"{today}_{args.project_code.lower()}_cnvkit_report.xlsx")
  (mos_files, cnr_files, cns_files, tumor_cns_files, gene_files, vcf_files) = gather_input_files(args)
  qc_task: asyncio.Task | None = None
  gene_summary_task: asyncio.Task | None = None
  arm_summary_task: asyncio.Task | None = None
  qc_metrics = nw.from_native(pl.DataFrame({}))
  cnv_variants = nw.from_native(pl.DataFrame({}))
  gene_summary = nw.from_native(pl.DataFrame({}))
  arm_summary = nw.from_native(pl.DataFrame({}))
  # Start Tasks with asyncio
  if mos_files and cnr_files and cns_files:
    qc_task = asyncio.create_task(
        get_qc_metrics(mos_files, cnr_files, cns_files, args.bucket_dir)
    )
  if gene_files:
    gene_summary_task = asyncio.create_task(summarize_gene_level(gene_files))
  if vcf_files:
    cnv_variants = await get_cnv_variants(vcf_files, tumor_cns_files)
  if cnv_variants:
    arm_summary_task = asyncio.create_task(summarize_chromosome_arm_level(args.genome_reference, cnv_variants.lazy()))
  # Await Tasks
  if qc_task is not None:
    qc_metrics = await qc_task
    (mtc_hist_mem_loc, mad_hist_mem_loc) = await graph_qc_metrics(qc_metrics)
  if gene_summary_task is not None:
    gene_summary = await gene_summary_task
  if arm_summary_task is not None:
    arm_summary = await arm_summary_task
  with xlsxwriter.Workbook(output_file_name) as workbook:
    qc_metrics.to_native().write_excel(workbook=workbook, worksheet="QC Metrics")
    cnv_variants.to_native().write_excel(workbook=workbook, worksheet="Copy Number Variants (CNV)")
    gene_summary.to_native().write_excel(workbook=workbook, worksheet="Gene Level Summary")
    arm_summary.to_native().write_excel(workbook=workbook, worksheet="Arm Level Summary")

if __name__ == "__main__":
  asyncio.run(main())
