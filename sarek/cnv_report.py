#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the CNV report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
from os.path import basename
from pathlib import Path
from glob import glob
import asyncio
import polars as pl
import narwhals.stable.v2 as nw
from lib.common import (
    append_common_args, append_bucket_arg, append_reference_arg, get_output_file_name, DEFAULT_EXCEL_OUTPUT
)
from lib.cnv.qc_sheet import get_qc_metrics, graph_qc_metrics
from lib.cnv.variant_sheet import get_cnv_variants
from lib.cnv.gene_sheet import summarize_gene_level
# from lib.cnv.chromosome_sheet import summarize_chromosome_arm_level  # Arm Level Summary tab disabled
from lib.cnv.armlevel_process import get_armlevel_calls
from lib.cnv.write_xlsx import write_xlsx


DEFAULT_CNV_REPORT_FILE_ENDING = "cnv_report.xlsx"
# Reference tabs copied into the report as-is (tab names preserved). Resolved relative to this
# script so the report can be generated from any working directory.
DOCS_DIR = Path(__file__).resolve().parent / "docs"
REFERENCE_SHEET_FILES = [
    DOCS_DIR / "cnv_arms_level_key.xlsx",
    DOCS_DIR / "cnv_arms_level_flags.xlsx",
]


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a CNV summary report.")
  append_common_args(parser, default_output_file=DEFAULT_EXCEL_OUTPUT)
  append_bucket_arg(parser)
  append_reference_arg(parser)
  # These 3 types of files are for the QC metrics sheet
  parser.add_argument("-m", "--mos-files", dest="mos_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .recal.mosdepth.summary.txt files with mean coverage in them.")
  parser.add_argument("-c", "--cnr-files", dest="cnr_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .cnr files with normalized log2 ratios in them.")
  parser.add_argument("-d", "--call-cns-files", dest="call_cns_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .call.cns files with deletion, amplification, or neutral calls in them.")
  parser.add_argument("-s", "--cns-files", dest="cns_files", default=[], action="extend", nargs="+",
                      help=("A list of paths to .cns files with segment calls in them. Use the most unfiltered ones"
                            + " i.e. the .cns files not .call.cns or .tumor_only.call.cns"))
  parser.add_argument("-t", "--tumor-cns-files", dest="tumor_cns_files", default=[], action="extend", nargs="+",
                      help=("A list of paths to .cns files with segment calls in them. Use the .tumor_only.call.cns."
                            + "These files have a cn column (copy number calls)."))
  # Gene summary sheet
  parser.add_argument("-g", "--gene-files", dest="gene_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .genemetrics.tsv files with gene summaries in them.")
  # Variant sheet
  parser.add_argument("-v", "--vcf-files", dest="vcf_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .cnvcall.vcf files with the actual variant calls in them.")
  validate_arguments(parser.parse_args(), parser)
  return parser.parse_args()


def validate_arguments(args, parser: argparse.ArgumentParser):
  """
  Validate that you entered in the arguments correctly. If they're not, it raises an error and the
  script exits.
  """
  if args.bucket_folder:
    if (args.mos_files or args.cnr_files or args.cns_files or args.gene_files or args.vcf_files):
      parser.error("No files should be specified if --bucket-folder is specified.")
  else:
    if not (args.mos_files or args.cnr_files or args.cns_files or args.gene_files or args.vcf_files):
      parser.error("Please specify a --bucket-folder or individual files (--mos-files, --cnr-files, --cns-files, "
                   + "--tumor-cns-files, --gene-files, --vcf-files).")


def gather_input_files(args: argparse.Namespace) -> tuple[list[str], list[str], list[str], list[str], list[str]]:
  """
  Gather all the files.
  :param args argparse.Namespace: All arguments for the script
  :returns tuple[5 lists]: 5 lists of file paths (mos_files, cnr_files, cns_files, gene_files, vcf_files) 
  """
  mos_files = []
  cnr_files = []
  cns_files = []
  call_cns_files = []
  tumor_cns_files = []
  gene_files = []
  vcf_files = []
  if not args.mos_files:
    mos_files = glob(f"{args.bucket_folder}/reports/mosdepth/*/*.recal.mosdepth.summary.txt")
  else:
    mos_files = args.mos_files

  if not mos_files:
    raise ValueError("Mos Depth files are required and were not found.")

  if not args.cnr_files:
    cnr_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cnr")
  else:
    cnr_files = args.cnr_files

  if not args.call_cns_files:
    call_cns_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.call.cns")
    # Remove tumor_only.call.cns.  We want sample.call.cns files only.
    call_cns_files = list(filter(
        lambda one_file: "tumor_only" not in basename(one_file),
        call_cns_files
    ))
  else:
    call_cns_files = args.call_cns_files

  if not args.cns_files:
    cns_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cns")
    # Remove call.cns files which are basically duplicates. In this case, we want the ones without call in filenames
    # because they have more segment calls which will lower the median absolute deviation (MAD) scores.
    cns_files = list(filter(
        lambda one_file: "call" not in basename(one_file) and "bintest" not in basename(one_file),
        cns_files
    ))
  else:
    cns_files = args.cns_files
  if not args.tumor_cns_files:
    tumor_cns_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.tumor_only.call.cns")
    somatic_cns_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.somatic.call.cns")
    germline_cns_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.germline.call.cns")
    tumor_cns_files = tumor_cns_files + somatic_cns_files + germline_cns_files
  else:
    tumor_cns_files = args.tumor_cns_files

  if not args.gene_files:
    gene_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.genemetrics.tsv")
  else:
    gene_files = args.gene_files

  if not args.vcf_files:
    vcf_files = glob(f"{args.bucket_folder}/variant_calling/cnvkit/*/*.cnvcall.vcf")
  else:
    vcf_files = args.vcf_files
  return (mos_files, cnr_files, cns_files, tumor_cns_files, gene_files, vcf_files, call_cns_files)


async def main():
  """
  Runs the whole script.
  """
  args = setup_arguments()
  (project_code, output_file_name) = get_output_file_name(args, DEFAULT_CNV_REPORT_FILE_ENDING, DEFAULT_EXCEL_OUTPUT)
  (mos_files, cnr_files, cns_files, tumor_cns_files, gene_files, vcf_files, call_cns_files) = gather_input_files(args)
  qc_task: asyncio.Task | None = None
  gene_summary_task: asyncio.Task | None = None
  # arm_summary_task: asyncio.Task | None = None  # Arm Level Summary tab disabled
  armlevel_table_task: asyncio.Task | None = None
  qc_metrics = nw.from_native(pl.DataFrame({}))
  cnv_variants = nw.from_native(pl.DataFrame({}))
  gene_summary = nw.from_native(pl.DataFrame({}))
  # arm_summary = nw.from_native(pl.DataFrame({}))  # Arm Level Summary tab disabled
  armlevel_table = nw.from_native(pl.DataFrame({}))
  mtc_hist_mem_loc = None
  mad_hist_mem_loc = None
  # Start Tasks with asyncio
  if mos_files and cnr_files and cns_files:
    qc_task = asyncio.create_task(
        get_qc_metrics(mos_files, cnr_files, cns_files, args.bucket_folder)
    )
  if cnr_files and cns_files and call_cns_files:
    armlevel_table_task = asyncio.create_task(
        get_armlevel_calls(cnr_files, cns_files, call_cns_files, args.bucket_folder)
    )
  
  if gene_files:
    gene_summary_task = asyncio.create_task(summarize_gene_level(gene_files))
  if vcf_files:
    cnv_variants = await get_cnv_variants(vcf_files, tumor_cns_files)
  # Arm Level Summary tab disabled
  # if cnv_variants:
  #   arm_summary_task = asyncio.create_task(summarize_chromosome_arm_level(args.reference_genome, cnv_variants))
  # Await Tasks
  if qc_task is not None:
    qc_metrics = await qc_task
    (mtc_hist_mem_loc, mad_hist_mem_loc) = await graph_qc_metrics(qc_metrics)
  if gene_summary_task is not None:
    gene_summary = await gene_summary_task
  # if arm_summary_task is not None:  # Arm Level Summary tab disabled
  #   arm_summary = await arm_summary_task
  if armlevel_table_task is not None:
    armlevel_table = await armlevel_table_task
  missing_reference_files = [str(one_file) for one_file in REFERENCE_SHEET_FILES if not one_file.exists()]
  if missing_reference_files:
    raise FileNotFoundError(f"Reference sheet files not found: {missing_reference_files}")
  write_xlsx(output_file_name, qc_metrics, cnv_variants, gene_summary,
             # arm_summary,  # Arm Level Summary tab disabled
             armlevel_table, mtc_hist_mem_loc, mad_hist_mem_loc, reference_files=REFERENCE_SHEET_FILES)


if __name__ == "__main__":
  asyncio.run(main())
