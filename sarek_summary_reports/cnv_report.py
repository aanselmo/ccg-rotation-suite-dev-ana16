#!/usr/local/miniforge3/bin/python

"""
Use this script to generate the CNV report for Sarek
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

import argparse
import gzip
from os import makedirs
from datetime import date
from os.path import dirname, abspath, basename, join
from glob import glob
import asyncio
import xlsxwriter
import polars as pl
import narwhals.stable.v2 as nw
from cnv_lib.qc_sheet import get_qc_metrics, graph_qc_metrics
from cnv_lib.variant_sheet import get_cnv_variants, CHROMOSOME_COL, START_COL, END_COL
from cnv_lib.gene_sheet import summarize_gene_level
from cnv_lib.chromosome_sheet import summarize_chromosome_arm_level, ARM_COL


AVAILABLE_GENOME_VERSIONS = ["GRCh37", "GRCh38"]


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all sarek files to create a summary report.")
  parser.add_argument("-o", "--output-file", type=str, default="output.xlsx",
                      dest="output_file_name", help="The name of the output file. Overrides project_code argument.")
  parser.add_argument("-p", "--project", dest="project_code", type=str, default="",
                      help="The project code for the report. Used to generate filename with the project_code in it.")
  parser.add_argument("-r", "--reference", "genome_reference", type=str, default="GRCh37",
                      choices=AVAILABLE_GENOME_VERSIONS,
                      help="The Genome Version. Options are 'GRCh37' or 'GRCh38'.")
  parser.add_argument("-b" "--bucket-folder", dest="bucket_folder", default="", type=str,
                      help=("The sarek base folder with all of the other files within it. You can either enter this"
                            + "or pass all the files individually."))
  parser.add_argument("-p", "--processes", type=int, required=False, default=3, dest="total_processes",
                      help=(
                          "The number of processes to allocate for parallel jobs. I haven't added threads because"
                          + " the documentation for google.cloud.storage says it's unlikely to change anything."
                          + " The default is 3 because there are 4 cpu cores on the nextflow4 server."
                      ))
  # These 3 types of files are for the QC metrics sheet
  parser.add_argument("-m" "--mos-files", "mos_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .recal.mosdepth.summary.txt files with mean coverage in them.")
  parser.add_argument("-r" "--cnr-files", "cnr_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .cnr files with normalized log2 ratios in them.")
  parser.add_argument("-s" "--cns-files", "cns_files", default=[], action="extend", nargs="+",
                      help=("A list of paths to .cns files with segment calls in them. Use the most unfiltered ones"
                            + " i.e. the .cns files not .call.cns or .tumor_only.call.cns"))
  parser.add_argument("-t" "--tumor-cns-files", "tumor_cns_files", default=[], action="extend", nargs="+",
                      help=("A list of paths to .cns files with segment calls in them. Use the .tumor_only.call.cns."
                            + "These files have a cn column (copy number calls)."))
  # Gene summary sheet
  parser.add_argument("-g" "--gene-files", "gene_files", default=[], action="extend", nargs="+",
                      help="A list of paths to .genemetrics.tsv files with gene summaries in them.")
  # Variant sheet
  parser.add_argument("-v" "--vcf-files", "vcf_files", default=[], action="extend", nargs="+",
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


def derive_chromosomal_arms(genome_version: str) -> nw.LazyFrame:
  """
  Derives the chromosomal arms from USCS cytoband file. Downloaded from:
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
  UCSC uses a half open system for coordinates for 0 based counting:
  https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
  :return nw.LazyFrame: A dataframe with the columns "chromosome", "arm", "start", "end"
  """
  l_genome = genome_version.lower()
  cytoband_data = nw.from_native(pl.LazyFrame({}))
  with gzip.open(f'{l_genome}CytoBand.txt.gz', 'rb') as cytoband_file:
    cytoband_data = pl.scan_csv(
        cytoband_file,
        separator="\t",
        has_header=False,
        new_columns=[CHROMOSOME_COL, START_COL, END_COL, "band", "stain"]
    )
    # Aggregate by arm substring
    cytoband_data = nw.from_native(cytoband_data)
    cytoband_data = cytoband_data.with_columns(arm=nw.col("band").str.slice(offset=0, length=1))  # first letter: p/q
    cytoband_data = cytoband_data.drop("band", "stain")
    cytoband_data = cytoband_data.group_by(CHROMOSOME_COL, ARM_COL).agg(nw.col(START_COL).min(), nw.col(END_COL).max())
  return cytoband_data.collect()


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
  chromosomal_arms = derive_chromosomal_arms(args.genome_version)
  (mos_files, cnr_files, cns_files, tumor_cns_files, gene_files, vcf_files) = gather_input_files(args)
  qc_task: asyncio.Task | None = None
  gene_summary_task: asyncio.Task | None = None
  arm_summary_task: asyncio.Task | None = None
  qc_metrics = nw.from_native(pl.LazyFrame({}))
  cnv_variants = nw.from_native(pl.LazyFrame({}))
  gene_summary = nw.from_native(pl.LazyFrame({}))
  arm_summary = nw.from_native(pl.LazyFrame({}))
  # Start Tasks with asyncio
  if mos_files and cnr_files and cns_files:
    qc_task = asyncio.create_task(
        get_qc_metrics(mos_files, cnr_files, cns_files, args.total_processes, args.bucket_dir)
    )
  if gene_files:
    gene_summary_task = asyncio.create_task(summarize_gene_level(gene_files))
  if vcf_files:
    cnv_variants = await get_cnv_variants(vcf_files, tumor_cns_files, args.total_processes)
  if cnv_variants:
    arm_summary_task = asyncio.create_task(summarize_chromosome_arm_level(chromosomal_arms, cnv_variants))
  # Await Tasks
  if qc_task is not None:
    qc_metrics = await qc_task
    (mtc_hist_mem_loc, mad_hist_mem_loc) = await graph_qc_metrics(qc_metrics)
  if gene_summary_task is not None:
    gene_summary = await gene_summary_task
  if arm_summary_task is not None:
    arm_summary = await arm_summary_task
  

if __name__ == "__main__":
  asyncio.run(main())
