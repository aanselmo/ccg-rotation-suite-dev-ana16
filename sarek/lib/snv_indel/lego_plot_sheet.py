#!/usr/local/miniforge3/bin/python

"""
Graph a lego plot for the Mutation Lego PLot sheet
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import argparse
import io
import sys
from os.path import join, exists
from os import getcwd, remove, mkdir
from tempfile import TemporaryDirectory
from glob import glob
from pysam.bcftools import merge as vcf_merge, norm as vcf_norm
from SigProfilerMatrixGenerator.scripts.SigProfilerMatrixGeneratorFunc import SigProfilerMatrixGeneratorFunc
from sigProfilerPlotting import plotSBS
from pdf2image import convert_from_path

UNFILTERED_VCF_FILE_ENDING = ".mutect2.vcf.gz"


def setup_arguments():
  parser = argparse.ArgumentParser(description="Process all vcf files to turn it into a Signature Plot. For testing.")
  parser.add_argument("-v", "--vcfs", dest="vcfs", default=[], action="extend",
                      nargs="+", help=(
                          f"A list of vcf.gz files with the actual variant calls in them. Preferably the"
                          + f" unfiltered ones ending in {UNFILTERED_VCF_FILE_ENDING}"))
  parser.add_argument("-o", "--output-dir", type=str, default=join(getcwd(), "test"),
                      dest="output_dir", help="The name of the output directory.")
  parser.add_argument("-e", "--bed-file", dest="bed_file", type=str, required=True,
                      help="The targeted bed file for the bait set. If it's Whole genome sequencing, pass 'WGS'.")
  return parser.parse_args()


def merge_vcfs_and_plot(vcf_files: list[str], project_code: str, output_dir: str, genome="GRCh37",
                        bed_file=None) -> str:
  for one_vcf in vcf_files:
    if not (one_vcf.endswith(".vcf") or one_vcf.endswith(".vcf.gz")):
      raise ValueError("The files given do not end in .vcf or .vcf.gz")
  merged_vcf = join(output_dir, f"{project_code}_all_mutect2.vcf")
  vcf_merge(
      "--output", merged_vcf,
      "--output-type", "v", "--force-single", "--force-samples",
      *vcf_files,
      catch_stdout=False,
  )
  # Normalize multiallelic sites because it crashes SigProfilerMatrixGeneratorFunc if there are too many
  vcf_norm(
    "--output", join(output_dir, f"{project_code}_all_mutect2_filtered.vcf"),
    "--output-type", "v",
    "--multiallelics", "-any",
    merged_vcf,
    catch_stdout=False
  )
  remove(merged_vcf)
  SigProfilerMatrixGeneratorFunc(project_code, genome, output_dir, plot=False, bed_file=bed_file)
  matrix_path = join(output_dir, "output", "SBS", f"{project_code}.SBS384.region")
  plotSBS(matrix_path, join(output_dir, "plot"), project_code, "384")
  output_plot_pdf = glob(join(output_dir, "plot", "*"))[0]
  return output_plot_pdf


def graph_lego_plot(vcf_files: list[str], project_code: str, genome="GRCh37", bed_file=None):
  png_plot_mem_loc = io.BytesIO()
  if not vcf_files:
    raise ValueError("No VCF file was passed.")
  with TemporaryDirectory() as tempdir:
    try:
      output_plot_pdf = merge_vcfs_and_plot(vcf_files, project_code, tempdir, genome, bed_file)
      images = convert_from_path(output_plot_pdf, dpi=50)
      images[0].save(png_plot_mem_loc, "PNG")
    except Exception as e:
      print("WARNING: Graphing the lego plot had an error. Continuing...", file=sys.stderr)
      print(e, file=sys.stderr)
      pass
  return png_plot_mem_loc


def main():
  """
  Main function here was written as a way to quickly test parsing of 1 file. Otherwise, because the snv_indel_report
  is parallelized/multi-process, the stack trace of the error can be obfuscated.
  """
  args = setup_arguments()
  if not exists(args.output_dir):
    mkdir(args.output_dir)
  merge_vcfs_and_plot(args.vcfs, "test", args.output_dir, bed_file=args.bed_file)
  

if __name__ == "__main__":
  main()