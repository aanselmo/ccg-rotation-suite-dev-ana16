#!/usr/local/miniforge3/bin/python

"""
Graph a lego plot for the Mutation Lego PLot sheet
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


import io
import sys
from os.path import join
from tempfile import TemporaryDirectory
from glob import glob
from pysam.bcftools import merge as vcf_merge
from SigProfilerMatrixGenerator.scripts.SigProfilerMatrixGeneratorFunc import SigProfilerMatrixGeneratorFunc
from sigProfilerPlotting import plotSBS
from pdf2image import convert_from_path

UNFILTERED_VCF_FILE_ENDING = ".mutect2.vcf.gz"


def graph_lego_plot(vcf_files: list[str], project_code: str, genome="GRCh37", bed_file=None):
  png_plot_mem_loc = io.BytesIO()
  if not vcf_files:
    raise ValueError("No VCF file was passed.")
  with TemporaryDirectory() as tempdir:
    for one_vcf in vcf_files:
      if not (one_vcf.endswith(".vcf") or one_vcf.endswith(".vcf.gz")):
        raise ValueError("The files given do not end in .vcf or .vcf.gz")
    try:
      vcf_merge(
          "--output", join(tempdir, f"{project_code}_all_mutect2.vcf"),
          "--output-type", "v", "--force-single",
          *vcf_files,
          catch_stdout=False,
      )
      SigProfilerMatrixGeneratorFunc(project_code, genome, tempdir, plot=False, bed_file=bed_file)
      matrix_path = join(tempdir, "output", "SBS", f"{project_code}.SBS384.region")
      plotSBS(matrix_path, join(tempdir, "plot"), project_code, "384")
      output_plot_pdf = glob(join(tempdir, "plot", "*"))[0]
      images = convert_from_path(output_plot_pdf, dpi=50)
      images[0].save(png_plot_mem_loc, "PNG")
    except Exception as e:
      print("WARNING: Graphing the lego plot had an error. Continuing...", file=sys.stderr)
      print(e, file=sys.stderr)
      pass
  return png_plot_mem_loc
