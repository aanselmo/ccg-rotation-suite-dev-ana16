#!/usr/local/miniforge3/bin/python

"""
Graph a lego plot for the Mutation Lego PLot sheet
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]


from os.path import basename, join
from pathlib import Path
from tempfile import TemporaryDirectory
import io
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen


async def graph_lego_plot(vcf_files: list[str], project_code: str, genome="GRCh37", bed_file=None):
  with TemporaryDirectory() as tempdir:
    for one_vcf in vcf_files:
      Path(join(tempdir, basename(one_vcf))).symlink_to(one_vcf)
    matrices = matGen.SigProfilerMatrixGeneratorFunc(project_code, genome, tempdir, plot=True, bed_file=bed_file)
    breakpoint()