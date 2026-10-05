#!/usr/local/miniforge3/bin/python
"""
This is a script that tests the write_xlsx.py script.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"


import pytest
from os.path import join
from pathlib import Path
from tempfile import NamedTemporaryFile
import math
import xlsxwriter
import polars as pl
import narwhals.stable.v2 as nw
import sarek.lib.snv_indel.vcf2df as vcf2df
import sarek.lib.snv_indel.snv_indel_sheet as snv_indel_sheet
import sarek.lib.snv_indel.variant_summary_sheet as variant_summary_sheet
import sarek.lib.snv_indel.lego_plot_sheet as lego_plot_sheet
import sarek.lib.snv_indel.analysis_legend_sheet as analysis_legend_sheet
import sarek.lib.snv_indel.write_snv_xlsx as write_snv_xlsx
from sarek.lib.excel_style import HEADER_STYLE


EXAMPLE_VCF_PATH = join(Path(__file__).parent, "data", "A008_X108_58797-1-C1.mutect2.filtered_VEP.ann.vcf.gz")
VCF_A016_PATH = join(Path(__file__).parent, "data", "A016_X116_BS07G52324-G1.mutect2.vcf.gz")
EXAMPLE_BED = join(Path(__file__).parent, "data", "POPv3_AllTargetedRegions_20251230.bed")


def test_write_snv_xlsx():
  all_variants = snv_indel_sheet.get_variants([EXAMPLE_VCF_PATH])
  variant_summary_table = variant_summary_sheet.get_variant_summary(all_variants)
  plot_memory_location = lego_plot_sheet.graph_lego_plot([VCF_A016_PATH], "SCB1118c", "GRCh37", bed_file=EXAMPLE_BED)
  legend_df = analysis_legend_sheet.get_legend("GRCh37")
  with NamedTemporaryFile() as tempfile:
    write_snv_xlsx.write_xlsx(
        tempfile.name,
        all_variants,
        variant_summary_table,
        nw.from_native(pl.DataFrame({})),
        nw.from_native(pl.DataFrame({})),
        legend_df,
        plot_memory_location
    )
    verify_summary_table = pl.read_excel(tempfile.name, sheet_name=write_snv_xlsx.VARIANT_SUMMARY_SHEET)
    verify_all_variants = pl.read_excel(tempfile.name, sheet_name=write_snv_xlsx.SNV_INDEL_SHEET)
    assert verify_summary_table.filter(
        pl.col(snv_indel_sheet.SAMPLE_COL) == "A008_X108_58797-1-C1"
    ).get_column(variant_summary_sheet.BOTH_COL)[0] == 1
    # Assert 1 variant has the correct alternate allele
    assert verify_all_variants.filter(
        pl.col(vcf2df.POSITION_COL) == 41266124
    ).get_column(vcf2df.ALT_COL)[0] == "G"


class TestOver1MebiRowsExcel:
  def test_many_samples(self):
    """
    Test when samples together have more than 1048576 rows.
    """
    all_variants = vcf2df.convert_vcf_to_data_frame(EXAMPLE_VCF_PATH)
    duplication_multiplier = math.ceil(write_snv_xlsx.MAX_ROWS_EXCEL / all_variants.shape[0] * 1.1)
    pl_all_variants = all_variants.to_native()
    mocked_sample_names = pl.select(pl.concat_str([
        pl.lit("S"),
        pl.int_range(0, duplication_multiplier).cast(pl.String),
    ]).alias(snv_indel_sheet.SAMPLE_COL))
    duplicated_variants = mocked_sample_names.join(pl_all_variants, how="cross")
    with NamedTemporaryFile() as tempfile:
      print(tempfile.name)
      with xlsxwriter.Workbook(tempfile.name, options={"constant_memory": True, "use_zip64": True}) as workbook:
        header_format = workbook.add_format(HEADER_STYLE)
        write_snv_xlsx.write_all_snv_indel_sheets(
            workbook, duplicated_variants, write_snv_xlsx.SNV_INDEL_SHEET, header_format
        )
      verify_variants = pl.read_excel(tempfile.name, sheet_name=f"{write_snv_xlsx.SNV_INDEL_SHEET} 0")
      # Assert 1 variant has the correct alternate allele
      assert verify_variants.filter(
          pl.col(vcf2df.POSITION_COL).eq(41266124)
          & pl.col(snv_indel_sheet.SAMPLE_COL).eq("S1")
      ).get_column("ALT")[0] == "G"

  def test_one_sample(self):
    """
    Test when one sample has more than 1048576 rows.
    """
    all_variants = vcf2df.convert_vcf_to_data_frame(EXAMPLE_VCF_PATH)
    duplication_multiplier = math.ceil(write_snv_xlsx.MAX_ROWS_EXCEL / all_variants.shape[0] * 1.1)
    pl_all_variants = all_variants.to_native()
    mocked_sample = pl.DataFrame(data={snv_indel_sheet.SAMPLE_COL: "A008"})
    duplicated_variants = mocked_sample.join(pl_all_variants, how="cross")
    duplicated_variants = duplicated_variants.select(pl.all().repeat_by(duplication_multiplier).list.explode())
    with NamedTemporaryFile() as tempfile:
      print(tempfile.name)
      with xlsxwriter.Workbook(tempfile.name, options={"constant_memory": True, "use_zip64": True}) as workbook:
        header_format = workbook.add_format(HEADER_STYLE)
        write_snv_xlsx.write_all_snv_indel_sheets(
            workbook, duplicated_variants, write_snv_xlsx.SNV_INDEL_SHEET, header_format
        )
      verify_variants = pl.read_excel(tempfile.name, sheet_name=f"{write_snv_xlsx.SNV_INDEL_SHEET}_a")
      # Assert 1 variant has the correct alternate allele
      assert verify_variants.filter(
          pl.col(vcf2df.POSITION_COL).eq(41266124)
          & pl.col(snv_indel_sheet.SAMPLE_COL).eq("A008")
      ).get_column("ALT")[0] == "G"
