#!/opt/conda/bin/python
"""
Write the SNV Indel data to CSV files.
"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel"

import io
import os
import narwhals as nw


def write_csv(
    file_name: str,
    all_variants: nw.DataFrame,
    variant_summary_table: nw.DataFrame,
    mouse_contamination_table: nw.DataFrame,
    all_hamas: nw.DataFrame,
    legend: nw.DataFrame,
    lego_plot: io.BytesIO,
):
  """
  :param lego_plot io.BytesIO: lego_plot.getbuffer().nbytes == 0 when there is no plot.
  """
  pl_variants = all_variants.to_native()
  pl_summary = variant_summary_table.to_native()
  pl_legend = legend.to_native()
  pl_mouse_contamination = mouse_contamination_table.to_native()
  pl_all_hamas = all_hamas.to_native()

  base, _ = os.path.splitext(file_name)

  pl_variants.write_csv(f"{base}_snvs_and_indels.csv")
  pl_summary.write_csv(f"{base}_variant_summary_table.csv")
  pl_legend.write_csv(f"{base}_legend.csv")
  if not pl_all_hamas.is_empty():
    pl_all_hamas.write_csv(f"{base}_mouse_alleles_hamas.csv")
  if not pl_mouse_contamination.is_empty():
    pl_mouse_contamination.write_csv(f"{base}_mouse_contamination.csv")

  if lego_plot.getbuffer().nbytes > 0:
    with open(f"{base}_lego_plot.png", "wb") as f:
      f.write(lego_plot.getvalue())
