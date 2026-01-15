#!/usr/local/miniforge3/bin/python

"""
Create a MultiQC report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Neil A. Patel"]

from os import makedirs
from pathlib import Path
from os.path import join, basename
from glob import glob
from datetime import date
from tempfile import TemporaryFile
from multiqc import run, BaseMultiqcModule
import polars as pl
import narwhals.stable.v2 as nw

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Custom Target QC',
            anchor='custom_target_qc',
            href='',
            info='Custom target coverage metrics'
        )

        # Load data
        df = pd.read_csv('data/my_metrics.tsv', sep='\t')

        df.set_index('sample', inplace=True)

        # ---------- Table ----------
        self.add_section(
            name='Target QC table',
            description='Per-sample target coverage metrics',
            content=self.render_table(
                headers={
                    'pct_on_target': {'title': '% on target'},
                    'pct_ge_30x': {'title': '% ≥30×'},
                    'uncovered_pct': {'title': '% uncovered'}
                },
                data=df.to_dict(orient='index')
            )
        )

        # ---------- Plotly plot ----------
        self.add_section(
            name='Coverage vs capture efficiency',
            plot=self.plot.scatter(
                df,
                x='pct_on_target',
                y='pct_ge_30x',
                title='% on target vs % ≥30×',
                xlab='% on target',
                ylab='% ≥30×',
                hover_label='index'
            )
        )


def create_multiqc_report(summary_table: nw.DataFrame, project_report_dir: str, output_file_name: str):
  with TemporaryFile("w") as summary_tempfile:
    summary_table.to_native().write_csv(summary_tempfile, separator="\t")
    script_dir = Path(__file__).resolve().parent
    run(
      analysis_dir=".",
      outdir=project_report_dir,
      config=join(script_dir, "multiqc.yml"),
      module=["multiqc_myreport"],
    )