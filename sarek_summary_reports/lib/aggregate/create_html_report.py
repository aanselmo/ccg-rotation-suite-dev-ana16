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
import polars as pl
import narwhals.stable.v2 as nw
import report_creator as rc
import plotly.express as px


def create_html_report(summary_table: nw.DataFrame, project_code: str, output_file_name: str):
  with rc.ReportCreator(
      title=f"{project_code} Aggregation Metrics Report",
      description="A summary of aggregation metrics and a fingerprinting result report.",
      footer="Created by CCG's Bioinformaticians (Neil A. Patel, & Anthony Anselmo)",
  ) as report:
    view = rc.Block(
        rc.Text(
            """It was the best of times, it was the worst of times, it was the age of wisdom, it was the age of foolishness, it was the epoch of belief, it was the epoch of incredulity, it was the season of light, it was the season of darkness, it was the spring of hope, it was the winter of despair.""",
            label="Charles Dickens, A Tale of Two Cities",
        ),
        rc.Group(
            rc.Metric(
                heading="Answer to Life, The Universe, and Everything",
                value="42",
            ),
            rc.Metric(
                heading="Author",
                value="Douglas Adams",
            ),
        ),
        rc.Bar(
            px.data.medals_long(),
            x="nation",
            y="count",
            dimension="medal",
            label="Bar Chart - Olympic Medals",
        ),
        rc.Scatter(
            px.data.iris(),
            x="sepal_width",
            y="sepal_length",
            dimension="species",
            marginal="histogram",
            label="Scatter Plot - Iris",
        ),
    )

    report.save(view, output_file_name)
