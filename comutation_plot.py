#!/usr/local/miniconda3/bin/python
"""
This script has the purpose of creating a co-mutation plot. See the following URL for reference:
https://iovs.arvojournals.org/article.aspx?articleid=2650631. There will be genes on the y-axis and sample names on
the y-axis.

Input data is a CSV that is formatted as the following:
sample_name, gene, mutation_Type
sample1, KRAS, missense
sample1, KRAS, synonymous
sample2, TP53, missense

Set priority.

Set colors.

"""
__docformat__ = "reStructuredText"
__author__ = "Neil A. Patel and Bruce Wollison"
__credits__ = ["Neil A. Patel", "Bruce Wollison"]

# Come with python by default.
import argparse
import logging

# Install via pip.
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def setup_arguments():
  """
  This is the function that sets up the flags and the arguements you can pass to the script.
  :author: Neil A. Patel
  """
  aParser = argparse.ArgumentParser("Assigns reaction conditions to samples based on sample UDFs.")
  aParser.add_argument('-f', '--file', action='store', dest='file_input_data', required=True, 
      help='Please supply a csv file that has "sample_name, gene, mutation_type" in the title.')
  aParser.add_argument('-o', '--order', action='store', dest='order_priority', required=False,
      help='This is an option if you want to supply how colors overwrite each other in one square of the heatmap.'
      + 'For example, let\'s say you have mutation_types A, B, C, and D. You want C\'s color to overwrite the' 
      + 'other colors, you would supply 3,1,2,4 as a comma separated list with no spaces.'
    )
  aParser.add_argument('-g', '--log', action='store', dest='log_file_name', required=False,
      help="This is the log file.")

  return aParser.parse_args()

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
  """
  Create a heatmap from a numpy array and two lists of labels.
  :param data: A 2D numpy array of shape (N,M)
  :param row_labels: A list or array of length N with the labels for the rows
  :param col_labels: A list or array of length M with the labels for the columns
  :param ax: (optional) A matplotlib.axes.Axes instance to which the heatmap is plotted. If not provided, use 
      current axes or create a new one.
  :param cbar_kw: (optional) A dictionary with arguments to :meth:`matplotlib.Figure.colorbar`.
  :param cbarlabel: (optional) The label for the colorbar.
  :param **kwargs: (optional) All other arguments are directly passed on to the imshow call.
  """
  if not ax:
      ax = plt.gca()

  # Plot the heatmap
  im = ax.imshow(data, **kwargs)

  # Create colorbar
  cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
  cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

  # We want to show all ticks...
  ax.set_xticks(np.arange(data.shape[1]))
  ax.set_yticks(np.arange(data.shape[0]))
  # ... and label them with the respective list entries.
  ax.set_xticklabels(col_labels)
  ax.set_yticklabels(row_labels)

  # Let the horizontal axes labeling appear on top.
  ax.tick_params(top=True, bottom=False,
                 labeltop=True, labelbottom=False)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
           rotation_mode="anchor")

  ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
  ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
  ax.grid(which="minor", color="black", linestyle='-')
  ax.tick_params(which="minor", bottom=False, left=False)

  return im, cbar


def main():
  args = setup_arguments()
  LOGGING_FORMAT = "[%(filename)s:%(lineno)s] %(asctime)s %(message)s"
  # This is a conditional variable assignment. There is a default log file of LogGenerateCoMutationPlot.log 
  # if it is not supplied. 
  log_file = args.log_file_name if args.log_file_name else "LogGenerateCoMutationPlot.log" 
  logging.basicConfig(filename=log_file, level=logging.DEBUG, format=LOGGING_FORMAT)
  df_gene_to_mutation_types = pd.read_csv(args.input_file_name)

  # If there is 2 types of mutations for 1 gene, overwrite one with the other. This is done by alphabetical
  # order by default (e.g. indel, nonsense, missense, synonymous, ...). You can also supply the order with the
  fig, ax3 = plt.subplots()

  data = np.random.randn(6, 6)
  y = ["Prod. {}".format(i) for i in range(10, 70, 10)]
  x = ["Cycle {}".format(i) for i in range(1, 7)]

  qrates = np.array(list("ABCDEFG"))

  # The data itself is categorical. Here we use a :class:`matplotlib.colors.BoundaryNorm` to get the data into classes
  # and use this to colorize the plot, but also to obtain the class labels from an array of classes.
  norm = matplotlib.colors.BoundaryNorm(np.linspace(-3.5, 3.5, 8), 7)
  fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: qrates[::-1][norm(x)])

  im, _ = heatmap(data, y, x, ax=ax3,
                  cmap=plt.get_cmap("tab10", 7), norm=norm,
                  cbar_kw=dict(ticks=np.arange(-3, 4), format=fmt),
                  cbarlabel="Quality Rating")

  plt.show()


if __name__ == '__main__':
  main()
