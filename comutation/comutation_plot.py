#!/usr/local/miniconda3/bin/pytho:
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

COL_SAMPLE_NAME = 'sample_name'
COL_GENE = 'gene'
COL_MUTATION_TYPE = 'mutation_type'

def setup_arguments():
  """
  This is the function that sets up the flags and the arguements you can pass to the script.
  :author: Neil A. Patel
  """
  a_parser = argparse.ArgumentParser("Get's the arguments for the comutation_plot tool.")
  a_parser.add_argument('-f', '--file', action='store', dest='file_input_data', required=True, 
      help='Please supply a csv file that has "sample_name, gene, mutation_type" in the title.')
  a_parser.add_argument('-o', '--order', action='store', dest='order_priority', required=False,
      help='This is an option if you want to supply how colors overwrite each other in one square of the heatmap.'
      + 'For example, let\'s say you have mutation_types A, B, C, and D. You want C\'s color to overwrite the' 
      + 'other colors, you would supply 3,1,2,4 as a comma separated list with no spaces.'
    )
  a_parser.add_argument('-g', '--log', action='store', dest='log_file_name', required=False, 
      default="LogGenerateCoMutationPlot.log", help="This is the log file.")

  return a_parser.parse_args()

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

def get_prioritized_mutation_types(list_mutation_types, order_priority):
  """
  If there is 2 types of mutations for 1 gene, overwrite one with the other. This is done by alphabetical
  order by default (e.g. indel, nonsense, missense, synonymous, ...). You can also supply the order with the
  order_priority.
  """
  mutation_types_in_priority = list(set(list_mutation_types))  # We need the unique values.
  mutation_types_in_priority.sort()  # Put in alphabetical order.

  if order_priority:
    order_priority = order_priority.split(",")
    # Convert the string list to an int list so we can use them as indexes for the mutation_types_in_priority list.
    order_priority = list(map(int, order_priority))
    compare_priority_length = len(order_priority) - len(mutation_types_in_priority)
    if compare_priority_length > 0:
      ERROR_TO_PRINT = ("There are {0} more numbers supplied in the --order flag than there are mutation_types in"
          + " the csv.".format(compare_priority_length))
      raise IndexError(ERROR_TO_PRINT)
      logging.error(ERROR_TO_PRINT)
    elif compare_priority_length < 0:
      ERROR_TO_PRINT = ("There are {0} fewer numbers supplied in the --order flag than there are mutation_types in"
          + " the csv.".format(-1 * compare_priority_length))
      raise IndexError(ERROR_TO_PRINT)
      logging.error(ERROR_TO_PRINT)
    else:
      mutation_types_in_priority.sort(key=order_priority)


def main():
  LOGGING_FORMAT = "[%(filename)s:%(lineno)s] %(asctime)s %(message)s"
  MARK_BOTH_DUPLICATES_AS_TRUE = False

  args = setup_arguments()
  # This is a conditional variable assignment. There is a default log file of LogGenerateCoMutationPlot.log 
  # if it is not supplied. 
  logging.basicConfig(filename=log_file, level=logging.DEBUG, format=LOGGING_FORMAT)
  df_gene_to_mutation_types = pd.read_csv(args.input_file_name)
  logging.info("Starting generation of co-mutation plot: " + str(args))

  bool_series_duplicate_rows = df_gene_to_mutation_types.duplicated(subset=[COL_SAMPLE_NAME, COL_GENE],
      keep=MARK_BOTH_DUPLICATES_AS_TRUE)
  mutation_types_in_priority = df_gene_to_mutation_types.loc[:, COL_MUTATION_TYPE].tolist()
  mutation_types_in_priority = get_prioritized_mutation_types(mutation_types_in_priority, args.order_priority)

  # First, get the clean ones that don't have duplicates.
  df_mutation_types_to_plot = df_gene_to_mutation_types[~bool_series_duplicate_rows.values, :]
  # Then, loop through and add the priority ones of the duplicates.
  num_mut_types = len(mutation_types_in_priority)
  mutation_type_indices_to_loop =  num_mut_types - 1
  for mutation_type_index in mutation_type_indices_to_loop:
    mutation_type = mutation_types_in_priority[mutation_type_index]
    df_mutation_types_to_plot.append(df_gene_to_mutation_types.loc[bool_series_duplicate_rows 
        & df_gene_to_mutation_types.loc[:, COL_MUTATION_TYPE] == mutation_type])
    bool_series_duplicate_rows = df_mutation_types_to_plot.duplicated(subset=[COL_SAMPLE_NAME, COL_GENE],
      keep='first')
    df_mutation_types_to_plot = df_mutation_types_to_plot[bool_series_duplicate_rows, :]
  
  fig, ax3 = plt.subplots()

  y_labels_genes = list(set(df_mutation_types_to_plot.loc[:, COL_GENE].tolist()))
  x_labels_samples = list(set(df_mutation_types_to_plot.loc[:, COL_SAMPLE_NAME].tolist()))
  a_mut_types = np.array(mutation_types_in_priority)
  for mutation_type in mutation_types_in_priority:
    df_mutation_types_to_plot.loc[mutation_type, COL_MUTATION_TYPE] = mutation_types_in_priority.index(mutation_type)
  array_2d_mutation_types = df_mutation_types_to_plot.pivot(index=COL_GENE, columns=COL_SAMPLE_NAME, 
      values=COL_MUTATION_TYPE)



  # The data itself is categorical. Here we use a :class:`matplotlib.colors.BoundaryNorm` to get the data into classes
  # and use this to colorize the plot, but also to obtain the class labels from an array of classes.
  norm = matplotlib.colors.BoundaryNorm(np.linspace(0.5, num_mut_types+0.5, 
      num_mut_types+1, num_mut_types))
  fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: a_mut_types[::-1][norm(x)])

  im, _ = heatmap(array_2d_mutation_types, y_labels_genes, x_labels_samples, ax=ax3,
                  cmap=plt.get_cmap("tab10", num_mut_types), norm=norm,
                  cbar_kw=dict(ticks=np.arange(0, num_mut_types+1), format=fmt),
                  cbarlabel="Mutational Types")

  plt.show()
  logging.info("End generation of co-mutation plot.")

if __name__ == '__main__':
  main()
