#!/usr/local/miniconda3/bin/pytho:
"""
This script has the purpose of creating a co-mutation plot. See the following URL for reference:
https://iovs.arvojournals.org/article.aspx?articleid=2650631. There will be genes on the y-axis and sample names on
the y-axis.

Input data is a CSV that is formatted as the following:
sample_name,gene,mutation_type
sample1,KRAS,missense
sample1,KRAS,synonymous
sample2,TP53,missense

Unfortunately, you cannot have multiple colors in each square of the gridplot. Therefore, you have to set priority.
Please see the -h flag for help on the --order argument.
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
import dill

COL_SAMPLE_NAME = 'sample_name'
COL_GENE = 'gene'
COL_MUTATION_TYPE = 'mutation_type'

MUT_TYPE_NO_MUTATION = 'no_mutation'
MARK_BOTH_DUPLICATES_AS_TRUE = False
COLOR_MAP = 'tab10'

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


def heatmap(data, row_labels, col_labels, ax=None, cmap=None,
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
  im = ax.imshow(data, cmap=cmap, **kwargs)

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
      mutation_order_dict = dict(zip(mutation_types_in_priority, order_priority))
      mutation_types_in_priority.sort(key=mutation_order_dict.get)
  logging.debug("Mutation Types: " + str(mutation_types_in_priority) + "; Order: " 
      + str(list(range(0, len(mutation_types_in_priority)))))
  return mutation_types_in_priority


def remove_lower_priority_duplicates(df_gene_to_mutation_types, mutation_types_in_priority):
  """
  Let's say the mutation_types are "non, mis, syn" and the priorities are "3, 1, 2". That means that the priority
  order is "syn, non, mis". That means that the square for the color in the heatmap will be "syn" even if there are
  mutations for non and mis in the the same gene. The lower priorities mutation_types are removed here in this
  function so they are not graphed.
  """
  bool_series_duplicate_rows = df_gene_to_mutation_types.duplicated(subset=[COL_SAMPLE_NAME, COL_GENE],
      keep=MARK_BOTH_DUPLICATES_AS_TRUE)
  # First, get the clean ones that don't have duplicates.
  df_mutation_types_to_plot = df_gene_to_mutation_types.loc[~bool_series_duplicate_rows.values, :]
  # Then, loop through and add the priority ones of the duplicates.
  num_mut_types = len(mutation_types_in_priority)
  mutation_type_indices_to_loop =  num_mut_types - 1
  for mutation_type_index in range(0, mutation_type_indices_to_loop):
    mutation_type = mutation_types_in_priority[mutation_type_index]
    df_mutation_types_to_plot = df_mutation_types_to_plot.append(df_gene_to_mutation_types.loc[
        bool_series_duplicate_rows & (df_gene_to_mutation_types.loc[:, COL_MUTATION_TYPE] == mutation_type)])
    if not df_mutation_types_to_plot.empty:
      bool_series_duplicate_rows = df_mutation_types_to_plot.duplicated(subset=[COL_SAMPLE_NAME, COL_GENE],
        keep='last')
      # Don't Keep the last, actually. I want the bool_series to look like F, F, F, F, T so I can then just flip it, 
      # and remove the ones I don't need (which are the lower priority duplicates).
      df_mutation_types_to_plot = df_mutation_types_to_plot.loc[~bool_series_duplicate_rows, :]
  return df_mutation_types_to_plot


def convert_to_2d_df_for_heatmap(df_mutation_types_to_plot, mutation_types_in_priority):
  """
  The csv data is converted to a 2D dataframe that looks similar to the heatmap. The outputed array looks like the 
  following:
        Sample1 Sample2 Sample3
  TP53  1       0       3
  BRAF  3       2       0
  The numbers in the core of the data are the indexes in the mutation_types_in_priority. Also, any gene/sample pair 
  that had no mutation is filled in with 0. This no_mutation is added to the mutation_types_in_priority also and 
  returned.
  """
  for mutation_type in mutation_types_in_priority:
    df_mutation_types_to_plot.loc[(df_mutation_types_to_plot.loc[:, COL_MUTATION_TYPE] == mutation_type), 
        COL_MUTATION_TYPE] = mutation_types_in_priority.index(mutation_type)
  df_2d_mutation_types = df_mutation_types_to_plot.pivot(index=COL_GENE, columns=COL_SAMPLE_NAME, 
      values=COL_MUTATION_TYPE)
  df_2d_mutation_types = df_2d_mutation_types + 1
  mutation_types_in_priority = [MUT_TYPE_NO_MUTATION] + mutation_types_in_priority
  df_2d_mutation_types = df_2d_mutation_types.fillna(value=mutation_types_in_priority.index(MUT_TYPE_NO_MUTATION))
  return (df_2d_mutation_types, mutation_types_in_priority)

def main():
  LOGGING_FORMAT = "[%(filename)s:%(lineno)s] %(asctime)s %(message)s"
  args = setup_arguments()
  logging.basicConfig(filename=args.log_file_name, level=logging.DEBUG, format=LOGGING_FORMAT)
  df_gene_to_mutation_types = pd.read_csv(args.file_input_data)

  logging.info("Starting generation of co-mutation plot: " + str(args))
  mutation_types = df_gene_to_mutation_types.loc[:, COL_MUTATION_TYPE].tolist()
  logging.info("Retrieving priority of mutation_types.")
  mutation_types_in_priority = get_prioritized_mutation_types(mutation_types, args.order_priority)
  logging.info("Removing lower priority duplicates.")
  df_mutation_types_to_plot = remove_lower_priority_duplicates(df_gene_to_mutation_types, mutation_types_in_priority)
  (df_2d_mutation_types, mutation_types_in_priority) = convert_to_2d_df_for_heatmap(
      df_mutation_types_to_plot, mutation_types_in_priority)
  num_mut_types = len(mutation_types_in_priority)

  logging.info("Graphing the heatmap.")
  fig, ax3 = plt.subplots()
  y_labels_genes = list(df_2d_mutation_types.index.tolist())
  x_labels_samples = list(df_2d_mutation_types.columns.tolist())
  array_2d_mutation_types = df_2d_mutation_types.values
  a_mut_types = np.array(mutation_types_in_priority)
  
  # The data itself is categorical. Here we use a :class:`matplotlib.colors.BoundaryNorm` to get the data into classes
  # and use this to colorize the plot, but also to obtain the class labels from an array of classes.
  boundary_between_colors = matplotlib.colors.BoundaryNorm(np.linspace(-0.5, num_mut_types-0.5, 
      num_mut_types+1), num_mut_types)
  fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: a_mut_types[boundary_between_colors(x)])
  plot_color_map = plt.get_cmap(COLOR_MAP, num_mut_types)
  plot_color_map.colors[mutation_types_in_priority.index(MUT_TYPE_NO_MUTATION)] = np.array([1, 1, 1, 1])  # Set white.

  im, _ = heatmap(array_2d_mutation_types, y_labels_genes, x_labels_samples, ax=ax3,
                  cmap=plot_color_map, norm=boundary_between_colors,
                  cbar_kw=dict(ticks=np.arange(0, num_mut_types+1), format=fmt),
                  cbarlabel="Mutational Types")

  # I'm pickling the figure so you can load it and modify it later. Use code like the following to load it:
  # import dill
  # import matplotlib.pyplot as plt
  # ax = dill.load(open('test_file_for_comutation_plot_chart.pickle', 'rb'))
  # plt.show()
  dill.dump(ax3, open((args.file_input_data[:-4] + "_chart.pickle"), 'wb'))
  plt.savefig(args.file_input_data[:-4] + "_chart.png")
  logging.info("End generation of co-mutation plot.")

if __name__ == '__main__':
  main()
