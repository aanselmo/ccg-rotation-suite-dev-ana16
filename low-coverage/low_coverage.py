#!/usr/local/miniconda3/bin/python

import os
import pandas as pd
import numpy as np
from pyensembl import EnsemblRelease

def get_first_ref(cell_string):
  """
  Get's the first ref| of the name column so you get a just 1 gene name.
  """
  split_string = cell_string.split("ref|")
  gene = cell_string
  if len(split_string) >= 2:
    gene = split_string[1].replace(",", "")
  return gene


class GeneLocator:
  LATEST_GRCH37_RELEASE = 75
  data = EnsemblRelease(LATEST_GRCH37_RELEASE)

  def get_gene(self, chrom, pos):
    """
    Get's the gene locus based on chromosome and position.
    """
    gene_names = self.data.gene_names_at_locus(contig=chrom, position=pos)
    if gene_names:
      return gene_names[0]
    else:
      return ""


def main():
  PROJECT_FOLDER = "/home/nap24/Documents/Pipeline/ccgd/AN0008a/GC_AN0008a_SE51038a_20210817/SE51038a"
  CHROM_COL = "chrom"
  START_COL = "start"
  END_COL = "end"
  LENGTH_COL = "length"
  NAME_COL = "name"
  GC_CONTENT_COL = "%gc"
  join_columns = [CHROM_COL, START_COL, END_COL, LENGTH_COL, NAME_COL, GC_CONTENT_COL]
  MEAN_COV_COL = "mean_coverage"
  data_type_dict = {CHROM_COL : "str", START_COL : "int64", END_COL : "int64", LENGTH_COL : "int32",
      NAME_COL : "str", GC_CONTENT_COL : "float", MEAN_COV_COL : "float"}
  MEAN_OF_MEAN_COV_COL = "mean_of_mean_coverage"

  combined_data = pd.DataFrame()
  for index, filename in enumerate(os.listdir(PROJECT_FOLDER)):
    interval_data = pd.read_csv(os.path.join(PROJECT_FOLDER, filename), sep="\t", dtype=data_type_dict)
    new_mean_cov_col = MEAN_COV_COL + "_" + filename.replace("_full_targets_hsIntervals.txt", "")
    interval_data.rename(columns={MEAN_COV_COL : new_mean_cov_col}, inplace=True)
    if index == 0:
      combined_data = interval_data.loc[:, join_columns] 
    combined_data = combined_data.merge(interval_data.loc[:, (join_columns + [new_mean_cov_col])], on=join_columns)

  data_columns = combined_data.columns.tolist()
  coverage_columns = []
  for column in data_columns:
    if MEAN_COV_COL in column:
      coverage_columns.append(column)

  mean_of_mean_coverage = combined_data.loc[:, coverage_columns].mean(axis=1)
  combined_data.insert(combined_data.shape[1], MEAN_OF_MEAN_COV_COL, mean_of_mean_coverage,
      allow_duplicates=False)
  low_coverage_regions = combined_data.where(combined_data.mean_of_mean_coverage < 2).dropna()
  low_coverage_regions = low_coverage_regions.astype(combined_data.dtypes.to_dict())
  #low_coverage_regions.loc[:, NAME_COL] = low_coverage_regions.loc[:, NAME_COL].apply(lambda x: get_first_ref(x))
  gene_locator = GeneLocator()
  low_coverage_regions.loc[:, NAME_COL] = low_coverage_regions.apply(
      lambda x: gene_locator.get_gene(x.chrom, x.start), axis=1)
  low_coverage_regions.loc[:, (join_columns + [MEAN_OF_MEAN_COV_COL])].to_csv(
      "SE51038a_2_low_coverage.csv", index=False)


if __name__ == "__main__":
  main()
