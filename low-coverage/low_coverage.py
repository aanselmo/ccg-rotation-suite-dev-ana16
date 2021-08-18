#!/usr/local/miniconda3/bin/python

import os
import pandas as pd

PROJECT_FOLDER = "/home/nap24/Documents/Pipeline/ccgd/AN0008a/GC_AN0008a_SE51038a_20210817/AN0008a"
CHROM_COL = "chrom"
START_COL = "start"
END_COL = "end"
LENGTH_COL = "length"
NAME_COL = "name"
GC_CONTENT_COL = "%gc"
join_columns = [CHROM_COL, START_COL, END_COL, LENGTH_COL, NAME_COL, GC_CONTENT_COL]
MEAN_COV_COL = "mean_coverage"

combined_data = pd.DataFrame()
for index, filename in enumerate(os.listdir(PROJECT_FOLDER)):
  interval_data = pd.read_csv(os.path.join(PROJECT_FOLDER, filename), sep="\t")
  if index == 0:
    combined_data = interval_data.loc[:, join_columns] 
  interval_data.loc[:, join_columns.append(MEAN_COV_COL)]

