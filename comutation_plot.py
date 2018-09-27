#!/usr/local/miniconda3/bin/python
"""
This script has the purpose of creating a co-mutation plot. See the following URL for reference:
https://iovs.arvojournals.org/article.aspx?articleid=2650631. There will be genes on the y-axis and sample names on
the y-axis.

Input data is a CSV that is formatted as the following:
Sample_name, Gene, Mutation_Type
sample1, KRAS, missense
sample1, KRAS, synonymous
sample2, TP53, missense

Set priority.

Set colors.

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

