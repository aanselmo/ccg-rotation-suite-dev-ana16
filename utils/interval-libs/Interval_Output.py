#!/usr/bin/env python
# version 1.0
# tested with python 2.7.13

__author__ = "Bruce Wollison"
__date__ = "June, 2018"


'''
    Description: This script contains functions that will take in Baits class and write the output
    to interval lists in an organized way
'''

##From the Make_Intervals.py script. Refactoring into here
def outputIntervalFiles(baitset_output_folder):
    # check to see if the bait set directory already exists in the output
    # directory, warn if it exists
    if not os.path.exists(baitset_output_folder):
        os.mkdir(baitset_output_folder)
    else:
        logger.warn("Bait set output folder " + baitset_output_folder + " already exists")

#EOF