#!/usr/bin/env python
# version 1.0
# tested with python 2.7.13

__author__ = "Bruce Wollison"
__date__ = "June, 2018"


'''
    Description: This script converts Agilent design files to a Baits class object
'''

import os
import Intervals_Class as ic

# AgilentToIntervals
def AgilentToIntervals(agilentParentDir, folderToFeature=None):
    """
    bmw35 2018-07-02
    AgilentToIntervals takes a start directory and iterates over all the bed files
    in the child directories and adds them to an Interval class object
    :param agilentParentDir: string for the parent directory to search through
    :param folderToFeature:
    :return: returns Interval class object
    """

    # interval information is stored in '.bed' files
    INTERVAL_FILE_PATTERN = ".bed"

    # initiate the Interval class that we will add to
    currIntervals = ic.Intervals()

    # iterate over parent dir and find all files in the design directory
    for root, dirs, files in os.walk(agilentParentDir):
        for file in files:
            if file.endswith(INTERVAL_FILE_PATTERN):
                print((os.path.join(root, file)))

# EOF