#!/usr/bin/env python
# version 1.0
# tested with python 2.7.13

__author__ = "Bruce Wollison"
__date__ = "June, 2018"

'''
    Description: This script defines the Intervals class for manipulating intervals and creating various interval 
    lists out of
'''

# library imports
import pybedtools

class Interval:
    """
    bmw35 06-26-2018
    class Interval abstracts groups of genomic intervals that define bait sets
    Ideally we want to store all information in a bait set into one object that
    has methods for the various different formats we need. Should be created and called
    by Vendor_to_Interval.py functions, and manipulated for outputs by Interval_Output.py.
    
    Methods should be able to:
    -merge with other Interval objects (useful for building)
    -intersect with other Interval objects (ex: return the intervals )
    -expand or contract intervals (useful for adding overhangs to baits)
    -redefine features of intervals (may be necessary if we leave interval features undefined)
    -combine different features
    -combine different coverage_types
    -return raw intervals based on how they are stored (may have overlaps)
    -return intervals without overlaps
    
    We will mostly be using pybedtools which handles most of this functionality laready
    
    Bait sets have a Designed region and a Covered region. Designed is what we want the
    vendor to bait for and Covered is what they actually baited for. Designed may or may not exist.
    We expect Covered to always exist. While a vendor may not tell us what they actually covered
    and only give us back the Designed regions we will consider this the 'Covered' for the sake of regions
    """

    # globals for this class

    # constructor
    # initialize the Interval class with nulls for Designed and Covered baits
    def __init__(self):
        self.Designed = None
        self.Covered = None

    # return Not Covered intervals
    # the difference between what we wanted intervals for and what we got intervals for
    def returnNotCovered(self):
        if self.Designed && self.Covered:
            return self.Covered - self.Designed

    # return combined intervals
    # if there is a Designed and Covered file then merge them and return the intervals
    def combineDesignedAndCovered(self):
        if self.Designed && self.Covered:
            return self.Covered + self.Designed
        else:
            # Desgined doesn't exist, report error
            ### DESIGN DECISION, do we return the covered only and just call this combine?
            ### actual use case may work better that way
            return "Error"

    # merge features from another bed formated file
    # maps the feature (designated by featureNum) from the provided featureBed and sets those features
    # in the current interval object.
    # useful for using a reference file to pull out Exonic/Intronic regions for feature mapping
    # or for getting gene names from another ref
    def mergeFeatures(self, featureBed, featureNum):

    # expands intervals by a given size (size in bp)
    def expandIntervals(self, expandSizeBP):

    # contracts intervals by a given size (size in bp)
    def contractIntervals(self, contractSizeBP):

    # add regions to design interval file
    def addDesignIntervals(self, intervalsToAdd):
        self.Designed = self.Designed + pybedtools.BedTool(intervalsToAdd)

    # add regions to design interval file
    def addCoveredIntervals(self, intervalsToAdd):
        self.Covered = self.Covered + pybedtools.BedTool(intervalsToAdd)

# EOF