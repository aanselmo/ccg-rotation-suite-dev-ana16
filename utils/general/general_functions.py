#!/usr/bin/env python
# version 1.0

__author__ = "Bruce Wollison"
__date__ = "Dec, 2017"

'''
    Description: This script contains functions that we will use regularly
    in multiple scripts. These functions 
'''

# Purpose: This function will alter the paths we get for files from the DB
# to the appropriate relative path based on the config file.
# Design: More generally it just looks for the last common substring in str2
# deliminated by '/' characters and replaces the first hald of str1 and
# the second half of str2
# EX: str1 "/mounted/Project/"
# str2 "/data/Project/file.txt"
# output: "/mounted/Project/file1.txt"
def replaceWCommonSubStr(str1, str2):
    path1 = str1.split("/")
    path2 = str2.split("/")
    returnStr = ""
    commonSubStr = ""
    # search for a common substring
    for substr2 in path2:
        for substr1 in path1[::-1]:
            if substr1 == substr2:
                if commonSubStr == "":
                    commonSubStr = substr1
    # if we found a common substring join on it
    if commonSubStr != "":
        returnPathPart1 = "/".join(path1[:path1.index(commonSubStr)])
        returnPathPart2 = "/".join(path2[path2.index(commonSubStr):])
        returnStr = "/".join([returnPathPart1, returnPathPart2])
    # return str2 if no common substr exists
    else:
        returnStr = str2
    return returnStr


# Purpose: replaces rbind() from r. Pandas seems to be inflexible to how I like to use
# dataframes. I would like to be able to take N columns from DF1 and bind them
# to N columns from DF2 even if the column names are dissimilar.
# use cases include cases where we have tumor/normal data in 2 tables w diff columns
# names that have the same data structure but we want
# them in one single dataframe in order to iterate easily
# default behavior should use the columns from dataframe1 if no new columns are specified
def rbind(dataframe1, dataframe2, newCols=None):
    if newCols == None:
        newCols = dataframe1.columns
    dataframe1.columns = newCols
    dataframe2.columns = newCols
    return dataframe1.append(dataframe2)

# EOF