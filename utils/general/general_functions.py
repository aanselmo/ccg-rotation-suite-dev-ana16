#!/usr/bin/env python
# version 1.0

__author__ = "Bruce Wollison"
__date__ = "Dec, 2017"

'''
    Description: This script contains functions that we will use regularly
    in multiple scripts. These functions 
'''
import logging, sys
from logging.handlers import RotatingFileHandler

def getLogger(logFileName):
    """
    bmw35 12-29-2017
    Logger should be same format for all scripts, and there isn't a reason to
    rewrite the same basic configuration everytime, so we put it here
    we want to write to a log file and stdout, while rotating the log file when necessary
    :param logFileName: string of the file name to write the log file to
    :return: logging object
    """
    # get logging object
    logger = logging.getLogger('Logger')
    # add rotating file handler if size gets too big
    MAX_LOGFILE_DISK_SIZE = 50000000
    fileHandle = RotatingFileHandler(filename=logFileName,
                                     maxBytes=MAX_LOGFILE_DISK_SIZE, backupCount=5)
    formatter = logging.Formatter('%(asctime)s: %(module)s - %(levelname)s: %(message)s')
    fileHandle.setFormatter(formatter)
    logger.addHandler(fileHandle)
    # set and format the stream handler
    logger.addHandler(logging.StreamHandler(sys.stdout))
    # set logging level
    logger.setLevel(logging.INFO)
    return logger

def replaceWCommonSubStr(str1, str2):
    """
    bmw35 12-29-2017
    Purpose: alter the paths we get for files from the DB
    to the appropriate relative path based on the config file.
    Design: More generally it just looks for the last common substring in str2
    that is also in common with str1 (with substrings deliminated by '/' characters)
    and returns a string with the first half of str1 and the second half of str2
    EX: str1 "/mounted/Project/"
    str2 "/data/Project/file.txt"
    output: "/mounted/Project/file1.txt"
    :param str1: first '/' delimited string
    :param str2: second '/' delimited string
    :return: string of the left str1 chars and the right of str2 chars with the common
      substr in the middle
    """
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


def rbind(dataframe1, dataframe2, newCols=None):
    """
    bmw35 12-29-2017
    Purpose: replaces rbind() from r. Pandas seems to be inflexible to how I like to use
    dataframes. I would like to be able to take N columns from DF1 and bind them
    to N columns from DF2 even if the column names are dissimilar.
    use cases include cases where we have tumor/normal data in 2 tables w diff columns
    names that have the same data structure but we want
    them in one single dataframe in order to iterate easily
    default behavior should use the columns from dataframe1 if no new columns are specified
    :param dataframe1: first dataframe. Used as colnames if newCols is not specified
    :param dataframe2: second dataframe to append to the bottom of dataframe1
    :param newCols: optional. list of colnames to give the returned dataframe
    :return: dataframe of dataframe1 and dataframe2 bound together
    """
    if newCols == None:
        newCols = dataframe1.columns
    dataframe1.columns = newCols
    dataframe2.columns = newCols
    return dataframe1.append(dataframe2)

# EOF
