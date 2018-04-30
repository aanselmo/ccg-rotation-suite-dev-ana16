#!/usr/bin/env python
# alt shebang for cluster: !/apps/python-2.7.10/bin/python
# version 1.0
# tested with python 2.7.13
# sorry :(

__author__ = "Bruce Wollison"
__date__ = "April, 2018"

'''
    Description: This script will generate a breakmer report for a given project
'''

#######################
### library imports ###
#######################

import sys, os, argparse, shutil, pd
sys.path.append('../sql')
sys.path.append('../excel')
import db_connect as db
import excel_write as xls

##########################################
### globals, parsing and configuration ###
##########################################
# 'id' in the riker DB that points to the current breakmer run tool
# likely won't change often or ever
DEFAULT_BREAKMER_TOOL_ID = '13'

### handle command line interface with argparse ###
parser = argparse.ArgumentParser(description='This script generates reports for BreaKmer based on '\
                                             'the given project code.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--project', '-p', help='Project code to generate the report off of',
                    type=str, required=True)
parser.add_argument('--ignore', '-i', help='Creates report regardless of what errors would otherwise'\
                                           ' cause the script to die.',
                    type=bool, default=False)
parser.add_argument('--version', action='version', version='Breakmer_report.py ver 1.0')
args = parser.parse_args()

################################
### mysql queires we'll need ###
################################

# get the project id from the project code
def projectIDFromProjectCode(DBconn, project_code):
    querySTMT =  "SELECT project_id FROM project WHERE project_code = %(proj_code)s"
    queryResults = DBconn.DBquery(querySTMT, {'proj_code' : project_code})
    return queryResults

# get a panda dataframe using the project ID
# returns all complete aggregate pairs for a project and their associated completed breakmer runs
# should look at the whole panda df. If there are nulls for the breakmer path that means it was
# not run or completed on those aggregate pairs, worth mentioning
def breakmerDirectoriesFromProjectID(DBconn, project_id):
    querySTMT =  "SELECT * FROM aggregate_pair ap "\ 
                 "left join analysis_tool_run atr on atr.aggregate_pair_id = ap.aggregate_pair_id "\ 
                 "where atr.analysis_tool_id = %(breakmer_id)s and atr.status_id = 1 "\ 
                 "and atr.annotation_status_id = 1 and ap.qc_status_id = 1 "\
                 "and ap.aggregate_pair_project_id = %(proj_id)s"
    queryResults = DBconn.DBquery(querySTMT, {'proj_id' : project_id, 'breakmer_id': DEFAULT_BREAKMER_TOOL_ID})
    return queryResults

########################
### it's main() time ###
########################

if __name__ == "__main__":
    print 'RobustCNV_PON_automation.py started'

    # open configuration file and read the necessary info
    # we want to get the project name so we can give meaningful names to
    # the folders we make, and we need to know where to get the normals from
    # additionally we will need be able to write PON configs
    configFile = open(args.config).readlines()
    for lines in configFile:
        if 'tumorCoveragePath' in lines:
            tumorCoveragePath = lines.replace('"', "'").split("'")[1]
        if 'normalCoveragePath' in lines:
            normalCoveragePath = lines.replace('"', "'").split("'")[1]
        if 'experiment' in lines:
            projectCode = lines.replace('"', "'").split("'")[1]

    # now we can define the folders we want to make. We want to make one folder to
    # to hold all iterations that we plan on running.
    outputFolder =  os.getcwd() + '/' + projectCode + '_RobustCNVLeaveXOut'
    tumorFolder = outputFolder + '/tumorCoverages'
    normalFolder = outputFolder + '/normalCoverages'
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
        os.mkdir(tumorFolder)
        os.mkdir(normalFolder)

    # create a log of the parameters used for this run in the folder
    logFileName = os.getcwd() + '/' + projectCode + '_RobustCNVLeaveXOut/log.txt'
    logFile = open(logFileName, 'w')
    logFile.writelines("Reading in configuration file " + args.config)

    # copy the tumors over (if applicable)
    if args.tumors:
        tumorFiles = [tumorCoveragePath + '/' + f for f in os.listdir(tumorCoveragePath)]
        for tumorFile in tumorFiles:
            shutil.copy(tumorFile, tumorCoveragePath)

    # make note of which samples are the normals and track whether they were run
    # use lists of lists, [[samp_name1, bool], [samp_name2, bool], ... [samp_nameN, bool]]
    normalFiles = [normalCoveragePath + '/' + f for f in os.listdir(normalCoveragePath)]
    normalSamples = []
    for normalFile in normalFiles:
        # check that file is a robustCNV coverage file
        if sampPattern in normalFile:
            newSampleName = normalFile.split(fileNameSeperator)[0].split('/')[-1]
            normalSamples.append(newSampleName)

    # do sanity check here: make sure there are as many normals as there are requested
    # iterations, otherwise warn and exit
    if len(normalSamples) < args.iterations:
        logFile.writelines(["exiting, only ", str(len(normalSamples)), " when requesting ",
                            str(args.iterations), " iterations"])
        print("Fewer normals than iterations, exiting. Try running again with more normals or fewer iterations")
        sys.exit()

    # define the normals for each iteration by creating a random order of indexes
    # for pulling samples from normalSamples and splitting it up
    iterationInds = range(0, len(normalSamples))
    random.shuffle(iterationInds)
    iterationInds = split_list(iterationInds, args.iterations)

    # copy the normals into the tumor folder now
    for normalFile in normalFiles:
        shutil.copy(normalFile, tumorFolder)

    logFile.write("Script finished\n")
    print 'RobustCNV_PON_automation.py finished'

# EOF