#!/usr/bin/env python
# alt shebang for cluster: !/apps/python-2.7.10/bin/python
# version 1.0
# tested with python 2.7.13
# sorry :(

__author__ = "Bruce Wollison"
__date__ = "May, 2018"

'''
    Description: This script will generate a breakmer report for a given project
    
    Workflow is broken into 3 steps:
    1) handle inputs, query the database, get all the information we need to verify that the
        run is finished and we have what we need
    2) create symlinks from the output of the breakmer run so we can run the process_cohort.py script
        and run it to generate the report data files
    3) generate the excel report based on the data files. 
'''

#######################
### library imports ###
#######################

import sys, os, argparse, yaml, shutil, fnmatch
import pandas as pd

# add special modules to the path
sys.path.append('utils/sql')
sys.path.append('utils/excel')
sys.path.append('utils/general')
import dbConnect as db
import excel_write as xls
import general_functions as genFun
import common_queries as com_queries

##########################################
### globals, parsing and configuration ###
##########################################
# 'id' in the riker DB that points to the current breakmer run tool
# likely won't change often or ever
DEFAULT_BREAKMER_TOOL_ID = '13'
# default config file needed for connecting to the DB info
DEFAULT_CONFIGURATION_FILE = 'ccgd_rotation_scripts.config'
DEFAULT_OUTPUT_DIRECTORY = os.getcwd() + '/BreaKmer_Reports'
REARRANGE_FILE_PATTERN = "_Breakmer_0_0_6_rearrangement_svs.out"
INDEL_FILE_PATTERN = "_Breakmer_0_0_6_indel_svs.out"
yes = {'yes', 'y', 'ye', ''}
no = {'no', 'n'}

### handle command line interface with argparse ###
parser = argparse.ArgumentParser(description='This script generates reports for BreaKmer based on '\
                                             'the given project code.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--project_code', '-p', help='Project code to generate the report off of',
                    type=str, required=True)
parser.add_argument('--ignore', '-i', help='Creates report regardless of what errors would otherwise'\
                                           ' cause the script to die.',
                    type=bool, default=False)
parser.add_argument('--config', '-c', help='configuration file to use',
                    default=DEFAULT_CONFIGURATION_FILE, type=str)
parser.add_argument('--version', action='version', version='Breakmer_report.py ver 1.0')
args = parser.parse_args()

### import configurations from configuration file ###
print "Reading in configuration file " + args.config
try:
    with open(args.config, 'r') as ymlfile:
        config = yaml.load(ymlfile)
except IOError:
    # whoa whats happening bros? File might not exist, or some other IO error
    print "Error reading configuration file!"
    sys.exit()
DBconfig = config['rikerDB']
genConfig = config['general']
projectCode = args.project_code

#####################################
## small script specific functions ##
#####################################


################################
### mysql queires we'll need ###
################################

# get a panda dataframe using the project ID
# returns all complete aggregate pairs for a project and their associated completed breakmer runs
# should look at the whole panda df. If there are nulls for the breakmer path that means it was
# not run or completed on those aggregate pairs, worth mentioning
def breaKmerDirectoriesFromProjectID(DBconn, project_id, breakmer_tool_id):
    querySTMT =  "SELECT ap.aggregate_pair_id, ap.aggregate_pair_name, atr.* FROM aggregate_pair ap " +\
                 "left join analysis_tool_run atr on atr.aggregate_pair_id = ap.aggregate_pair_id " + \
                 "and atr.annotation_status_id = 1 and " +\
                 "atr.analysis_tool_id = %(breakmer_id)s and atr.status_id = 1 " + \
                 "where ap.qc_status_id = 1 and ap.aggregate_pair_project_id = %(proj_id)s"
    queryResults = DBconn.DBquery(querySTMT, {'proj_id' : project_id, 'breakmer_id': breakmer_tool_id})
    return queryResults

# Need to query DB and get back the project code with the description
# ex: SCB0020a gives back SCB0020a_RadioGenomicsB1
def projectDescriptionQuery(DBconn, project_code):
    querySTMT = "SELECT concat(project_code, '_', project_name) as projectDescription FROM project " +\
                "where project_code = %(proj_code)s;"
    queryResults = DBconn.DBquery(querySTMT, {'proj_code': project_code})
    return queryResults


########################
### it's main() time ###
########################

if __name__ == "__main__":
    print 'BreaKmer_report.py started'

    ################# WORKFLOW PT 1 ################################
    ### Check the DB to see what we need to do with this project ###

    # make a DB connection
    dbConn = db.DBconnect(DBconfig)

    # get the necessary project/sample information
    # start by getting the project_id from the project code
    projectIDQueryResults = com_queries.projectIDFromProjectCode(dbConn, projectCode)

    # make sure something was returned, AKA project code exists in DB
    if projectIDQueryResults.empty:
        print "Error, project code " + projectCode + " does not exist in " +\
              DBconfig['schema'] + "!!! Exiting..."
        sys.exit()
    project_id = projectIDQueryResults.project_id[0]

    # get project description
    projectDescriptionQueryResults = projectDescriptionQuery(dbConn, projectCode)
    projectDescription = projectDescriptionQueryResults['projectDescription'][0]

    # now query for the main breakmer files table
    breaKmerDirectoriesQueryResults = breaKmerDirectoriesFromProjectID(dbConn, project_id, DEFAULT_BREAKMER_TOOL_ID)

    # make sure something was returned, AKA aggregate pairs exist exists in DB
    if breaKmerDirectoriesQueryResults.empty:
        print "Error, no aggregate pairs exist for project  " + projectCode + " in " +\
              DBconfig['schema'] + "!!! Exiting..."
        sys.exit()

    # ag pairs exist, but the breakmer run may not be successful, or finished, so lets make sure:
    unFinishedAgPairs = breaKmerDirectoriesQueryResults[ \
        (breaKmerDirectoriesQueryResults['status_id'] != 1) | \
        (breaKmerDirectoriesQueryResults['annotation_status_id'] != 1)]
    # if subsetting on unfinished breakmer runs returns nothing, then all breakmer runs are done
    # kinda a round-a-bout way of doing it really
    if unFinishedAgPairs.empty:
        print "All aggregate pairs have finished processing."
    else:
        print "Some aggregates have not finished processing. Specifically ag pairs: " + \
              (", ").join(unFinishedAgPairs.aggregate_pair_name) +\
              "\nContinue generating report with currently finished ags? Y/N"
        # raw_input returns the empty string for "enter"
        choice = raw_input().lower()
        if choice in no:
            print "nohup ./AggregatePairAnalysisforAnalysisTools.pm -AggregatePairId " + \
                (", ").join(unFinishedAgPairs.aggregate_pair_id) +\
                " -config ./Riker.ccgd.rcsgc.config &> log/AnalysisTools_" + \
                  projectCode + "_breakmer.txt &"
            sys.exit()

    # check that all the files are available on the file system
    noLongerAvailableBreaKmerRuns = breaKmerDirectoriesQueryResults[ \
        (breaKmerDirectoriesQueryResults['data_location_id'] != 1) & \
        (breaKmerDirectoriesQueryResults['status_id'] == 1) & \
        (breaKmerDirectoriesQueryResults['annotation_status_id'] == 1)]

    if not noLongerAvailableBreaKmerRuns.empty:
        # seems like we got some backed up data we are trying to make a report from
        # warn that we are missing things on the file system we need for the full report
        print "Some breaKmer outputs have already been backed up. Specifically ag pairs: " + \
              (", ").join(noLongerAvailableBreaKmerRuns.aggregate_pair_name) + \
              "\nContinue generating report with currently finished ags? Y/N"
        # raw_input returns the empty string for "enter"
        choice = raw_input().lower()
        if choice in no:
            sys.exit()

    # important to abstract out of this step anything specific to the query
    # in the future we will likely want multiple paths through workflow pt1, one where
    # we don't hit Riker but go off user input of starting directories when we
    # are running something outside of the DB
    # additionally we are also replacing the DB name with any local paths to the isilon
    # that we define in the config.
    breaKmerPerantDirectories = breaKmerDirectoriesQueryResults[['aggregate_pair_name','output_file']]
    newBreaKmerParentDir = []
    for parentDir in breaKmerPerantDirectories['output_file']:
        newDir = genFun.replaceWCommonSubStr(genConfig['isilon_project_directory'], parentDir)
        newBreaKmerParentDir.append(newDir)
    breaKmerPerantDirectories['output_file'] = newBreaKmerParentDir

    ################# WORKFLOW PT 2 #############################
    ### replaces Ryans link_breakmer_files.rcprime.py script  ###
    ### will also create some files to simplify the workflow ####

    # create new folder for this project
    projectOutputDirectory = DEFAULT_OUTPUT_DIRECTORY + "/" + projectCode + "_breaKmerData"
    if not os.path.exists(projectOutputDirectory):
        os.mkdir(projectOutputDirectory)
    else:
        # warn that the directory already exists and that we would want to
        # overwrite everything.
        print "This script appears to already have been run on this project. " + \
              "\nDelete current BreaKmer results directory to Continue? Y/N"
        # raw_input returns the empty string for "enter"
        choice = raw_input().lower()
        if choice in yes:
            shutil.rmtree(projectOutputDirectory)
            os.mkdir(projectOutputDirectory)
        else:
            sys.exit()

    # need to make links for the following files: "*rearrangement_svs.out", "*indel_svs.out", "output/"
    # that way Ryans other script can find them on the file system and extract the needed info
    for aggPairName,directoryToLink in zip(breaKmerPerantDirectories['aggregate_pair_name'], \
                                           breaKmerPerantDirectories['output_file']):
        print "creating symlinks for aggregate pair" + aggPairName +\
            " using directory " + directoryToLink
        # for each parent dir path, create one of each symlink needed (indel, rearrangement, outputfolder)
        # first collect all the file paths we need fro the breakmer output
        breaKmerOutputDirectoryPath = directoryToLink + "/output"
        # NOTE: output directory must exist. Indel/rearrangement sys.out may not!
        if os.path.exists(breaKmerOutputDirectoryPath):
            allBreakmerFiles = os.listdir(breaKmerOutputDirectoryPath)
        else:
            print "Warning!!!! Did not find the following expected output directory:\n" + \
                  breaKmerOutputDirectoryPath + "\nCheck this before continuing"
            sys.exit()
        rearrangementOutPath = [returnPath for returnPath in allBreakmerFiles if REARRANGE_FILE_PATTERN in returnPath]
        indelOutPath = [returnPath for returnPath in allBreakmerFiles if INDEL_FILE_PATTERN in returnPath]

        # create sym links
        outputSymLinkDest = projectOutputDirectory + "/" + aggPairName + "_breakmer_out"
        os.symlink(breaKmerOutputDirectoryPath, outputSymLinkDest)

        # for these we assume there is only 1 of each file type in each directory
        # this should be true for the breaKmer generated output
        # but could be technically wrong if there is nefarious behavior going on
        if len(indelOutPath) > 0:
            indelSymLinkDest = projectOutputDirectory + "/" + aggPairName + "_Breakmer_0_0_6_indel_svs.out"
            os.symlink(("/").join([breaKmerOutputDirectoryPath, indelOutPath[0]]), indelSymLinkDest)

        if len(rearrangementOutPath) > 0:
            rearrangementSymLinkDest = projectOutputDirectory + "/" + aggPairName +\
                                       "_Breakmer_0_0_6_rearrangement_svs.out"
            os.symlink(("/").join([breaKmerOutputDirectoryPath, rearrangementOutPath[0]]), rearrangementSymLinkDest)

    # create the filters file, normal file, and tumors file. These are neede for running the
    # process_cohort_results.py script. After we do this, lets run that script to generate the desired
    # output
    filterListFile = projectOutputDirectory + "/" +  "filters.txt"
    open(filterListFile, 'a').close()
    # Note: we have decided in CCGD not to make use of normal samples for filtering
    # breaKmer calls, so all samples are treated as tumors. However to wrap to the process_cohort/py
    # script we still need to create the empty file
    normalSampleList = projectOutputDirectory + "/" +  projectCode +  "_normals.txt"
    open(normalSampleList, 'a').close()
    # we need to write the tumor file list, that is usually done by the link files script
    # it is just a list of aggregate_pair_names and paths
    tumorSampleList = projectOutputDirectory +  "/" +  projectCode + "_tumors.txt"
    with open(tumorSampleList, 'w') as tumorSampleListFile:
        for aggPair in breaKmerPerantDirectories['aggregate_pair_name']:
            tumorSampleListFile.write(("\t").join([projectDescription, aggPair]))
            tumorSampleListFile.write("\n")

    ################# WORKFLOW PT 3 #############################
    ### we have all the data from the runs now, lets produce  ###
    ################## the report of the data ###################


# EOF