#!/usr/bin/env python
# alt shebang for cluster: !/apps/python-2.7.10/bin/python
# version 1.0
# tested with python 2.7.13
# sorry :(

__author__ = "Bruce Wollison"
__date__ = "April, 2018"

'''
    Description: This script will run RobustCNV multiple times on by employing a leave-X-out
    strategy on the normals directory. This should be used for new bait sets where we are unsure
    about the quality of the normals.
    
    Currently does not support replacing samples each iterations (bootstrapping).
'''

#######################
### library imports ###
#######################

import sys, os, argparse, shutil, random

##########################################
### globals, parsing and configuration ###
##########################################
DEFAULT_ROBUST_LOCATION = "/ifs/rcgroups/ccgd/bmw35/RobustCNV_projects/robustCNV-GATK-MUTECT_Kluge.R"

### handle command line interface with argparse ###
parser = argparse.ArgumentParser(description='This script runs RobustCNV multiple times using a '\
                                             'Leave-X-Out approach on the normals for verifying a PON.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--RobustCNV', '-r', help='location of RobustCNV script to run',
                    type=str, default=DEFAULT_ROBUST_LOCATION)
parser.add_argument('--config', '-c', help='Robust configuration file to use (must exist)',
                    type=str, required=True)
parser.add_argument('--iterations', '-i', help='Specifies how many iterations of the PON testing to perform.',
                    type=int, default=5)
parser.add_argument('--tumors', '-t', help='if `True` then tumor samples are included in the test '\
                                           'cohort each iteration',
                    type=bool, default=False)
parser.add_argument('--version', action='version', version='RobustCNV_PON_automation.py ver 1.0')
args = parser.parse_args()

# filenameing pattern for the RobustCNV coverages
sampPattern = ".sample_interval_summary"
# first substring the Riker pipeline adds to the file after the
# aggregate name, this is used to seperate the aggregate name from
# the file name
fileNameSeperator = ".dedup"

####################################################
### define any 'small' script-specific functions ###
####################################################

# splitting lists into lists of lists
# taking from stack overflow:
# https://stackoverflow.com/questions/752308/split-list-into-smaller-lists
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
             for i in range(wanted_parts) ]

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

    # loop over however many leave-X-out iterations we want
    for currIterationNum in range(0,args.iterations):

        # name for new iteration
        currntIterationFolder = projectCode + '_PONIteration_' + str(currIterationNum)

        # get the new set of samples to iterate with
        # and move them into the normals folder from the tumor folder
        currPON = [normalSamples[ind] for ind in iterationInds[currIterationNum]]
        logFile.writelines('\nStarting Robust iteration ' + str(currIterationNum) + ' using the following as the PON:\n')
        logFile.writelines('\n'.join(currPON))
        for PONmember in currPON:
            PONtoMove = [f for f in os.listdir(tumorFolder) if PONmember in f]
            for PONmember in PONtoMove:
                oldFileDest = tumorFolder + '/' + PONmember
                newFileDest = normalFolder + '/' + PONmember
                shutil.move(oldFileDest, newFileDest)

        # create new config to run
        # we're just going to copy the old config line by line and replace the lines we need altered
        currIterationConfigFile = outputFolder + '/robustCNV_PON_Test_' + projectCode + '_iteration' \
                                  + str(currIterationNum) + '.txt'
        with open(args.config, 'r') as originalConfig, open(currIterationConfigFile, 'w') as newConfig:
            for line in originalConfig:
                if 'tumorCoveragePath' in line:
                    newConfig.write("tumorCoveragePath = '" + tumorFolder + "'\n")
                elif 'normalCoveragePath' in line:
                     newConfig.write("normalCoveragePath = '" + normalFolder + "'\n")
                elif 'experiment' in line:
                     newConfig.write("experiment = '" + currntIterationFolder + "'\n")
                else:
                    newConfig.write(line)

        # run RobustCNV again using the new config on the new PON
        robustRunCommand = 'Rscript ' + args.RobustCNV + ' ' + currIterationConfigFile
        logFile.writelines('\nRunning robust command:' + robustRunCommand)
        os.system(robustRunCommand)

        # after it returns move it to the RobustCNVLeaveXOut folder
        if '/' in args.RobustCNV:
            robustParentFolder = '/'.join(args.RobustCNV.split('/')[:-1]) + '/'
        else:
            robustParentFolder = ''
        RobustIterationOutput = robustParentFolder + currntIterationFolder + '_RobustCNV_results'
        shutil.move(RobustIterationOutput, outputFolder)

        # move anything in the normals folder to the tumor folder
        for oldPONfile in os.listdir(normalFolder):
            oldFileDest = normalFolder + '/' + oldPONfile
            newFileDest = tumorFolder + '/' + oldPONfile
            shutil.move(oldFileDest, newFileDest)

    logFile.write("Script finished\n")
    print 'RobustCNV_PON_automation.py finished'

# EOF