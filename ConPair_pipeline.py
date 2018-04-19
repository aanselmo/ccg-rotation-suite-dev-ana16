#!/usr/bin/env python
# version 1.0
# tested with python 2.7.13

__author__ = "Bruce Wollison"
__date__ = "Nov, 2017"

'''
    Description: This will run ConPair for a specified project code
    This will require (aka workflow outline):
        1) hitting the RikerDB to find the samples associated
        with the specified project code, get the paths for the sample pairs
        to test
    
        2) running the Conpair pipeline - this includes generating the pileups 
        for each sample, concordance for sample pairs, and contamination for
        sample pairs (recommended to use suprocess module instead of os.system())
        
        3) organizing the concordance and contamination results for all sample
        pairs and writing them to an excel file
    
    Notes:
        ConPair only works on tumor/normal pairs
        ConPair takes a while to run - specifically the pileup and contamination steps
            because of this we should be verbose. It could take hours and the user
            should know the program is processing
        This is version 1.0, the output excel file may be subject to modification as we
            get a better idea of what we want in a report
            
    How to run (minimal command line):
        python ConPair_pipeline.py -p [project_code]
'''

#######################
### library imports ###
#######################

import pandas as pd
import sys, imp, os, argparse, yaml

# add special modules to the path
sys.path.append('utils/sql')
sys.path.append('utils/excel')
sys.path.append('utils/general')
import db_connect as db
import excel_write as xls
import common_queries as com_queries
import general_functions as genFun

##########################################
### globals, parsing and configuration ###
##########################################
DEFAULT_LOGGER_FILE = "ConPairLogFile.log"
DEFAULT_CONFIGURATION_FILE = 'ccgd_rotation_scripts.config'
DEFAULT_OUTPUT_DIRECTORY = os.getcwd() + '/ConPair_Reports'

# set up logger
logger = genFun.getLogger(DEFAULT_LOGGER_FILE)
logger.info("Starting new ConPair Pipeline run with command: "  + " ".join(sys.argv))

### handle command line interface with argparse ###
parser = argparse.ArgumentParser(description='Run ConPair pipeline and create report.',
                                 epilog='Warning: This script will often take a long time to complete',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--project_code', '-p', help='project code of project to run the ConPair pipeline on',
                    type=str, required=True)
parser.add_argument('--output_folder', '-o', help='folder to output the results to',
                    default=DEFAULT_OUTPUT_DIRECTORY, type=str)
parser.add_argument('--config', '-c', help='configuration file to use',
                    default=DEFAULT_CONFIGURATION_FILE, type=str)
parser.add_argument('--version', action='version', version='ConPair_pipeline.py ver 1.0')
args = parser.parse_args()

### import configurations from configuration file ###
logger.info( "Reading in configuration file " + args.config)
try:
    with open(args.config, 'r') as ymlfile:
        config = yaml.load(ymlfile)
except IOError:
    # whoa whats happening bros? File might not exist, or some other IO error
    logger.error("Error reading configuration file!")
    sys.exit()

# parse the needed variables from the config file into readable variables
genConfig = config['general']
DBconfig = config['rikerDB']
logger.info( "configuration file read in successfully" )

# add the conpair module directory to the pythonpath variable
# then import the ConPair scripts
ConPairConfig = config['conpair']
sys.path.append(ConPairConfig['conpair_modules_path'])
sys.path.append(ConPairConfig['conpair_directory']+'/scripts')
import run_gatk_pileup_for_sample as run_pileup
import verify_concordance as run_concordance
import estimate_tumor_normal_contamination as run_contamination

# lets make a few more variable names that will be useful later
# starting with output folders of the ConPair analysis
project_folder = args.output_folder + "/" + args.project_code
pileup_output_dir = project_folder + "/pileups"
concordance_output_dir = project_folder + "/concordance"
contamination_output_dir = project_folder + "/contamination"
conpair_report_fileName = project_folder + "/" + args.project_code + "_report.xlsx"

##################################################
### define any small script-specific functions ###
##################################################

def concordance_parser(concordanceDF, tumorSample, normalSample):
    """
    bmw35 12-29-2017
    concordance output by ConPair is specific
    so we need to parse it in very specific ways
    :param concordanceDF: dataframe with only 1 column for each row produced by
      ConPair concordance output
    :param tumorSample: string for the tumor sample name
    :param normalSample: string for the normal sample name
    :return: dataframe with the relevant concordance information parsed into a usable format
    """
    concordanceValue = concordanceDF[0][0].split(': ')[1]
    markersUsed = concordanceDF[0][1].split('Based on ')[1].split(' markers ')[0]
    markerThreshold = concordanceDF[0][1].split(': ')[1].split(' reads')[0]
    minMapQuality = concordanceDF[0][2].split(': ')[1]
    minBaseQuality = concordanceDF[0][3].split(': ')[1]
    return pd.DataFrame({'tumor sample' : [tumorSample],
                         'paired normal': [normalSample],
                         'concordance' : [concordanceValue],
                         'markers used' : [markersUsed],
                         'min coverage for markers' : [markerThreshold],
                         'min map quality' : [minMapQuality],
                         'min base quality' : [minBaseQuality]})

def contamination_parser(contaminationDF, tumorSample, normalSample):
    """
    bmw35 12-29-2017
    contamination table is similar to the concordance table, just with less
    specific data to parse
    :param contaminationDF: dataframe with only 1 column for each row produced by
      ConPair contamination output
    :param tumorSample: string for the tumor sample name
    :param normalSample: string for the normal sample name
    :return: 
    """
    # contamination tool doesn't report out contamination values
    # for tumors when there are too few loci with adequate coverage,
    # here we assume that if there are contamination values the normal
    # returns, since the normal loci are used for determining the tumor
    # loci, therefore the tumor loci should be a subset of the normals
    normalSampleContamination = contaminationDF[0][0].split(': ')[1]
    # if there isn't 2 rows, there must be no tumor contamination information
    FULL_CONTAMINATION_DF_ROWS = 2
    if len(contaminationDF.index) == FULL_CONTAMINATION_DF_ROWS:
        tumorSampleContamination = contaminationDF[0][1].split(': ')[1]
    else:
        tumorSampleContamination = 'NA'
    return pd.DataFrame({'tumor sample' : [tumorSample],
                                       'paired normal' : [normalSample],
                                       'tumor sample contamination' : [tumorSampleContamination],
                                       'normal sample contamination' : [normalSampleContamination]})

########################
### it's main() time ###
########################

if __name__ == "__main__":
    ###### WORKFLOW PART 1 ######
    ### start with flight check and gathering the necessary information
    ### from the DB
    logger.info("Checking database to see what samples need to be run through ConPair")
    # set up DB connection
    dbConn = db.DBconnect(DBconfig)

    # first check the project exists in the DB and that some samples
    projectIDQueryResults = com_queries.projectIDFromProjectCode(dbConn, args.project_code)
    # make sure something was returned, AKA project code exists in DB
    if projectIDQueryResults.empty:
        logger.error("Error, project code " + args.project_code + " does not exist in "
                 + DBconfig['schema'] + "!!! Exiting...")
        sys.exit()
    project_id = projectIDQueryResults.project_id[0]
    # make sure samples exist with pairs (tumor/normals) per individual
    projectSamplePairsQuery = com_queries.projectSamplePairsFromProjectID(dbConn, project_id)
    if projectSamplePairsQuery.empty:
        logger.error("Error, project " + args.project_code + " has no sample pairs to run ConPair on according to "
                 + DBconfig['schema'] + "!!! Exiting...")
        sys.exit()
    # get the aggregation infct scormation so we know where the bam files are
    aggregatePairTableQuery = com_queries.aggregatePairTableFromProjectID(dbConn, project_id)
    if aggregatePairTableQuery.empty:
        logger.error("Error, project " + args.project_code + " has no aggregate pairs to run ConPair on according to "
                 + DBconfig['schema'] + "!!! Exiting...")
        sys.exit()

    # check that the output directory exists, and create the necessary subdirs
    if os.path.exists(args.output_folder):
        # make the necessary subdirs if needed
        for neededDirectory in [project_folder, pileup_output_dir, contamination_output_dir, concordance_output_dir]:
            if not os.path.exists(neededDirectory):
                os.mkdir(neededDirectory)
    else:
        # output folder doesn't exist, user needs to create it or pick a diff output folder
        logger.error("Error, output directory " + args.output_folder
                 + " does not exist!!!\nYou need to make sure it exists before running.")
        sys.exit()

    ### needed directories exist, and we
    ### now need to gather all the information
    ### to run the commands
    # reorganize sample_id and bam_file_path info
    aggPairs = genFun.rbind(aggregatePairTableQuery[['control_sample_id','control_sample_path']],
          aggregatePairTableQuery[['case_sample_id','case_sample_path']], ['sample_id','bam_file_path']).dropna()

    # iterate over list altering the DB given path to the bam file
    # to the actual path on this system
    for bamPath in aggPairs.bam_file_path:
        newPath = genFun.replaceWCommonSubStr(genConfig['isilon_project_directory'], bamPath)
        aggPairs.bam_file_path.replace(bamPath, newPath, inplace=True)

    # reorganize tumor/normal samples into 2 columns of sample/sample_id
    allSamples = genFun.rbind(projectSamplePairsQuery[['normal_sample', 'normal_sample_id']],
                               projectSamplePairsQuery[['tumor_sample', 'tumor_sample_id']],
                               ['sample_name','sample_id']).dropna()

    # merge samplePair info with the sample ID info
    allSamples = allSamples.merge(aggPairs, on="sample_id", how="left")

    # lets do 1 final data check here
    # there is a conceivable case where the sample table indicates pairs
    # of samples, but the aggregate table only has some of those samples
    # currently available (ex: all tumors aggregated, but none of the normals)
    # so lets just check the merge
    if allSamples['bam_file_path'].isnull().all():
        logger.error("Error, project " + args.project_code + " has no tumor/normal pairs with bam files according to "
                 + DBconfig['schema'] + "!!! Exiting...")
        sys.exit()

    # create the pileup output file path
    pileup_output_file = pileup_output_dir + "/" + allSamples.sample_name + "_ConPairPileup.pu"
    allSamples['pileup_output_file'] = pileup_output_file

    # now that we have the pileup output filenames we can generate the
    # table for all pairs and populate it with the pileup paths
    samplePairs = projectSamplePairsQuery
    samplePairs['normal_pileup_file'] = samplePairs.merge(allSamples[['sample_id','pileup_output_file']],
                                    how='left', right_on='sample_id', left_on='normal_sample_id')['pileup_output_file']
    samplePairs['tumor_pileup_file'] = samplePairs.merge(allSamples[['sample_id','pileup_output_file']],
                                    how='left', right_on='sample_id', left_on='tumor_sample_id')['pileup_output_file']
    # add concordance/contamination files to the table for easy iterating over
    # samplePairs table to generate the commands
    samplePairs['concordance_output_file'] = concordance_output_dir + "/" + samplePairs['tumor_sample'] + \
                                             "_v_" + samplePairs['normal_sample'] + "_ConPairConcordance.tsv"
    samplePairs['contamination_output_file'] = contamination_output_dir + "/" + samplePairs['tumor_sample'] + \
                                             "_v_" + samplePairs['normal_sample'] + "_ConPairContamination.tsv"
    # column to track what we will report out, set True everytime we verify that
    # we can get a concordance/contamination because the tumor/normal pileups exist
    samplePairs['report_out'] = False

    ###### WORKFLOW PART 2 ######
    ### invoking ConPair commands
    ### starting with pileups
    logger.info("Running ConPair pipeline, starting with pileups")
    # iterate over allSamples table, since that contains all samples
    # we need to analyze then run the appropriate pileup
    for index, sampleRow in allSamples.iterrows():
        pileupOutputFile = sampleRow['pileup_output_file']
        bamFilePath = sampleRow['bam_file_path']
        # check that the bam file exists, it may be inaccessible
        # this should not be the case, however we frequently do hacks to the riker
        # database so this is a fine sanity check
        if pd.isnull(bamFilePath):
            logger.warn("Skipping generating pileup for file " + pileupOutputFile + " because bam path is NaN")
            # mark the pileup output to indicate that pileup failed
            allSamples['index', 'pileup_output_file'] = ""
        elif os.path.isfile(bamFilePath):
            # if we don't want to overwrite the an existing pileup we should also
            # check to make sure it doesn't already exist
            if os.path.isfile(pileupOutputFile) and ConPairConfig['rewrite_previous_files'] == False:
                logger.warn("Skipping generating pileup for file " + pileupOutputFile + " because it already exists")
            else:
                run_pileup.run_gatk_pileup_for_sample(ConPairConfig['conpair_directory'],
                                                      ConPairConfig['gatk_jar_path'], bamFilePath, pileupOutputFile,
                                                      ConPairConfig['java_path'])
        else:
            logger.warn("WARNING! File " + bamFilePath + " does not exist!!!")
    # TODO: Currently no check to see that any pairs have pileups in order to move forward
    # the program gives a lot of warnings but will still generate an empty report.

    # next we need to run the concordance/contamination tools for all sample pairs
    for index, samplePair in samplePairs.iterrows():
        contaminationOutFile = samplePair['contamination_output_file']
        concordanceOutFile = samplePair['concordance_output_file']
        normalPileupFile = samplePair['normal_pileup_file']
        tumorPileupFile = samplePair['tumor_pileup_file']
        tumorSampleName = samplePair['tumor_sample']
        # check that the pileup file exists, it may be inaccessible
        if os.path.isfile(normalPileupFile) and os.path.isfile(tumorPileupFile):
            # if we don't want to overwrite the an existing pileup we should also
            # check to make sure it doesn't already exist
            # run ConPair concordance tool
            if os.path.isfile(concordanceOutFile) and ConPairConfig['rewrite_previous_files'] == False:
                logger.warn("Skipping generating concordance for file " + concordanceOutFile + " because it already exists")
            else:
                run_concordance.verify_concordance(ConPairConfig['conpair_directory'], tumorPileupFile, normalPileupFile,
                                                   concordanceOutFile)
            # run ConPair contamination tool
            if os.path.isfile(contaminationOutFile) and ConPairConfig['rewrite_previous_files'] == False:
                logger.warn("Skipping generating contamination for file " + contaminationOutFile + " because it already exists")
            else:
                run_contamination.estimate_tumor_normal_contamination(ConPairConfig['conpair_directory'],
                                                                      tumorPileupFile, normalPileupFile, contaminationOutFile)
            # update samplePairs table to mark this as a sample we can report out
            # but first check that the files were actually created
            if os.path.isfile(contaminationOutFile) and os.path.isfile(concordanceOutFile):
                samplePairs.loc[index, 'report_out'] = True
        else:
            # make sure neither sample failed pileup
            if normalPileupFile == "" or tumorPileupFile == "":
                logger.warn("Warning! Skipped concordance/contamination because one of the pileups for pair ",
                      + tumorSampleName + " could not be generated!")
            else:
                logger.warn("WARNING! Pileup file " + normalPileupFile + " and/or " + tumorPileupFile + " does not exist!!!")

    ###### WORKFLOW PART 3 ######
    ### gather all the ConPair output into an excel file for reporting
    logger.info("Collecting concordance/contamination results and putting it into an excel file")
    # create empty data frames for storing the values
    concordanceTable = pd.DataFrame({'tumor sample' : [],
                         'paired normal': [],
                         'concordance' : [],
                         'markers used' : [],
                         'min coverage for markers' : [],
                         'min map quality' : [],
                         'min base quality' : []})

    contaminationTable = pd.DataFrame({'tumor sample' : [],
                                       'paired normal' : [],
                                       'tumor sample contamination' : [],
                                       'normal sample contamination' : []})

    # import all the concordance/contaminaion outputs
    for index, samplePair in samplePairs.iterrows():
        # check that there wasn't an issue with generating
        # the concordance/contamination data
        if samplePair['report_out'] == True:
            contaminationFile = samplePair['contamination_output_file']
            concordanceFile = samplePair['concordance_output_file']
            normalSample = samplePair['normal_sample']
            tumorSample = samplePair['tumor_sample']
            # collect the concordance data
            # parsing concordance files is not generalizable, ConPair has its own
            # format, so we are just going to make a specific function to sanitize their output
            currConcordanceColumn = concordance_parser(pd.read_csv(concordanceFile, sep="\t", header=None),
                                                 tumorSample, normalSample)
            concordanceTable = concordanceTable.append(currConcordanceColumn)
            # contamination table has similar issues with parsing that the concordance
            # table does
            currContaminationColumn = contamination_parser(pd.read_csv(contaminationFile, sep="\t", header=None),
                                                 tumorSample, normalSample)
            contaminationTable = contaminationTable.append(currContaminationColumn)

    # merge the concordance/contamination tables on the tumor/normal columns
    # that way we can read the conrodance, marker sites, etc in the same table
    # as the contamination values
    ConPairReportDataFrame = concordanceTable.merge(contaminationTable, on=['tumor sample','paired normal'],
                                                    how='left')

    # create excel file and write the data to it
    conPairReport = xls.XLwrite(conpair_report_fileName,'ConPair_Report')
    conPairReport.writeTable(ConPairReportDataFrame,
                             styles={'normal sample contamination':'0.00%', 'tumor sample contamination':'0.00%'})
    conPairReport.save()

# EOF
