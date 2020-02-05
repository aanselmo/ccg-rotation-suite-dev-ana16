#!/usr/bin/env python

"""
  WARNING: This script relies on an old script utils/BreaKmer/generate_report.R which
  requires an R instance with the libraries RMySQL and XLConnect. Additionally this
  script has DB connetions with credentials stored in plain-text. Distribute with caution.

  Description: This script will generate a breakmer report for a given project

  Workflow is broken into 3 steps:
  1) handle inputs, query the database, get all the information we need to verify that the
    run is finished and we have what we need
  2) create symlinks from the output of the breakmer run so we can run the process_cohort.py script
    and run it to generate the report data files
  3) generate the excel report based on the data files.

  NOTE: this script is a stop-gap measure meant to avoid 'manually' generating a script through multiple
  time consuming steps. It is not meant to be fully polished production ready and maybe has a few shortcuts here and
  there. Off the top of my head, I think there is no check to see if there are any indel/rearrangements, and the
  script may not successfully generate and empty file if there are no events found.
"""

__author__ = "Bruce Wollison, Neil A. Patel"
__credits__ = ["Bruce Wollison", "Neil A. Patel"]

import sys
import os
import argparse
import yaml
import shutil
import pandas as pd
import datetime

import utils.sql.dbConnect as db
import utils.general.general_functions as genFun
import utils.sql.common_queries as com_queries
import utils.BreaKmer.process_cohort_results as process_cohort
# using openpyxl directly because we are editing, not creating, reports:
import openpyxl as xl
from openpyxl.styles import Font, PatternFill

# 'id' in the riker DB that points to the current breakmer run tool
# likely won't change often or ever
DEFAULT_BREAKMER_TOOL_ID = '13'
# default config file needed for connecting to the DB info
DEFAULT_CONFIGURATION_FILE = 'ccgd_rotation_scripts.config'
DEFAULT_OUTPUT_DIRECTORY = './BreaKmer_Reports'
REARRANGE_FILE_PATTERN = "_Breakmer_0_0_6_rearrangement_svs.out"
INDEL_FILE_PATTERN = "_Breakmer_0_0_6_indel_svs.out"
YES = {'yes', 'y', 'ye', ''}
NO = {'no', 'n'}


def setup_arguments():
  """
  This is the function that sets up the flags and the arguements you can pass to the script.
  :author: Bruce Wollison, Neil A. Patel
  """
  parser = argparse.ArgumentParser(description=('This script generates reports for BreaKmer based on '
                                   + 'the given project code.'),
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--project_code', '-p', help='Project code to generate the report off of',
                      type=str, required=True)
  parser.add_argument('--ignore', '-i', help=('Creates report regardless of what errors would otherwise'
                      + 'cause the script to die.'),
                      type=bool, default=False)
  parser.add_argument('--config', '-c', help='configuration file to use',
                      default=DEFAULT_CONFIGURATION_FILE, type=str)
  args = parser.parse_args()
  return args


def breaKmerDirectoriesFromProjectID(DBconn, project_id, breakmer_tool_id):
  """
  get a panda dataframe using the project ID
  returns all complete aggregate pairs for a project and their associated completed breakmer runs
  should look at the whole panda df. If there are nulls for the breakmer path that means it was
  not run or completed on those aggregate pairs, worth mentioning
  """
  querySTMT =  "SELECT ap.aggregate_pair_id, ap.aggregate_pair_name, atr.* FROM aggregate_pair ap " +\
         "left join analysis_tool_run atr on atr.aggregate_pair_id = ap.aggregate_pair_id " + \
         "and atr.annotation_status_id = 1 and " +\
         "atr.analysis_tool_id = %(breakmer_id)s and atr.status_id = 1 " + \
         "where ap.qc_status_id = 1 and ap.aggregate_pair_project_id = %(proj_id)s"
  queryResults = DBconn.DBquery(querySTMT, {'proj_id' : project_id, 'breakmer_id': breakmer_tool_id})
  return queryResults


def get_project_information(projectCode, dbConn):
  """
  Gets basic project information from the project_code.
  :return: a tuple of project_id, project_name, and number of samples in the project.
  """
  projectIDQueryResults = com_queries.get_project_id_from_code(dbConn, projectCode)
  if projectIDQueryResults.empty:
    raise("Error, project code " + projectCode + " does not exist in " + DBconfig['schema'] + "!!! Exiting...")
  else:
    project_id = projectIDQueryResults.project_id[0]
  projectDescriptionQueryResults = com_queries.get_project_name_from_code(dbConn, projectCode)
  project_name = projectDescriptionQueryResults['project_name'][0]
  numberOfSamplesQueryResults = com_queries.count_number_of_samples_in_project(dbConn, projectCode)
  numberOfSamples = numberOfSamplesQueryResults['sampleCount'][0]
  return (project_id, project_name, numberOfSamples)


def prompt_and_filter_unfinished_agg_pairs(breaKmerDirectoriesQueryResults): 
  unFinishedAgPairs = breaKmerDirectoriesQueryResults[
      (breaKmerDirectoriesQueryResults['status_id'] != 1) | 
      (breaKmerDirectoriesQueryResults['annotation_status_id'] != 1)]
  if unFinishedAgPairs.empty:
    print("All aggregate pairs have finished processing.")
  else:
    print("Some aggregates have not finished processing. Specifically ag pairs: "
          + (", ").join(unFinishedAgPairs['aggregate_pair_name'])
          + "\nContinue generating report with currently finished ags? Y/N")
    # raw_input returns the empty string for "enter"
    choice = input().lower()
    if choice in NO:
      # note: This is a little convoluted, the query returns two aggregate_pair_id columns, the first is the
      # one we want because the second is expected to be empty in this case (since it represents ag pairs without
      # a processed breakmer analysis tool run). So we need .iloc[:,0] in pandas to access the first column with
      # ag pair name, then need to convert all the values to strings from numbers
      agPairIDString = (",").join([str(agPairID) for agPairID in unFinishedAgPairs.aggregate_pair_id.iloc[:,0]])
      print("nohup ./AggregatePairAnalysisforAnalysisTools.pm -AggregatePairId "
          + agPairIDString + " -AnalysisTool Breakmer_0.0.6 "
          + " -config ./Riker.ccgd.argos.config &> log/AnalysisTools_"
          + projectCode + "_breakmer.txt &")
      sys.exit()
    else:
      # we need to remove the unfinished aggregate pairs from any downstream processing
      # we will remove these from the breaKmerDirectoriesQueryResults
      # TODO: this is a copy+paste solution. Should take the inds of unFinishedAgPairs and reuse
      # the list of inds here
      breaKmerDirectoriesQueryResults = breaKmerDirectoriesQueryResults[
          (breaKmerDirectoriesQueryResults['status_id'] == 1) |
          (breaKmerDirectoriesQueryResults['annotation_status_id'] == 1)]
  return breaKmerDirectoriesQueryResults


def prompt_and_filter_missing_breakmer_dirs(breaKmerDirectoriesQueryResults):
  # check that all the files are available on the file system
  noLongerAvailableBreaKmerRuns = breaKmerDirectoriesQueryResults[ \
    (breaKmerDirectoriesQueryResults['data_location_id'] != 1) & \
    (breaKmerDirectoriesQueryResults['status_id'] == 1) & \
    (breaKmerDirectoriesQueryResults['annotation_status_id'] == 1)]

  if not noLongerAvailableBreaKmerRuns.empty:
    print("Some breaKmer outputs have already been backed up. Specifically ag pairs: " + \
        (", ").join(noLongerAvailableBreaKmerRuns.aggregate_pair_name) + \
        "\nContinue generating report with currently finished ags? Y/N")
    # raw_input returns the empty string for "enter"
    choice = input().lower()
    if choice in NO:
      sys.exit()
    else:
      # we need to remove the not available aggregate pairs from any downstream processing
      # we will remove these from the breaKmerDirectoriesQueryResults
      # TODO: this is a copy+paste solution. Should take the inds of noLongerAvailableBreakmerRuns and reuse
      # the list of inds here
      breaKmerDirectoriesQueryResults = breaKmerDirectoriesQueryResults[ \
        (breaKmerDirectoriesQueryResults['data_location_id'] == 1) & \
        (breaKmerDirectoriesQueryResults['status_id'] == 1) & \
        (breaKmerDirectoriesQueryResults['annotation_status_id'] == 1)]
  return breaKmerDirectoriesQueryResults


def prompt_and_create_output_dir(projectCode):
  """
  Creates the directory where links will be made to breakmer files. Prompts the user whether to delete the directory
  if it already exists.
  """
  projectOutputDirectory = DEFAULT_OUTPUT_DIRECTORY + "/" + projectCode + "_breaKmerData"
  if not os.path.exists(projectOutputDirectory):
    os.mkdir(projectOutputDirectory)
  else:
    # warn that the directory already exists and that we would want to
    # overwrite everything.
    print("This script appears to already have been run on this project. " + \
        "\nDelete current BreaKmer results directory to Continue? Y/N")
    # raw_input returns the empty string for "enter"
    choice = input().lower()
    if choice in YES:
      shutil.rmtree(projectOutputDirectory)
      os.mkdir(projectOutputDirectory)
    else:
      sys.exit()
  return projectOutputDirectory


def update_mount_point_prefix(breaKmerDirectoriesQueryResults, mount_point):
  """
  The database has the location of output files for breakmer. However, depending on where you mount the storage unit
  for these files on a different server or your desktop, the prefix/parent directories will change for the path.
  This function updates the paths for these breakmer output folders.
  
  important to abstract out of this step anything specific to the query
  in the future we will likely want multiple paths through workflow pt1, one where
  we don't hit Riker but go off user input of starting directories when we
  are running something outside of the DB
  additionally we are also replacing the DB name with any local paths to the isilon
  that we define in the config.
  """
  breakmerParentDirectories = breaKmerDirectoriesQueryResults[['aggregate_pair_name','output_file']]
  newBreaKmerParentDir = []
  for parentDir in breakmerParentDirectories['output_file']:
    newDir = genFun.replaceWCommonSubStr(mount_point, parentDir)
    newBreaKmerParentDir.append(newDir)
  breakmerParentDirectories['output_file'] = newBreaKmerParentDir
  return breakmerParentDirectories


def create_breakmer_file_links(projectOutputDirectory, breakmerParentDirectories):
  """
  need to make links for the following files: "*rearrangement_svs.out", "*indel_svs.out", "output/"
  that way Ryans other script can find them on the file system and extract the needed info
  """
  for aggPairName,directoryToLink in zip(breakmerParentDirectories['aggregate_pair_name'], \
                                   breakmerParentDirectories['output_file']):
    print("creating symlinks for aggregate pair" + aggPairName +\
      " using directory " + directoryToLink)
    # for each parent dir path, create one of each symlink needed (indel, rearrangement, outputfolder)
    # first collect all the file paths we need fro the breakmer output
    breaKmerOutputDirectoryPath = directoryToLink + "/output"
    # NOTE: output directory must exist. Indel/rearrangement sys.out may not!
    if os.path.exists(breaKmerOutputDirectoryPath):
      allBreakmerFiles = os.listdir(breaKmerOutputDirectoryPath)
    else:
      print("Warning!!!! Did not find the following expected output directory:\n" + \
          breaKmerOutputDirectoryPath + "\nCheck this before continuing")
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


def create_sample_list_files(projectOutputDirectory, projectCode, breakmerParentDirectories, projectDescription):
  """
  Create the filters file, normal file, and tumors file. These are neede for running the
  process_cohort_results.py script. After we do this, lets run that script to generate the desired
  output. The normal file contains the list of normal samples (samples that we can use to filter out variants
  and that we do not want the variants listed out for. Historically, we have left this blank so the collaborator
  can see the variants in the normals and do their own filtering.
  """
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
    for aggPair in breakmerParentDirectories['aggregate_pair_name']:
      tumorSampleListFile.write(("\t").join([projectDescription, aggPair]))
      tumorSampleListFile.write("\n")


def add_event_counts_to_summary(resultsSummaryFilePath, summaryWS, numberOfSamples):
  """
  now we just need to add the summary information to a summary tab
  start by reading in the relevant info.
  """
  inversion_count = indel_count =  rearrangement_count =  sample_count = \
    sample_event_count = '0'
  # read the file and check each line which aspect is corresponds to
  with open(resultsSummaryFilePath) as summaryFile:
    for line in summaryFile:
      if 'samples with event' in line :
        sample_event_count = line.split('\t')[1]
      elif 'samples' in line:
        sample_count = line.split('\t')[1]
      elif 'indel' in line:
        indel_count = line.split('\t')[1]
      elif 'rearrangement events' in line:
        rearrangement_count = line.split('\t')[1]
      elif 'rearrangement' in line:
        inversion_count = line.split('\t')[1]
      else:
        print('unrecognized line in the summary file!\n')
  summaryWS['B2'] = str(sample_count).rstrip("\n")
  summaryWS['B3'] = str(sample_event_count).rstrip("\n")
  summaryWS['B4'] = numberOfSamples
  summaryWS['B5'] = str(rearrangement_count).rstrip("\n")
  summaryWS['B6'] = str(inversion_count).rstrip("\n")
  summaryWS['B7'] = str(indel_count).rstrip("\n")
  return summaryWS


def remove_bkpt_bed_files():
  """
  clean up the bkpt_*.bed files that appear. Going to remove them individually
  because if we remove things programatically or using patterns I am worried we'll
  delete something we want for whatever reason
  """
  brkptRemoveFiles = ['brkpt_cluster_anno_intersect.bed', 'brkpt_cytoband_anno.bed',
                'brkpt_exon_anno_intersect.bed', 'brkpt_gene_anno_downstream.bed',
                'brkpt_gene_anno_intersect.bed', 'brkpt_gene_anno_upstream.bed',
                'brkpts.normal.bed', 'brkpts.tumor.bed']
  [os.remove(removeFile) for removeFile in brkptRemoveFiles]


def main():
  args = setup_arguments()
  print("Reading in configuration file " + args.config)
  try:
    with open(args.config, 'r') as ymlfile:
      config = yaml.load(ymlfile)
  except IOError:
    print("Error reading configuration file!")
    sys.exit()
  DBconfig = config['rikerDB']
  genConfig = config['general']
  bedToolsPath = genConfig['bedtools_path']
  RPath = genConfig['R_path']
  mount_point = genConfig['isilon_project_directory']
  projectCode = args.project_code
  
  print('BreaKmer_report.py started')
  dbConn = db.DBconnect(DBconfig)
  project_id, Description, numberOfSamples = get_project_information(projectCode, dbConn)

  projectDescription = "_".join([projectCode, Description])
  revProjectDescription = "_".join([Description, projectCode])
  
  breaKmerDirectoriesQueryResults = breaKmerDirectoriesFromProjectID(dbConn, project_id, DEFAULT_BREAKMER_TOOL_ID)
  if breaKmerDirectoriesQueryResults.empty:
    print("Error, no aggregate pairs exist for project  " + projectCode + " in " +
          DBconfig['schema'] + "!!! Exiting...")
    sys.exit()

  breaKmerDirectoriesQueryResults = prompt_and_filter_unfinished_agg_pairs(breaKmerDirectoriesQueryResults)
  breaKmerDirectoriesQueryResults = prompt_and_filter_missing_breakmer_dirs(breaKmerDirectoriesQueryResults)
  breakmerParentDirectories = update_mount_point_prefix(breaKmerDirectoriesQueryResults, mount_point)
  projectOutputDirectory = prompt_and_create_output_dir(projectCode)
  create_breakmer_file_links(projectOutputDirectory, breakmerParentDirectories)
  create_sample_list_files(projectOutputDirectory, projectCode, breakmerParentDirectories, projectDescription) 

  # This will go through the breakmer data on the file system and generate the necessary data for reporting
  process_cohort.processCohortResults(projectDescription, projectCode, projectOutputDirectory, bedToolsPath)
  remove_bkpt_bed_files()

  # move the summary file to the appropriate project specific folder
  resultsSummaryFilePath = ''.join([projectOutputDirectory, "/", 'BreaKmer_Summary.txt'])
  os.rename('BreaKmer_Summary.txt', resultsSummaryFilePath)

  Rcmnd = RPath + " utils/BreaKmer/./generate_report.R " + projectCode + ' ' + projectOutputDirectory + '/'
  print('running R command: ' + Rcmnd)
  os.system(Rcmnd)

  now = datetime.datetime.now()
  fullReport = (projectOutputDirectory + '/' + now.strftime("%y%m%d") + '_' + revProjectDescription
      + '_Breakmer_report.xlsx')
  fullReportXL = xl.load_workbook(fullReport)
  summaryWS = fullReportXL.worksheets[0]
  summaryWS = add_event_counts_to_summary(resultsSummaryFilePath, summaryWS, numberOfSamples)

  # bold the first row of all sheets except for the circos and summary ones
  for currSheet, sheetName in zip(fullReportXL, fullReportXL.sheetnames):
    if 'circos' not in sheetName and 'Summary' not in sheetName:
      for currCell in currSheet["1:1"]:
        currCell.font = Font(name='Calibri', bold=True)
        currCell.fill = PatternFill(fill_type=None)

  # save our final workbook now that it is complete
  fullReportXL.save(fullReport)


if __name__ == "__main__":
  main()
