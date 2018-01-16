#!/usr/bin/env python
# version 1.0

__author__ = "Bruce Wollison"
__date__ = "Jan, 2018"

'''
    Description: This script will generate interval files that are used in
    Riker Pipeline from Agilent Bait Set Design files, as well as the XML 
    configuration that needs to be inserted into the config
    Steps
      1) Made bed files from input folders using bedtools
      2) Create interval files from set of input folders using Picard Tools
      3) Create XML output file for entering into the config

    Notes:
        Adapted from AgilentBaitSetAnnotation written by M Ducar
        Original program was written in Java, has some 'gotchas', and doesn't automate
        everything, this script should replace it functionally and make it easier
        to run without background information
        
    EX Picard call:
    java -jar picard.jar IntervalListTools 
    INPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/exons/Gutierrez_POPv3_Plus_NoTRL_3042751_p0_BaitedRegions.interval_list 
    INPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/snps/Gutierrez_POPv3_Plus_NoTRL_3042751_p1_BaitedRegions.interval_list 
    INPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/snps/Gutierrez_POPv3_Plus_NoTRL_3042751_p2_BaitedRegions.interval_list 
    OUTPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/161228_Gutierrez_POPv3_Plus_NoTRL_3042751_BaitedRegions.interval_list 
    UNIQUE=true SORT=true COMMENT="Baited regions for bait set Gutierrez_POPv3_Plus_NoTRL_3042751"
'''

#######################
### library imports ###
#######################

import pandas as pd
import sys, imp, os, argparse, yaml
from pybedtools import BedTool

# add special modules to the path
import general_functions as genFun


##########################################
### globals, parsing and configuration ###
##########################################
DEFAULT_LOGGER_FILE = "AgilentIntervalFile.log"
DEFAULT_CONFIGURATION_FILE = 'ccgd_rotation_scripts.config'
DEFAULT_OUTPUT_DIRECTORY = os.getcwd() + '/Agilent_Interval_Files'

# set up logger
logger = genFun.getLogger(DEFAULT_LOGGER_FILE)
logger.info("Starting to make Agilent itnerval files with command: " + " ".jol.com / in (sys.argv))

### handle command line interface with argparse ###
parser = argparse.ArgumentParser(description='Make Riker pipeline ready interval files from Agilent Bait Design files.',
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
logger.info("Reading in configuration file " + args.config)
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
logger.info("configuration file read in successfully")

# add the conpair module directory to the pythonpath variable
# then import the ConPair scripts
ConPairConfig = config['conpair']

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

########################
### it's main() time ###
########################

if __name__ == "__main__":

# EOF
