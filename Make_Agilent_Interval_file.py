#!/usr/bin/env python
# version 1.0
# tested with python 2.7.13

__author__ = "Bruce Wollison"
__date__ = "Jan, 2018"

'''
    STATUS: UNFINISHED
'''

'''
    Description: This script will generate interval files that are used in
    Riker Pipeline from Agilent Bait Set Design files, as well as the XML 
    configuration that needs to be inserted into the config
    Steps
      1) Made bed files from input folders using bedtools
      2) Convert these to picard interval files
      3) Create interval files from set of input folders using Picard Tools
      4) Create XML output file for entering into the config

    Notes:
        Adapted from AgilentBaitSetAnnotation written by M Ducar
        Original program was written in Java, has some 'gotchas', and doesn't automate
        everything, this script should replace it functionally and make it easier
        to run without background information
        
    EX Final Picard call:
    java -jar picard.jar IntervalListTools 
    INPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/exons/Gutierrez_POPv3_Plus_NoTRL_3042751_p0_BaitedRegions.interval_list 
    INPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/snps/Gutierrez_POPv3_Plus_NoTRL_3042751_p1_BaitedRegions.interval_list 
    INPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/snps/Gutierrez_POPv3_Plus_NoTRL_3042751_p2_BaitedRegions.interval_list 
    OUTPUT=/home/bmw35/Desktop/Gutierrez_POPv3_Plus_NoTRL/161228_Gutierrez_POPv3_Plus_NoTRL_3042751_BaitedRegions.interval_list 
    UNIQUE=true SORT=true COMMENT="Baited regions for bait set Gutierrez_POPv3_Plus_NoTRL_3042751"
    
    Step1 makes covered, notcovered and baitedRegion interval_list files
    Covered are $bedtools intersect -a [...Covered.bed] -b [...Regions.bed]
    NotCovered are $bedtools subtract -a [...Regions.bed] -b [...Covered.bed] 
    BaitedRegion are just the Covered.bed
    EX:bedtools intersect -a Weinstock_FL_v2.0-CARD10_exons_Covered.bed -b Weinstock_FL_v2.0-CARD10_exons_Regions.bed
    This is pretty straightforward with pybedtools, ex:
    Covered = BedTool('Weinstock_FL_v2.0-CARD10_exons_Covered.bed')
    Region = BedTool('Weinstock_FL_v2.0-CARD10_exons_Regions.bed')
    print(Covered.intersect(Region))
    It is really that simple.
    
    Step2 read in bed, convert its coordinates from base 0 to base 1, give interval classification
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
logger.info("Starting to make Agilent interval files with command: " + " ".join (sys.argv))

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

# lets make a few more variable names that will be useful later
# starting with output folders of the ConPair analysis
project_folder = args.output_folder + "/" + args.project_code
conpair_report_fileName = project_folder + "/" + args.project_code + "_report.xlsx"


##################################################
### define any small script-specific functions ###
##################################################

########################
### it's main() time ###
########################

if __name__ == "__main__":

# EOF
