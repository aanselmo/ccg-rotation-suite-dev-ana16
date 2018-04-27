#!/usr/bin/env python
# version 1.0
# tested with python 2.7.13

__author__ = "Bruce Wollison"
__date__ = "Jan, 2018"

'''
    STATUS: UNFINISHED
    TODO: Include making a gc content file for the interval file
'''

'''
    Description: This script will generate interval files that are used in
    Riker Pipeline from Agilent Bait Set Design files, as well as the XML 
    configuration that needs to be inserted into the config
    Steps
      1) Made bed files from input folders using bedtools
         a) This may be different for each vendor, so abstract it out making the bed files
            so Agilent specific aspects can be redone since NEB and IDT may be different. Subsequent
            steps should be the same no matter the vendor assuming we are starting from bed files,
            which would be sensible otherwise we may just use the interval list
      2) Convert these to picard interval files
         a) Assume we are in bed format, so change from base 0 to base 1
         b) create and organize Baited/Covered/NotCovered interval files for each feature type
      3) Create interval files from set of input folders using Picard Tools
         a) consists of taking the interval files for each feature (exon/intron/trl) and turning 
            them into a Baited, Covered, and NotCovered interval file for the whole bait set
      4) Create XML output file for entering into the config

    Notes:
        Adapted from AgilentBaitSetAnnotation written by M Ducar
        Original program was written in Java, has some 'gotchas', and doesn't automate
        everything, this script should replace it functionally and make it easier
        to run without background information
            
    ### STEP1 ###
    Assumptions (based on Agilent design files): 
      1) Each subdirectory has 1 set of design files
      2) Each subdirectory is either a single genomic type (ex: exon/intron/trl) or mix
      
    Step1 makes covered, notcovered and baitedRegion interval_list files
    Covered are $bedtools intersect -a [...Covered.bed] -b [...Regions.bed]
    NotCovered are $bedtools subtract -a [...Regions.bed] -b [...Covered.bed] 
    BaitedRegion are just the Covered.bed
    EX:bedtools intersect -a Weinstock_FL_v2.0-CARD10_exons_Covered.bed -b Weinstock_FL_v2.0-CARD10_exons_Regions.bed
    This is pretty straightforwar you can changed with pybedtools, ex:
    Covered = BedTool('Weinstock_FL_v2.0-CARD10_exons_Covered.bed')
    Region = BedTool('Weinstock_FL_v2.0-CARD10_exons_Regions.bed')
    print(Covered.intersect(Region))
    It is really that simple.
    
    Question: Do we want to have the uer specifcy as options which feature type each folder are, or ask them in a prompt?
    
    ### STEP2 ###
    Step2 read in bed, convert its coordinates from base 0 to base 1, give interval classification
    Parallel: Might want to have a decision tree here for BreaKmer bed files: if the baitset is under some
    amount of real estate we make all of them the bed file, otherwise maybe we only use some feature types
    
    ### STEP3 ###
    EX Final Picard call:
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
logger.info("Starting to make Agilent interval files with command: " + " ".join (sys.argv))

### handle command line interface with argparse ###
parser = argparse.ArgumentParser(description='This script Creates Pipeline ready interval files from '\
                                             'bed Design files.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vendor', '-v', help='bait set vendor (one of: Agilent)',
                    type=str, default='Agilent')
parser.add_argument('--baitset_name', '-n', help='name for the baitset',
                    type=str, required=True)
parser.add_argument('--input_directory', '-d', help='input directory for bait set bed files',
                    type=str, required=True)
parser.add_argument('--output_folder', '-o', help='folder to output the results to',
                    default=DEFAULT_OUTPUT_DIRECTORY, type=str)
parser.add_argument('--config', '-c', help='configuration file to use',
                    default=DEFAULT_CONFIGURATION_FILE, type=str)
parser.add_argument('--version', action='version', version='Make_Agilent_Interval_file.py ver 1.0')
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
logger.info("configuration file read in successfully")

# lets make a few more variable names that will be useful later
baitset_output_folder = args.output_folder + "/" + args.baitset_name
designsToFeatures = {}

##################################################
### define any small script-specific functions ###
##################################################

########################
### it's main() time ###
########################

if __name__ == "__main__":

    # check to see if the bait set directory already exists in the output
    # directory, warn if it exists
    if not os.path.exists(baitset_output_folder):
        os.mkdir(baitset_output_folder)
    else:
        logger.warn("Bait set output folder " + baitset_output_folder + " already exists")

    # for feature type create an output folder, but check and warn if it already exists
    for featureType in designsToFeatures:
        currFeatureOutputDir = baitset_output_folder + '/' + featureType + '_intervalFiles'
        if not os.path.exists(currFeatureOutputDir):
            os.mkdir(currFeatureOutputDir)
        else:
            logger.warn("Bait set feature output folder " + currFeatureOutputDir + " already exists")

    #### Agilent specific stuff here, abstract this to a specific
    # need to search for all subfolders in the input directory
    # also need to warn about any directory without design files


# EOF
