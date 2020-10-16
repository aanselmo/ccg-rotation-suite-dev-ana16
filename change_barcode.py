#!/usr/bin/env python

__author__ = "Neil A. Patel"
__date__ = "October 16, 2020"

import sys, os, argparse, yaml, shutil, fnmatch
import pandas as pd

import utils.sql.dbConnect as db
import datetime
# default config file needed for connecting to the DB info
DEFAULT_CONFIGURATION_FILE = 'ccgd_rotation_scripts.config'


def setup_arguments():
  """
  This is the function that sets up the flags and the arguements you can pass to the script.
  :author: Bruce Wollison, Neil A. Patel
  """
  parser = argparse.ArgumentParser(description=('This script generates reports for BreaKmer based on '
                                   + 'the given project code.'),
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--ignore', '-i', help=('Creates report regardless of what errors would otherwise'
                      + 'cause the script to die.'),
                      type=bool, default=False)
  parser.add_argument('--config', '-c', help='configuration file to use',
                      default=DEFAULT_CONFIGURATION_FILE, type=str)
  args = parser.parse_args()
  return args


def connect_to_riker_db(args_config):
    ### import configurations from configuration file ###
    print("Reading in configuration file " + args_config)
    try:
        with open(args_config, 'r') as ymlfile:
            config = yaml.load(ymlfile)
    except IOError:
        # whoa whats happening bros? File might not exist, or some other IO error
        print("Error reading configuration file!")
        sys.exit()
    DBconfig = config['rikerDB']
    return DBconfig


def get_barcode_set_member_name(DBconn, first_barcode, second_barcode):
    query_statement = ("SELECT * FROM seq_barcode_set_member sbsm "
        + "WHERE first_seq_barcode_id  = "
        + "(SELECT seq_barcode_id FROM seq_barcode WHERE seq_barcode.barcode_sequence = %(first_barcode)s) "
        + "AND second_seq_barcode_id = "
        + "(SELECT seq_barcode_id FROM seq_barcode WHERE seq_barcode.barcode_sequence = %(second_barcode)s);")

    query_results = DBconn.DBquery(query_statement, {'first_barcode' : first_barcode, 
        'second_barcode': second_barcode})
    return query_results


def get_good_name_and_id(member_name_options, good_name_contains):
   which_good = member_name_options.loc[:, 'seq_barcode_set_member_name'].str.contains(good_name_contains)
   try:
     good_name = member_name_options.loc[which_good, 'seq_barcode_set_member_name'].values[0]
     good_id = member_name_options.loc[which_good, 'seq_barcode_set_member_id'].values[0]
   except IndexError:
     good_name = 'NA'
     good_id = 'NA'
   return(good_name, good_id)


def get_update_query(member_id, library_name):
    return ("UPDATE library SET library.barcode_set_member_id = "
        + str(member_id)
        + " WHERE library_name = '" + library_name + "';")


if __name__ == "__main__":
    print('BreaKmer_report.py started')
    args = setup_arguments()
    DBconfig = connect_to_riker_db(args.config)
    dbConn = db.DBconnect(DBconfig)
    
    riker_export = pd.read_csv('riker_export_1105.csv')
    consise_export = riker_export[['Library Prep ID', 'External ID', 'Index ID', 'First Sequence', 'Second Sequence',
                                  'Method']]
    consise_export['good_name'] = ''
    consise_export['good_id'] = 0
    consise_export['update_query'] = '' 

    for index, row in consise_export.iterrows():
        member_names = get_barcode_set_member_name(dbConn, row['First Sequence'], row['Second Sequence'])
        print(row['Library Prep ID'])
        (good_name, good_id) = get_good_name_and_id(member_names, 'IDT_384_UDI_UMI')
        consise_export.loc[index, 'good_name'] = good_name
        consise_export.loc[index, 'good_id'] = good_id
        consise_export.loc[index, 'update_query'] = get_update_query(good_id, row['Library Prep ID'])
        
    
    consise_export.to_csv('riker_export_1105_change_to.csv')
# EOF
