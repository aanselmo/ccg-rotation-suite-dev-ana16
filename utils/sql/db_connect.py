#!/usr/bin/env python
# version 1.0

__author__ = "Bruce Wollison"
__date__ = "Nov, 2017"

'''
    Description: This script contains functions that connect to the riker mysql DB
    and run necessary queries. It should provide functions to get
    the needed info from the db so we can rewrite this module if needed
    for a different db/implementation
    
    Ideally we instantiate with the DB credentials and connect, leave the close 
    connection to the destructor, and then have a query wrapper function
    that calls queries based on a DB model and returns a panda data frame
    of the values or a specific value (meaning we want to make sure we handle 
    queries with 0 returned contents in a predictable way, while also being able to 
    pass on errors in a predictable way).
'''

# library imports
import MySQLdb as mysql
import pandas as pd
import numpy as np

### wrapping functions from the DB API here so we can keep them abstracted
### as opposed to having an implementation that is specific to a DB
### this way we can reimplement with a diff API or in a diff DB with
### hopefully minimal, if any, changing of the other scripts

# creating a new class for the DB connection that takes in the
# connection credentials from the config file as a Dict
class DBconnect:

    def __init__(self, DBConfigDict):
        """
        bmw35 12-29-2017
        Constructor
        :param DBConfigDict: takes a config dictionary that contains the needed 
          elements to start a connection
        :return: returns a DBconnect object, which serves as the database connection for the programs
           duration
        """
        # test to make sure dict items exist properly
        self.schema = DBConfigDict['schema']
        self.host = DBConfigDict['host']
        self.port = DBConfigDict['port']
        self.user = DBConfigDict['user']
        self.passwd = DBConfigDict['password']
        self.connection = mysql.connect(db=DBConfigDict['schema'], host=DBConfigDict['host'],
                            port=DBConfigDict['port'], user=DBConfigDict['user'], passwd=DBConfigDict['password'])
        # initiate connection to the db using the given credentials
        try:
            self.cursor = self.connection.cursor()
            print("successfully connected to database " + self.host)
        # Oops. Didn't work! Can't return a db object
        except ValueError as err:
            print("Something went wrong: {}".format(err))
            return None


    def DBquery(self, queryToRun):
        """
        bmw35 12-29-2017
        generalized function for performing a query and returning a panda data frame
        must be able to deal with any type of query return so that we have a common
        way to communicate empty results or errors
        :param queryToRun: string of the query that we want to hit the current DB connection
        :return: pandas Data Frame, returns the query results as a pandas dataframe
        """
        try:
            self.cursor.execute(queryToRun)
        # Issue with the query
        except mysql.DatabaseError as err:
            print("Something went wrong: {}".format(err))
            # can't return empty list because that could be interperted
            # as the query returning an empty table, return None type
            # and let whatever is getting the query deal with that
            return None
        queryResults = list(self.cursor)
        queryDescriptors = self.cursor.description
        # DB column name is in first element of each description array
        queryHeaders = [descriptor[0] for descriptor in queryDescriptors]
        # convert results to pandas dataframe object to return
        resultsToDF = pd.DataFrame.from_records(data=queryResults, columns=queryHeaders)
        return resultsToDF


    def __del__(self):
        """
        bmw35 12-29-2017
        Destructor. automatically close the connection when the program ends
        """
        self.cursor.close()
        self.connection.close()
        print("disconnected from " + self.host)

# EOF