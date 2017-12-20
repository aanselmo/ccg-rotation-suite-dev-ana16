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
        except ValueError as err:
            print("Something went wrong: {}".format(err))


    # general function for performing a query and returning a panda data frame
    # must be able to deal with any type of query return so that we have a common
    # way to communicate empty results or errors
    def DBquery(self, queryToRun):
        # debug line prinitng query:
        # print queryToRun

        self.cursor.execute(queryToRun)
        queryResults = list(self.cursor)
        queryDescriptors = self.cursor.description
        # DB column name is in first element of each description array
        queryHeaders = [descriptor[0] for descriptor in queryDescriptors]

        # debug line printing the return form the query:
        # print(queryHeaders)
        # print(queryResults)

        # convert results to pandas dataframe object to return
        resultsToDF = pd.DataFrame.from_records(data=queryResults, columns=queryHeaders)
        return resultsToDF


    # automatically close the connection when the program ends
    def __del__(self):
        self.cursor.close()
        self.connection.close()
        print("disconnected from " + self.host)

# EOF