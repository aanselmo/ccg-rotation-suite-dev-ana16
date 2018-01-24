#!/usr/bin/env python
# version 1.0

__author__ = "Bruce Wollison"
__date__ = "Dec, 2017"

'''
    Description: This script contains common queries that we will want
    to reuse across multiple scripts. All functions are the queries 
    needed to perform the stated functions. All queries take in the DB connection
    object, and the values we want to query on, then execute the query command
    using the mysql cursor executor to avoid any possibility of injections.
    params are always:
    :param DBconn - a database connection object that executes queries with the method DBquery()
    :param someOption - anything used in the query to specify the results, like the project code string
'''

# get the project id from the project code
def projectIDFromProjectCode(DBconn, project_code):
    querySTMT =  "SELECT project_id FROM project WHERE project_code = %(proj_code)s"
    queryResults = DBconn.DBquery(querySTMT, {'proj_code' : project_code})
    return queryResults

# return samples that have a normal/tumor pair
# for a specific project using the project_id
def projectSamplePairsFromProjectID(DBconn, project_id):
    querySTMT = "SELECT rec1.sample_collaborator_name as normal_sample, "\
                 "rec1.received_sample_id as normal_sample_id, "\
                 "rec2.sample_collaborator_name as tumor_sample, rec2.received_sample_id as tumor_sample_id "\
                 "FROM received_sample rec1 "\
                 "join received_sample rec2 on rec1.individual_id = rec2.individual_id "\
                 "where rec1.project_id = %(proj_id)s "\
                 "and rec2.project_id = %(proj_id)s "\
                 "and rec1.genomic_source_id = 2 and rec2.genomic_source_id = 1;"
    queryResults = DBconn.DBquery(querySTMT, {'proj_id' : project_id})
    # this query just returns all samples for the project:
    # querySTMT =  ("SELECT * FROM received_sample WHERE project_id = '" + str(project_id) + "'")
    return queryResults

# get the aggregate pair information for a given project_id
def aggregatePairTableFromProjectID(DBconn, project_id):
    querySTMT = "SELECT ap.*, ag1.received_sample_id as case_sample_id, ag1.aggregate_bam_file as case_sample_path, "\
                 "ag2.received_sample_id as control_sample_id, ag2.aggregate_bam_file as control_sample_path "\
                 "FROM aggregate_pair ap "\
                 "join aggregation_pipeline_run apr "\
                 "on ap.aggregation_pipeline_run_id = apr.aggregation_pipeline_run_id "\
                 "left join aggregate ag1 on ap.case_aggregate_id = ag1.aggregate_id "\
                 "left join aggregate ag2 on ap.control_aggregate_id = ag2.aggregate_id "\
                 "where aggregate_pair_project_id = %(proj_id)s "\
                 "and ap.qc_status_id = 1 "\
                 "and apr.aggregation_pipeline_run_status_id = 1;"
    queryResults = DBconn.DBquery(querySTMT, {'proj_id' : project_id})
    return queryResults

# EOF