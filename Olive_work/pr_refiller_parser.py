# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 11:48:39 2024

@author: rjovelin
"""



import sqlite3
import json
import requests
import gzip
import argparse
import time
import traceback
import os
import subprocess
import hashlib
# from utilities import connect_to_db
# from whole_genome import get_workflow_limskeys, find_WGS_blocks
# from whole_transcriptome import find_WT_blocks



# data -> list of donor data
# dict_keys(['cerberus_data',
#            'pinery_data',
#            'pinery_project_data',
#            'nabu_data',
#            'donor'])


# data[0]['donor']



# data[0]['nabu_data'] ->

# ========

# [[{'comment': 'GDR-839',
#    'file_id': 'vidarr:research/file/2f585f57205dbb2c52d96bc785b47b8ee3bd84ba931a5f8f2d4560b1520016c0',
#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25619611/AHCT_0029_01_LB01-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-TACTGCAG_R2.fastq.gz',
#    'file_swid': '25620329',
#    'fileqc_id': '58037',
#    'project': 'AHCT',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',
#    'username': 'rjovelin',
#    'workflow': 'bcl2fastq'}],
#  [{'comment': 'GDR-839',
#    'file_id': 'vidarr:research/file/31c12c06d91562e28142e4e95a394b6d9fc769ced2a2216d5e904096bb2aff04',
#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25618917/AHCT_0029_01_LB02-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-CGATTCCG_R2.fastq.gz',
#    'file_swid': '25619744',
#    'fileqc_id': '57741',
#    'project': 'AHCT',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',
#    'username': 'rjovelin',
#    'workflow': 'bcl2fastq'}],
#  [{'comment': 'GDR-839',
#    'file_id': 'vidarr:research/file/87056488f0a3a34da37114965e438008bacf7a1e2ca44a741c40eeecc1d81df0',
#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25619611/AHCT_0029_01_LB01-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-TACTGCAG_R1.fastq.gz',
#    'file_swid': '25620299',
#    'fileqc_id': '57996',
#    'project': 'AHCT',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',
#    'username': 'rjovelin',
#    'workflow': 'bcl2fastq'}],
#  [{'comment': 'GDR-839',
#    'file_id': 'vidarr:research/file/d586a196c2ed7bb95daaeb5eda6f01cdc91dc710d42b388eb24a2ede14257c02',
#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25618917/AHCT_0029_01_LB02-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-CGATTCCG_R1.fastq.gz',
#    'file_swid': '25619727',
#    'fileqc_id': '57569',
#    'project': 'AHCT',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',
#    'username': 'rjovelin',
#    'workflow': 'bcl2fastq'}]]
# =========

# i['cerberus_data'] -> list of lists

# [{'accession': 'vidarr:research/file/012cf6554e9b51f0dcab3186a134f02b26b78953534e56dbd14ba3158428c380',
#   'barcode': 'TTCTGGAA-CGATTCCG',
#   'external_donor_id': 'BMT-PMH-125',
#   'group_id': 'BMT-PMH-125',
#   'input_files': '[]',
#   'instrument_model': 'Illumina MiSeq',
#   'lane': '1',
#   'library_design': 'TS',
#   'library_name': 'AHCT_0029_01_LB02-01',
#   'library_type': '',
#   'lims': '{"id":"5893_1_LDI84239","provider":"pinery-miso","time":1676401595000,"version":"882752c33378b9d66793a619fc4918241285391bf490f0525325f90fd58db05b"}',
#   'md5': 'c55a9d4efd25ff5d081b4c3301453561',
#   'path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25547692/AHCT_0029_01_LB02-01_220418_M06816_0282_000000000-DG2YB_1_TTCTGGAA-CGATTCCG_R2.fastq.gz',
#   'project': 'AHCT',
#   'run': '220418_M06816_0282_000000000-DG2YB',
#   'stale': 'False',
#   'timestamp': '2022-04-19T13:43:09.757Z',
#   'tissue_origin': 'Pb',
#   'tissue_type': 'P',
#   'workflow': 'bcl2fastq',
#   'workflow_run_accession': 'vidarr:research/run/89fb17aed45e61d0caa9b23277fc973eb87309ec50de7469709daddfeb517a0e',
#   'workflow_run_attributes': '{}',
#   'workflow_run_labels': '{}',
#   'workflow_version': '[3,1,2]'}]



# =======

# i['pinery_data'] -> list of lists


# [{'barcode': 'TTCTGGAA-CGATTCCG',
#   'group_desc': '60',
#   'group_id': 'BMT-PMH-125',
#   'lane': '1',
#   'library_design': 'TS',
#   'library_name': 'AHCT_0029_01_LB02-01',
#   'library_type': '',
#   'lims': '{"id":"5884_1_LDI84239","provider":"pinery-miso-2.2","time":1676401595000,"version":"478c561cfe6b90bb6275692186e4cf9a87ea0c96de42f9954bcc77319ab67715"}',
#   'organism': 'Homo sapiens',
#   'path': '/oicr/data/archive/M06816/220412_M06816_0279_000000000-DG3HV',
#   'project': 'AHCT',
#   'run': '220412_M06816_0279_000000000-DG3HV',
#   'sex': 'None',
#   'start_date': '2022-04-12T00:00:00Z',
#   'tissue_origin': 'Pb',
#   'tissue_type': 'P'}]


# =====

# i['pinery_project_data']


# [[{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-2.2',
#    'sample_count': '416'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-stage-v8',
#    'sample_count': '367'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-v2',
#    'sample_count': '416'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-v7',
#    'sample_count': '416'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-v8',
#    'sample_count': '416'}]]




###############################################################

def connect_to_db(database):
    '''
    (str) -> sqlite3.Connection
    
    Returns a connection to SqLite database prov_report.db.
    This database contains information extracted from FPR
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    '''
    
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    return conn



def define_column_names():
    '''
    (None) -> dict

    Returns a dictionary with column names for each table in database
    '''

    # create dict to store column names for each table {table: [column names]}
    column_names = {'Workflows': ['wfrun_id', 'wf', 'wfv', 'project_id', 'case_id', 'attributes', 'file_count', 'lane_count', 'stale'],
                    'Parents': ['parents_id', 'children_id', 'project_id', 'case_id'],
                    'Projects': ['project_id', 'pipeline', 'last_updated', 'samples', 'library_types'],
                    'Files': ['file_swid', 'project_id', 'md5sum', 'workflow', 'version', 'wfrun_id', 'file', 'library_type', 'attributes', 'creation_date', 'limskey', 'stale', 'case_id'],
                    'FilesQC': ['file_swid', 'project_id', 'case_id', 'skip', 'user', 'date', 'status', 'reference', 'fresh', 'ticket'],
                    'Libraries': ['library', 'case_id', 'tissue_type', 'ext_id', 'tissue_origin',
                                  'library_type', 'group_id', 'group_id_description', 'project_id'],
                    'Workflow_Inputs': ['library', 'run', 'lane', 'wfrun_id', 'limskey', 'barcode', 'platform', 'project_id', 'case_id'],
                    'Samples': ['case_id', 'donor_id', 'species', 'sex', 'miso', 'created_date', 'modified_date', 'project_id', 'parent_project'],
                    'WGS_blocks': ['project_id', 'case_id', 'samples', 'anchor_wf', 'workflows', 'name', 'date', 'release_status', 'complete', 'clean', 'network'],
                    'WT_blocks': ['project_id', 'case_id', 'samples', 'anchor_wf', 'workflows', 'name', 'date', 'release_status', 'complete', 'clean', 'network'],
                    'Calculate_Contamination': ['sample_id', 'group_id', 'case_id', 'library_type', 'tissue_origin', 'tissue_type', 'contamination', 'merged_limskey'],
                    'Checksums': ['project_id', 'case_id', 'md5']
                    }
        
    return column_names


def define_column_types():
    '''
    (None) -> dict

    Returns a dictionary with column types for each table in database
    '''
    
    # create dict to store column names for each table {table: [column names]}
    column_types = {'Workflows': ['VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'TEXT', 'INT', 'INT', 'VARCHAR(128)'],
                    'Parents': ['VARCHAR(572)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)'],
                    'Projects': ['VARCHAR(128) PRIMARY KEY NOT NULL UNIQUE', 'VARCHAR(128)',
                                  'VARCHAR(256)', 'INT', 'INT'],
                    'Files': ['VARCHAR(572)', 'VARCHAR(128)',
                              'VARCHAR(256)', 'VARCHAR(128)', 'VARCHAR(128)',
                              'VARCHAR(572)', 'TEXT', 'VARCHAR(128)', 'TEXT', 'INT', 'VARCHAR(256)', 'VARCHAR(128)', 'VARCHAR(128)'],
                    'FilesQC': ['VARCHAR(572) PRIMARY KEY NOT NULL UNIQUE', 'VARCHAR(128)', 'VARCHAR(128)',
                                'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)'],
                    'Libraries': ['VARCHAR(256)', 'VARCHAR(128)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(128)'],
                    'Workflow_Inputs': ['VARCHAR(128)', 'VARCHAR(256)', 'INTEGER', 'VARCHAR(572)', 
                                        'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)'],
                    'Samples': ['VARCHAR(128) PRIMARY KEY NOT NULL', 'VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)'],
                    'WGS_blocks': ['VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(572)', 'TEXT', 'VARCHAR(256)', 'VARCHAR(128)', 'INT', 'INT', 'INT', 'TEXT'],
                    'WT_blocks': ['VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(572)', 'TEXT', 'VARCHAR(256)', 'VARCHAR(128)', 'INT', 'INT', 'INT', 'TEXT'],
                    'Calculate_Contamination': ['VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'FLOAT', 'VARCHAR(572)'],
                    'Checksums': ['VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)']
                    }
                    
    
    return column_types



def create_table(database, table):
    '''
    (str, str) -> None
    
    Creates a table in database
    
    Parameters
    ----------
    - database (str): Name of the database
    - table (str): Table name
    '''
    
    # get the column names
    column_names = define_column_names()[table]
    # get the column types
    column_types = define_column_types()[table]    
    
    # define table format including constraints    
    table_format = ', '.join(list(map(lambda x: ' '.join(x), list(zip(column_names, column_types)))))


    #if table  in ['Workflows', 'Parents', 'Files', 'FilesQC', 'Libraries', 'Workflow_Inputs', 'Samples', 'WGS_blocks', 'WT_blocks', 'Calculate_Contamination']:
    if table  in ['Workflows', 'Parents', 'Files', 'FilesQC', 'Libraries', 'Workflow_Inputs', 'Samples', 'WGS_blocks', 'WT_blocks', 'Checksums']:
        constraints = '''FOREIGN KEY (project_id)
            REFERENCES Projects (project_id)'''
        table_format = table_format + ', ' + constraints 
    
    if table == 'Parents':
        constraints = '''FOREIGN KEY (parents_id)
          REFERENCES Workflows (wfrun_id),
          FOREIGN KEY (children_id)
              REFERENCES Workflows (wfrun_id)''' 
        table_format = table_format + ', ' + constraints + ', PRIMARY KEY (parents_id, children_id, project_id, case_id)'
    
    if table == 'Worklows':
        table_format = table_format + ', PRIMARY KEY (wfrun_id, project_id)'
    
    if table in ['WGS_blocks', 'WT_blocks']:
        constraints = '''FOREIGN KEY (case_id)
          REFERENCES Samples (case_id)'''
        table_format = table_format + ', PRIMARY KEY (samples, anchor_wf)'
      
    if table == 'Files':
        constraints = '''FOREIGN KEY (wfrun_id)
            REFERENCES Workflows (wfrun_id),
            FOREIGN KEY (file_swid)
               REFERENCES FilesQC (file_swid)'''
        table_format = table_format + ', ' + constraints
    
    if table == 'Workflow_Inputs':
        constraints = '''FOREIGN KEY (wfrun_id)
            REFERENCES Workflows (wfrun_id),
            FOREIGN KEY (library)
              REFERENCES Libraries (library)'''
        table_format = table_format + ', ' + constraints
    
    if table == 'Samples':
        constraints = '''FOREIGN KEY (donor_id)
            REFERENCES Libraries (ext_id)'''
        table_format = table_format + ', ' + constraints

    if table == 'Libraries':
        constraints = '''FOREIGN KEY (case_id)
            REFERENCES Samples (case_id)'''
        table_format = table_format + ', ' + constraints

    # connect to database
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    # create table
    cmd = 'CREATE TABLE {0} ({1})'.format(table, table_format)
    cur.execute(cmd)
    conn.commit()
    conn.close()


def initiate_db(database):
    '''
    (str) -> None
    
    Create tables in database
    
    Parameters
    ----------
    - database (str): Path to the database file
    '''
    
    # check if table exists
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cur.fetchall()
    tables = [i[0] for i in tables]    
    conn.close()
    for i in ['Projects', 'Workflows', 'Parents', 'Files', 'FilesQC', 'Libraries',
              'Workflow_Inputs', 'Samples', 'WGS_blocks', 'WT_blocks', 
              'Calculate_Contamination', 'Checksums']:
        if i not in tables:
            create_table(database, i)




def insert_data(database, table, data, column_names):
    '''
    (str, str, list, list) -> None
    
    Inserts data into the database table with column names 
    
    Parameters
    ----------
    - database (str): Path to the database file
    - table (str): Table in database
    - data (list): List of data to be inserted
    - column_names (list): List of table column names
    '''
       

    # connect to db
    conn = sqlite3.connect(database)
    # add data
    vals = '(' + ','.join(['?'] * len(data[0])) + ')'
    conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), data)
    conn.commit()
    conn.close()



def delete_records(donors, database, table):
    '''
    (dict, str, str) -> None
    
    Remove all the rows from table with case_id in donors
    
    Parameters
    ----------
    - donors (dict): Dictionary with donors to remove from table
    - database (str): Path to the sqlite database
    - table (str): Table in database
    '''
    
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    id_list = list(donors.keys())
    query = "DELETE FROM {0} WHERE case_id IN ({1})".format(table, ", ".join("?" * len(id_list)))
    cur.execute(query, id_list)
    conn.commit()
    conn.close()



def add_checksums_info_to_db(database, donors_to_update, table = 'Checksums'):
    '''
    (str, dict, str) -> None
    
    Update table Checksums with the checksum of the donor info 
       
    Parameters
    ----------
    - database (str): Path to the database file
    - donors_to_update (dict): Dictionary with donors for which records needs to be updated
    - table (str): Name of Table in database. Default is Checksums
    '''
    
    if donors_to_update:
        delete_records(donors_to_update, database, table)
        
        # make a list of data to insert
        newdata = []
        
        # connect to db
        conn = sqlite3.connect(database, timeout=30)
        # get column names
        data = conn.execute("SELECT * FROM {0}".format(table))
        column_names = [column[0] for column in data.description]

        # order values according to column names
        for i in donors_to_update:
            L = [donors_to_update[i]['project_id'], i, donors_to_update[i]['md5']]
            newdata.append(L)
        vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
        conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
        conn.commit()
        conn.close()




def add_project_info_to_db(database, provenance_data, table = 'Projects'):
    '''
    (str, list, str) -> None
    
    Add project information into Projects table of database
       
    Parameters
    ----------
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - table (str): Name of Table in database. Default is Projects
    '''
    
    # collect project information
    project_info = collect_project_info(provenance_data)
    
    for project in project_info:
        # add time stamp to each project
        project_info[project]['last_updated'] = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
        # get the number of samples
        project_info[project]['samples'] = len(project_info[project]['samples'])
        # get the library types
        project_info[project]['library_types'] = ','.join(sorted(project_info[project]['library_types']))
        # add the project id
        project_info[project]['project_id'] = project
        
        
    # connect to db
    conn = sqlite3.connect(database, timeout=30)
    # get column names
    data = conn.execute("SELECT * FROM {0}".format(table))
    column_names = [column[0] for column in data.description]

    for project in project_info:
        # remove project info from table
        conn.execute('DELETE FROM {0} WHERE project_id = \"{1}\"'.format(table, project))
        conn.commit()
        # add project info in database
        # order values according to column names
        vals = [project_info[project][i] for i in column_names]
        # insert project info
        conn.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), tuple(vals)))
        conn.commit()
    
    conn.close()


def add_file_info_to_db(database, provenance_data, donors_to_update, table = 'Files'):
    '''
    (str, list, dict, str) -> None
    
    Inserts file information in database's Files table
       
    Parameters
    ----------
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - donors_to_update (dict): Dictionary with donors for which records needs to be updated
    - table (str): Table in database storing file information. Default is Files
    '''
        
    if donors_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(donors_to_update, database, table)
        print('deleted records in {0}'.format(table))
    
        # make a list of data to insert in the database
        newdata = []
            
        for donor_data in provenance_data:
            donor = donor_data['donor']
            # check if donor needs to be updated
            if donor in donors_to_update and donors_to_update[donor] != 'delete':
                file_info = collect_donor_file_info(donor_data)
                for file_swid in file_info:
                    file_info[file_swid]['limskey'] = ';'.join(sorted(list(set(file_info[file_swid]['limskey']))))
                    L = [file_info[file_swid][i] for i in column_names]
                    newdata.append(L)             
        
        # add data
        insert_data(database, table, newdata, column_names)
        

def add_library_info_to_db(database, provenance_data, donors_to_update, table = 'Libraries'):
    '''
    (str, list, dict, str) -> None
    
    Inserts library information in the Libraries table of the database
       
    Parameters
    ----------
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - donors_to_update (dict): Dictionary with donors for which records needs to be updated
    - table (str): Table in database storing file information. Default is Libraries
    '''                 

    if donors_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(donors_to_update, database, table)
        print('deleted records in {0}'.format(table))

        # make a list of data to insert in the database
        newdata = []
    
        for donor_data in provenance_data:
            donor = donor_data['donor']
            # check if donor needs to be updated
            if donor in donors_to_update and donors_to_update[donor] != 'delete':
                library_info = collect_donor_library_info(donor_data)
                
                for library in library_info:
                    for d in library_info[library]:
                        L = [d[i] for i in column_names]
                        newdata.append(L)             
  
        # add data
        insert_data(database, table, newdata, column_names)
        


# def add_fileQC_info_to_db(database, project, nabu_api, matched_ids, donors_to_update, table='FilesQC'):
#     '''
#     (str, str, str, dict, str, str) -> None
    
#     Inserts file QC information in database's FilesQC table
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - project (str): Name of project of interest
#     - nabu_api (str): URL of the nabu API
#     - matched_ids (dict): Dictionary of matched file swids and donor ids for each project in FPR
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): Table in database storing the file QC. Default is FilesQC
#     '''

#     # remove rows for donors to update
#     if donors_to_update:
#         delete_records(donors_to_update, database, table)
#         print('deleted records in FilesQC')

#         # collect QC info from nabu
#         D = collect_qc_info(project, database, nabu_api)
    
#         # check that data is recorded in nabu for project
#         if D:
#             # make a list of data to insert
#             newdata = []
#             # connect to db
#             conn = sqlite3.connect(database)
#             # get column names
#             column_names = define_column_names()[table]

#             # add data
#             for file_swid in D[project]:
#                 # check that file swid is recorded in FPR for the same project
#                 if file_swid in matched_ids[project]:
#                     if matched_ids[project][file_swid] in donors_to_update and donors_to_update[matched_ids[project][file_swid]] != 'delete':
#                         L = [D[project][file_swid][i] for i in column_names if i in D[project][file_swid]]
#                         L.insert(0, matched_ids[project][file_swid])
#                         L.insert(0, project)
#                         L.insert(0, file_swid)
#                         newdata.append(L)
#             vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#             conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#             conn.commit()
#             conn.close()




# def add_samples_info_to_db(database, project, pinery, table, donors_to_update, sample_info):
#     '''
#     (str, str, str, dict, dict, dict) -> None
    
#     Inserts samples data into Samples table of database    
    
#     Parameters
#     ----------
#     - database (str): Path to the databae file
#     - project (str): Name of project of interest
#     - pinery (str): Pinery API
#     - table (str): Name of table in database
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - sample_info (dict): Dictionary with sample information extracted from Pinery
#     '''
    
    
#     # remove rows for donors to update
#     if donors_to_update:
#         delete_records(donors_to_update, database, table)
#         print('deleted records in Samples')

#         # collect information about samples
#         samples = get_parent_sample_info(pinery, project, sample_info)
    
#         if samples:
#             # make a list of row data
#             newdata = []
#             # connect to db
#             conn = sqlite3.connect(database)
             
#             # get column names
#             data = conn.execute("SELECT * FROM {0} WHERE project_id = '{1}';".format(table, project))
#             column_names = [column[0] for column in data.description]

#             # add data into table
#             for i in samples:
#                 if i['case'] in donors_to_update and donors_to_update[i['case']] != 'delete':
#                    L = [i['case'], i['donor_id'], i['species'], i['sex'], i['miso'],
#                          i['created_date'], i['modified_date'], project, i['project']]          
#                    newdata.append(L)
            
#             vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#             conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#             conn.commit()
#             conn.close()








def add_workflows_to_db(database, provenance_data, donors_to_update, table = 'Workflows'):
    '''
    (str, list, dict, str) -> None
    
    Inserts or updates workflow information 
 
    Parameters
    ----------    
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - donors_to_update (dict): Dictionary with donors for which records needs to be updated
    - table (str): Table in database storing file information. Default is Workflows
    '''
    
    if donors_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(donors_to_update, database, table)
        print('deleted records in {0}'.format(table))

        # make a list of data to insert in the database
        newdata = []
        
        for donor_data in provenance_data:
            donor = donor_data['donor']
            # check if donor needs to be updated
            if donor in donors_to_update and donors_to_update[donor] != 'delete':
                workflow_info = collect_donor_workflow_info(donor_data)                
                for workflow in workflow_info:
                    L = [workflow_info[workflow][i] for i in column_names]
                    newdata.append(L)             
      
        # add data
        insert_data(database, table, newdata, column_names)
               
        

def add_workflow_inputs_to_db(database, provenance_data, donors_to_update, table = 'Workflow_Inputs'):
    '''
    (str, list, dict, str) -> None
    
    Inserts or updates workflow input information 
 
    Parameters
    ----------    
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - donors_to_update (dict): Dictionary with donors for which records needs to be updated
    - table (str): Table in database storing file information. Default is Workflow_Inputs
    '''

    if donors_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(donors_to_update, database, table)
        print('deleted records in {0}'.format(table))

        # make a list of data to insert in the database
        newdata = []
     
        for donor_data in provenance_data:
            donor = donor_data['donor']
            # check if donor needs to be updated
            if donor in donors_to_update and donors_to_update[donor] != 'delete':
                workflow_input_info = collect_donor_workflow_inputs(donor_data)
                for d in workflow_input_info:
                    L = [d[i] for i in column_names]
                    newdata.append(L)             

        # add data
        insert_data(database, table, newdata, column_names)
        


def add_workflows_relationships_to_db(database, provenance_data, donors_to_update, table = 'Parents'):
    '''
    (str, str, str, dict, str, str, str, str) -> None
    
    Inserts or updates workflow information and parent-children workflow relationships
 
    Parameters
    ----------   
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - donors_to_update (dict): Dictionary with donors for which records needs to be updated
    - table (str): Table in database storing file information. Default is Parents
    '''

    if donors_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(donors_to_update, database, table)
        print('deleted records in {0}'.format(table))
        
        # make a list of data to insert in the database
        newdata = []
        
        for donor_data in provenance_data:
            donor = donor_data['donor']
            project = get_project_name(donor_data)
            # check if donor needs to be updated
            if donor in donors_to_update and donors_to_update[donor] != 'delete':
                files = map_file_to_worklow(donor_data)
                workflow_inputs = get_workflow_inputs(donor_data)
                parents = identify_parent_children_workflows(workflow_inputs, files)
        
                # make a list of all workflows for the donor
                donor_workflows = list(set(list(files.values())))
                for workflow in donor_workflows:
                    for parent in parents[workflow]:
                        L = (os.path.basename(parent), os.path.basename(workflow), project, donor)
                        if L not in newdata:
                            newdata.append(L)
        
        # add data
        insert_data(database, table, newdata, column_names)
        


def is_project_active(donor_data):
    '''
    (dict) -> bool
    
    Returns True if the project is active and False otherwise
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data    
    '''
    
    active = [convert_to_bool(i['active']) for i in donor_data['pinery_project_data']]
    return all(active)



def remove_data_from_inactive_projects(provenance_data):
    '''
    (list) -> list

    Returns the list of donor data removing any data from inactive project
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    '''
    
    to_remove = [i for i in provenance_data if is_project_active(i) == False]
    if to_remove:
        for i in to_remove:
            provenance_data.remove(i)
    return provenance_data
    


def remove_donors_without_cerberus_data(provenance_data):
    '''
    (list) -> list
    
    Returns the list of donor data removing any donors that do not have data from cerberus
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    '''

    to_remove = [i for i in provenance_data if len(i['cerberus_data']) == 0]
    if to_remove:
        for i in to_remove:
            provenance_data.remove(i)
    return provenance_data


def compute_md5(d):
    '''
    (dict) -> str
    
    Returns the md5 checksum of a dictionary d
    
    Parameters
    ----------
    d (dict): Dictionary with information parsed from FPR
    '''
    
    return hashlib.md5(json.dumps(d, sort_keys=True).encode('utf-8')).hexdigest()


def compute_donor_md5sum(provenance_data):
    '''
    (list) -> dict
    
    Returns a dictionary with the md5sum of the each donor's data listed in
    provenance_data
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    '''

    D = {}
    for d in provenance_data:
        md5 = compute_md5(d)
        donor = d['donor']
        project = d['cerberus_data'][0]['project']
        assert donor not in D
        D[donor] = {'md5':md5, 'project_id':project}
    return D
    

def get_donors_md5sum(database, table = 'Checksums'):
    '''
    (str, str) -> dict

    Returns a dictionary of donors and checksum extracted from the table Checksums in the database
           
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - table (str): Table storing the donor checsum information. Default is Checksums 
    '''
        
    # connect to database, get recorded md5sums
    conn = connect_to_db(database)
    data = conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
    tables = [i['name'] for i in data]
    records = {}
    if table in tables:
        data = conn.execute('SELECT case_id, md5, project_id FROM {0}'.format(table)).fetchall()  
        for i in data:
            records[i['case_id']] = {'md5': i['md5'], 'project_id': i['project_id']}
    conn.close()
    
    return records
    

def donors_info_to_update(md5sums, recorded_md5sums):
    '''
    (dict, dict) -> dict

    Returns a dictionary of donors, checksum for which the information in the database needs to be updated
    (ie, the checksum in the database is different from the checksum of the production data,
     or the donor recorded in the database is no longer in production)
            
    Parameters
    ----------
    - md5sums (dict): Dictionary of donors, checksum for production data
    - recorded_md5sums (dict): Dictionary of recorded donors, checksum in the database
    '''
        
    donors = {}
    for donor in md5sums:
        # update if not already recorded
        if donor not in recorded_md5sums:
            donors[donor] = {'md5': md5sums[donor]['md5'], 'project_id': md5sums[donor]['project_id']}
        # update if md5sums are different
        else:
            if recorded_md5sums[donor]['md5'] != md5sums[donor]['md5']:
                donors[donor]['md5'] = md5sums[donor]['md5']
                assert donors[donor]['project_id'] == md5sums[donor]['project_id']
        
    # delete donors that are no longer recorded
    for donor in recorded_md5sums:
        if donor not in md5sums:
            donors[donor] = 'delete'
        
    return donors



def load_data(provenance_data_file):
    '''
    (str) -> list
    
    Returns the list of data contained in the provenance_data_file
    
    Parameters
    ----------
    - provenance_data_file (str): Path to the file with production data extracted from Shesmu
    '''

    infile = open(provenance_data_file, encoding='utf-8')
    provenance_data = json.load(infile)
    infile.close()
    
    return provenance_data



def get_project_name(donor_data):
    '''
    (dict) -> str
    
    Returns the name of the project extracted from the pinery project data of a donor 
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data    
    '''

    L = list(set([i['name'] for i in donor_data['pinery_project_data']]))
    assert len(L) == 1
    project = L[0]    
    return project


def get_pipeline(donor_data):
    '''
    (dict) -> str
    
    Returns the pipeline the project extracted from the pinery project data of a donor 
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data    
    '''

    L = list(set([i['pipeline'] for i in donor_data['pinery_project_data']]))
    assert len(L) == 1
    pipeline = L[0]    
    return pipeline
    



def get_donor_external_id(donor_data):
    '''
    (dict) -> str
    
    Returns the external id of a donor by extracting the external id from the
    cerberus data of that donor
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data    
    '''

    L = []

    for i in range(len(donor_data['cerberus_data'])):
        L.append(donor_data['cerberus_data'][i]['external_donor_id'])
    L = list(set(L))
    assert len(L) == 1
    return L[0]




def get_donor_sex(donor_data):
    '''
    (dict) -> str
    
    Returns the sex of a donor by extracting the sex from the
    pinery data of that donor
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data    
    '''
    
    sex = []
    for i in range(len(donor_data['pinery_data'])):
        sex.append(donor_data['pinery_data'][i]['sex'])       
    sex = list(set(sex))
    assert len(sex) == 1
    
    return sex[0]



def get_donor_species(donor_data):
    '''
    (dict) -> str
    
    Returns the species of a donor by extracting the species from the
    pinery data of that donor
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data    
    '''
    
    species = []
    for i in range(len(donor_data['pinery_data'])):
        species.append(donor_data['pinery_data'][i]['species'])       
    species = list(set(species))
    assert len(species) == 1
        
    return species[0]








def get_donor_library_types(donor_data):
    '''
    (dict) -> list
    
    Returns a list of 
    
    
    
    '''
    
    
    pass    
    
    
    
    
    
    
    
    

def get_donor_samples_and_library_types(donor_data):
    '''
    dict) -> dict

    Returns a dictionary mapping all samples with their corresponding library types
    for a given donor

    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data   
    '''
    
    D = {}
    for i in range(len(donor_data['cerberus_data'])):
        library_design = donor_data['cerberus_data'][i]['library_design']
        sample = '_'.join([donor_data['donor'],
                          donor_data['cerberus_data'][i]['tissue_type'], 
                          donor_data['cerberus_data'][i]['tissue_origin'],
                          library_design, donor_data['cerberus_data'][i]['group_id']])
        D[sample] = library_design
                           
    return D
            



def convert_to_bool(S):
    '''
    (str) -> bool
    
    Returns the boolean value of the string representation of a boolean
    
    Parameters
    ----------
    - S (str): String indicating True or False
    '''
    
    if S.lower() == 'true':
        B = True
    elif S.lower() == 'false':
        B = False
    return B
    
    
   
def get_file_timestamp(d):
    '''
    (dict) - > int
    
    Returns the time stamp of a file from the dictionary of the file from cerberus data
    for a donor
    
    Parameters
    ----------
    - d (dict): Dictionary representing a file information from cerberus
    '''    
    
    creation_date = d['timestamp']
    creation_date = ' '.join(creation_date.split('T')).replace('Z', '')
    creation_date = creation_date.split('.')
    if len(creation_date) > 1:
        creation_date = ' '.join(creation_date[:-1])
    else:
        creation_date = creation_date[0]
    pattern = '%Y-%m-%d %H:%M:%S'
    creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
    
    return creation_date
    
    
def collect_donor_file_info(donor_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the file information for a given donor
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data   
    '''    
        
    D = {}
    
    for i in range(len(donor_data['cerberus_data'])):
        file_swid = donor_data['cerberus_data'][i]['accession']
        project_id = donor_data['cerberus_data'][i]['project']
        md5sum = donor_data['cerberus_data'][i]['md5']
        workflow = donor_data['cerberus_data'][i]['workflow']
        case_id = donor_data['donor']
        file = donor_data['cerberus_data'][i]['path']
        library_type = donor_data['cerberus_data'][i]['library_design']
        stale = convert_to_bool(donor_data['cerberus_data'][i]['stale'])
        wfrun_id = os.path.basename(donor_data['cerberus_data'][i]['workflow_run_accession'])
        version = '.'.join(map(lambda x: str(x), json.loads(donor_data['cerberus_data'][i]['workflow_version'])))
        limskey = json.loads(donor_data['cerberus_data'][i]['lims'])['id']
        creation_date = get_file_timestamp(donor_data['cerberus_data'][i])
        if 'file_attributes' in donor_data['cerberus_data'][i]:
            attributes = json.loads(donor_data['cerberus_data'][i]['file_attributes'])
            for k in attributes:
                attributes[k] = attributes[k][0]
            if len(attributes) == 0:
                attributes = ''
            else:
                attributes = json.dumps(attributes)
        else:
            attributes = ''
    
        # collect file info if not stale             
        if stale == False:
            if file_swid not in D:
                D[file_swid] = {'file_swid': file_swid, 'project_id': project_id,
                                'md5sum': md5sum, 'workflow': workflow, 'case_id': case_id,
                                'file': file, 'library_type': library_type,
                                'stale': stale, 'wfrun_id': wfrun_id, 'version': version,
                                'limskey': [limskey], 'creation_date': creation_date,
                                'attributes': attributes}
            else:
                D[file_swid]['limskey'].append(limskey)
                
    return D
    


def collect_donor_library_info(donor_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the library information for a given donor
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data   
    '''    

    D = {}
    
    for i in range(len(donor_data['pinery_data'])):
        library = donor_data['pinery_data'][i]['library_name']
        case_id = donor_data['donor']
        ext_id = donor_data['cerberus_data'][0]['external_donor_id']
        project_id = donor_data['pinery_data'][i]['project']
        group_id = donor_data['pinery_data'][i]['group_id']
        group_id_description = donor_data['pinery_data'][i]['group_desc']
        library_type = donor_data['pinery_data'][i]['library_design']
        tissue_type = donor_data['pinery_data'][i]['tissue_type']
        tissue_origin = donor_data['pinery_data'][i]['tissue_origin']
        
        d = {'library': library, 'case_id': case_id, 'ext_id': ext_id,
                      'project_id': project_id, 'group_id': group_id,
                      'group_id_description': group_id_description,
                      'library_type': library_type, 'tissue_type': tissue_type,
                      'tissue_origin': tissue_origin}
        
        if library not in D:
            D[library] = [d] 
        else:
            if d not in D[library]:
                D[library].append(d)
                   
    return D       
        

        



def collect_donor_file_qc_info(donor_data):
    '''
    
    
    
    

    Parameters
    ----------
    donor_info : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    D = {}
        
    for i in range(len(donor_data['nabu_data'])):
        ticket = donor_data['nabu_data'][i]['comment']
        file_swid = donor_data['nabu_data'][i]['file_id']
        project = donor_data['nabu_data'][i]['project']
        user = donor_data['nabu_data'][i]['username']
        
        
        pass
        
                # 'case_id',
                # 'skip',
                # 'date',
                # 'status',
                # 'reference',
                # 'fresh',
                



#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25619611/AHCT_0029_01_LB01-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-TACTGCAG_R2.fastq.gz',
#    'file_swid': '25620329',
#    'fileqc_id': '58037',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',

#    'workflow': 'bcl2fastq'}],
#  [{'comment': 'GDR-839',
#    'file_id': 'vidarr:research/file/31c12c06d91562e28142e4e95a394b6d9fc769ced2a2216d5e904096bb2aff04',
#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25618917/AHCT_0029_01_LB02-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-CGATTCCG_R2.fastq.gz',
#    'file_swid': '25619744',
#    'fileqc_id': '57741',
#    'project': 'AHCT',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',
#    'username': 'rjovelin',
#    'workflow': 'bcl2fastq'}],
#  [{'comment': 'GDR-839',
#    'file_id': 'vidarr:research/file/87056488f0a3a34da37114965e438008bacf7a1e2ca44a741c40eeecc1d81df0',
#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25619611/AHCT_0029_01_LB01-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-TACTGCAG_R1.fastq.gz',
#    'file_swid': '25620299',
#    'fileqc_id': '57996',
#    'project': 'AHCT',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',
#    'username': 'rjovelin',
#    'workflow': 'bcl2fastq'}],
#  [{'comment': 'GDR-839',
#    'file_id': 'vidarr:research/file/d586a196c2ed7bb95daaeb5eda6f01cdc91dc710d42b388eb24a2ede14257c02',
#    'file_path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25618917/AHCT_0029_01_LB02-01_220503_A00469_0305_AHYMN2DRXY_1_TTCTGGAA-CGATTCCG_R1.fastq.gz',
#    'file_swid': '25619727',
#    'fileqc_id': '57569',
#    'project': 'AHCT',
#    'qc_date': '2022-10-28T17:16:11.457Z',
#    'qc_passed': 'True',
#    'username': 'rjovelin',
#    'workflow': 'bcl2fastq'}]]



def collect_donor_sample_info(donor_data):
    '''
    
    
    '''
    
    
    D = {}
        
    for i in range(len(donor_data['nabu_data'])):
        case_id = donor_data['donor']       
        donor_id = get_donor_external_id(donor_data)
        sex = get_donor_sex(donor_data)
        species = get_donor_species(donor_data)
        project_id = get_project_name(donor_data)
    
         
    
    # 'miso',
    # 'created_date',
    # 'modified_date',
    # 'parent_project'
    








# i['cerberus_data'] -> list of lists

# [{'accession': 'vidarr:research/file/012cf6554e9b51f0dcab3186a134f02b26b78953534e56dbd14ba3158428c380',
#   'barcode': 'TTCTGGAA-CGATTCCG',
#   'external_donor_id': 'BMT-PMH-125',
#   'group_id': 'BMT-PMH-125',
#   'input_files': '[]',
#   'instrument_model': 'Illumina MiSeq',
#   'lane': '1',
#   'library_design': 'TS',
#   'library_name': 'AHCT_0029_01_LB02-01',
#   'library_type': '',
#   'lims': '{"id":"5893_1_LDI84239","provider":"pinery-miso","time":1676401595000,"version":"882752c33378b9d66793a619fc4918241285391bf490f0525325f90fd58db05b"}',
#   'md5': 'c55a9d4efd25ff5d081b4c3301453561',
#   'path': '/oicr/data/archive/seqware/seqware_analysis_12/hsqwprod/seqware-results/bcl2fastq_3.1.2/25547692/AHCT_0029_01_LB02-01_220418_M06816_0282_000000000-DG2YB_1_TTCTGGAA-CGATTCCG_R2.fastq.gz',
#   'project': 'AHCT',
#   'run': '220418_M06816_0282_000000000-DG2YB',
#   'stale': 'False',
#   'timestamp': '2022-04-19T13:43:09.757Z',
#   'tissue_origin': 'Pb',
#   'tissue_type': 'P',
#   'workflow': 'bcl2fastq',
#   'workflow_run_accession': 'vidarr:research/run/89fb17aed45e61d0caa9b23277fc973eb87309ec50de7469709daddfeb517a0e',
#   'workflow_run_attributes': '{}',
#   'workflow_run_labels': '{}',
#   'workflow_version': '[3,1,2]'}]



# =======

# i['pinery_data'] -> list of lists


# [{'barcode': 'TTCTGGAA-CGATTCCG',
#   'group_desc': '60',
#   'group_id': 'BMT-PMH-125',
#   'lane': '1',
#   'library_design': 'TS',
#   'library_name': 'AHCT_0029_01_LB02-01',
#   'library_type': '',
#   'lims': '{"id":"5884_1_LDI84239","provider":"pinery-miso-2.2","time":1676401595000,"version":"478c561cfe6b90bb6275692186e4cf9a87ea0c96de42f9954bcc77319ab67715"}',
#   'organism': 'Homo sapiens',
#   'path': '/oicr/data/archive/M06816/220412_M06816_0279_000000000-DG3HV',
#   'project': 'AHCT',
#   'run': '220412_M06816_0279_000000000-DG3HV',
#   'start_date': '2022-04-12T00:00:00Z',
#   'tissue_origin': 'Pb',
#   'tissue_type': 'P'}]


# =====

# i['pinery_project_data']


# [[{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-2.2',
#    'sample_count': '416'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-stage-v8',
#    'sample_count': '367'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-v2',
#    'sample_count': '416'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-v7',
#    'sample_count': '416'}],
#  [{'active': 'True',
#    'name': 'YOCRC',
#    'pipeline': 'Accredited',
#    'provider': 'pinery-miso-v8',
#    'sample_count': '416'}]]


  
    

def collect_donor_workflow_info(donor_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the workflow information for a given donor
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data   
    '''

    D = {}
    
    for i in range(len(donor_data['cerberus_data'])):
        project_id = donor_data['cerberus_data'][i]['project']
        wfrun_id = os.path.basename(donor_data['cerberus_data'][i]['workflow_run_accession'])
        wf = donor_data['cerberus_data'][i]['workflow']
        wfv = '.'.join(map(lambda x: str(x), json.loads(donor_data['cerberus_data'][i]['workflow_version'])))
        case_id = donor_data['donor']
        stale = convert_to_bool(donor_data['cerberus_data'][i]['stale'])
        file_swid = donor_data['cerberus_data'][i]['accession']     
        limskey = json.loads(donor_data['cerberus_data'][i]['lims'])['id']
        attributes = json.loads(donor_data['cerberus_data'][i]['workflow_run_attributes'])
        attributes = json.dumps(attributes)
        
        d = {'project_id': project_id, 'wfrun_id': wfrun_id, 'wf': wf,
             'wfv': wfv, 'case_id': case_id, 'stale': stale, 'file_count': [file_swid],
             'lane_count': [limskey], 'attributes': attributes}
             
        if wfrun_id not in D:
            D[wfrun_id] = d
        else:
            D[wfrun_id]['file_count'].append(file_swid)
            D[wfrun_id]['lane_count'].append(limskey)
        
    for i in D:
        D[i]['file_count'] = len(list(set(D[i]['file_count'])))
        D[i]['lane_count'] = len(list(set(D[i]['lane_count'])))
        
    return D            




def collect_donor_workflow_inputs(donor_data):
    '''
    (dict) -> list
    
    Returns a list of dictionaries with the workflow input information for a given donor
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data   
    '''
        
    L = []
    
    for i in range(len(donor_data['cerberus_data'])):
        case_id = donor_data['donor']
        library = donor_data['cerberus_data'][i]['library_name']
        project_id = donor_data['cerberus_data'][i]['project']
        barcode = donor_data['cerberus_data'][i]['barcode']
        platform = '_'.join(donor_data['cerberus_data'][i]['instrument_model'].split())
        lane = donor_data['cerberus_data'][i]['lane']
        wfrun_id = os.path.basename(donor_data['cerberus_data'][i]['workflow_run_accession'])
        run = donor_data['cerberus_data'][i]['run']
        limskey = json.loads(donor_data['cerberus_data'][i]['lims'])['id']
        
        d = {'case_id': case_id, 'library': library, 'project_id': project_id,
             'barcode': barcode, 'platform': platform, 'lane': lane, 'wfrun_id': wfrun_id,
             'run': run, 'limskey': limskey}
        
        if d not in L:
            L.append(d)
    
    return L










def collect_project_info(provenance_data):
    '''
    (list) -> dict
    
    Returns a dictionary with project level information for each active project in provenance_data
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    '''


    D = {}


    for i in range(len(provenance_data)):
        project_name = get_project_name(provenance_data[i])
        if project_name not in D:
            D[project_name] = {'samples': [], 'library_types': []}
        samples = get_donor_samples_and_library_types(provenance_data[i])
        pipeline = get_pipeline(provenance_data[i])
        D[project_name]['pipeline'] = pipeline
        D[project_name]['samples'].extend(list(samples.keys()))
        D[project_name]['library_types'].extend(list(samples.values()))
        D[project_name]['samples'] = list(set(D[project_name]['samples']))
        D[project_name]['library_types'] = list(set(D[project_name]['library_types']))
    
    return D    
        
        
        
        
        
        
def map_file_to_worklow(donor_data):
    '''
    (dict) -> dict

    Returns a dictionary of file swids matched to their workflow run id for a donor    
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data   
    '''
    
    D = {}
       
    for i in range(len(donor_data['cerberus_data'])):
        wfrun_id = donor_data['cerberus_data'][i]['workflow_run_accession']
        file_swid = donor_data['cerberus_data'][i]['accession'] 
        if file_swid in D:
            assert D[file_swid] == wfrun_id
        else:
            D[file_swid] = wfrun_id
    
    return D
    
            
def get_workflow_inputs(donor_data):
    '''
    (dict) -> dict

    Returns a dictionary of workflows and their input files for a donor    
    
    Parameters
    ----------
    - donor_data (dict): Dictionary with a single donor data   
    '''
    
    D = {}
       
    for i in range(len(donor_data['cerberus_data'])):
        wfrun_id = donor_data['cerberus_data'][i]['workflow_run_accession']
        input_files = json.loads(donor_data['cerberus_data'][i]['input_files']) 
        if wfrun_id in D:
            assert D[wfrun_id] == input_files
        else:
            D[wfrun_id] = input_files
    
    return D
        
        
def identify_parent_children_workflows(workflow_inputs, files):
    '''
    (dict, dict) -> dict     
    
    Returns a dictionary of children: parents workflows relationsips for a donor
        
    Parameters
    ----------
    - workflow_inouts (dict): Dictionary of workflows and their input files for a donor
    - files (dict): Dictionary of file swids matched to their workflow run id for a donor    
    '''
    
    # parents record child-parent workflow relationships
    D = {}
    
    for workflow in workflow_inputs:
        if workflow_inputs[workflow]:
            parent_workflows = sorted(list(set([files[i] for i in workflow_inputs[workflow] if i in files])))
        else:
            parent_workflows = ['NA']
        if workflow not in D:
            D[workflow] = parent_workflows
        else:
            assert D[workflow] == parent_workflows
    
    return D
        



def generate_database(database, provenance_data_file):
    '''
    
    
    '''
    
    
    # create database if file doesn't exist
    if os.path.isfile(database) == False:
        initiate_db(database)
    print('initiated database')
    
    # load data from file
    provenance_data = load_data(provenance_data_file)
    print('loaded data')
    
    # clean up data 
    # remove data from inactive projects
    provenance_data = remove_data_from_inactive_projects(provenance_data)
    print('removed inactive projects')
    provenance_data = remove_donors_without_cerberus_data(provenance_data)
    print('removed donors without cerberus data')

    # collect the md5sum of each donor's data
    md5sums = compute_donor_md5sum(provenance_data)
    print('computed md5sums')
    # collect the recorded md5sums of the donor data from the database
    recorded_md5sums = get_donors_md5sum(database, table = 'Checksums')
    print('pulled md5sums from database')
    # determine the donors that need an update
    donors_to_update = donors_info_to_update(md5sums, recorded_md5sums)
    print('determined donors to update')
    
    # add project information
    add_project_info_to_db(database, provenance_data, 'Projects')
    print('added project info to database')
    
    # add file information to database
    add_file_info_to_db(database, provenance_data, donors_to_update, 'Files')
    print('added file info to database')
    
    add_library_info_to_db(database, provenance_data, donors_to_update, 'Libraries')
    print('added library info to database')
    
    
    add_workflows_to_db(database, provenance_data, donors_to_update, 'Workflows')
    print('added workflow info to database')
         
    
    
    add_workflow_inputs_to_db(database, provenance_data, donors_to_update, 'Workflow_Inputs')
    print('added workflow inputs to database')
    
    
    add_workflows_relationships_to_db(database, provenance_data, donors_to_update, 'Parents')
    print('added workflow relationships to database')
    



    
    # update the checksums for donors
    add_checksums_info_to_db(database, donors_to_update, 'Checksums')
    print('added md5sums to database')
    
    



generate_database('test2.db', 'provenance_reporter.json')    















    
    
        
        
        
            
                






#############################################################





# def extract_project_info(pinery):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with project information pulled down from Pinary API

#     Parameters
#     ----------
#     - pinery (str): Pinery API, http://pinery.gsi.oicr.on.ca
#     '''
    
#     project_provenance = pinery + '/sample/projects'
    
#     headers = {'accept': 'application/json',}
#     response = requests.get(project_provenance, headers=headers)
    
#     if response.ok:
#         L = response.json()
#     else:
#         L = []
    
#     D = {}
    
#     if L:
#         for i in L:
#             name = i['name']
#             assert name not in D
#             D[name] = {'project_id': name}
#             for j in ['pipeline', 'description', 'active', 'contact_name', 'contact_email']:
#                 if j in i:
#                     D[name][j] = i[j] 
#                 else:
#                     D[name][j] = ''
#                 if j == 'active' and j in i:
#                     if i[j]:
#                         D[name][j] = 'Active'
#                     else:
#                         D[name][j] = 'Completed'
#     return D                            


# def collect_info(data, names, keys):
#     '''
#     (dict, list, list) -> dict

#     Returns a dictionary with keys by extracting specific information using names from data
#     Note, the value of any name in data may be a string or single-element-list

#     Parameters
#     ----------
#     - data (dict): Dictionary with data to be extracted
#     - names (list): List of keys of interest in data
#     - keys (list): List of renamed keys in output dictionary
#     '''
    
#     d = {}
#     for j in range(len(names)):
#         if names[j] in data:
#             if type(data[names[j]]) == list:
#                 d[keys[j]] = data[names[j]][0]         
#             else:
#                 d[keys[j]] = data[names[j]]
#         else:
#             d[keys[j]] = ''

#     return d




# def get_QC_status_from_nabu(project, workflow, api):
#     '''
#     (str, str, str) -> dict
    
#     Returns a dictionary with qc information extracted from nabu for files generated by workflow in project
        
#     Parameters
#     ----------
#     - project (str): Project of interest
#     - workflow (str): Workflow of interest
#     - api (str): URL of the Nabu api
#     '''

#     # get end-point
#     api += 'get-fileqcs' if api[-1] == '/' else '/get-fileqcs'
    
#     D = {project: {}}
    
#     # check each fastq-generating workflow
#     headers = {'accept': 'application/json','Content-Type': 'application/json'}
#     json_data = {"project": "{0}".format(project), "workflow": workflow}
#     response = requests.post(api, headers=headers, json=json_data)
#     # check response code
#     if response.status_code == 200:
#         L = response.json()['fileqcs']
#         if L:
#             for i in L:
#                 qc = collect_info(i, ['skip', 'username', 'date', 'qcstatus', 'ref', 'stalestatus', 'comment'], ['skip', 'user', 'date', 'status', 'reference', 'fresh', 'ticket']) 
#                 file_swid = i['fileid']    
#                 D[project][file_swid] = qc
#     return D    
    


# def get_project_workflows(project, database, workflow_table = 'Workflows'):
#     '''
#     (str, str, str) -> list
    
#     Returns a list of all workflows for a given project

#     Parameters
#     ----------
#     - project (str): Project of interest
#     - database (str): Path to the sqlite database
#     - workflow_table (str): Name of table with workflow information. Default: Workflows
#     '''

#     # make a list of all workflows for a given project
#     conn = sqlite3.connect(database)
#     conn.row_factory = sqlite3.Row
#     data = conn.execute('SELECT * FROM {0} WHERE project_id="{1}"'.format(workflow_table, project)).fetchall()
#     conn.close()
#     workflows = set()
#     for i in data:
#         i = dict(i)
#         workflows.add(i['wf'])
#     workflows = list(workflows)
    
#     return workflows


# def add_missing_QC_status(D, project, database, table = 'Files'):
#     '''
#     (dict, str, str, str) -> dict
    
#     - D (dict): Dictionary with QC information for project files
#     - project (str): Name of project of interest
#     - database (str): Path to the database file
#     - table (str): Table in database storing file information. Default is Files
#     - nabu_api (str): URL of the nabu API
#     '''
    
#     # QC may not be available for all files
#     # retrieve file swids from database instead of parsing again fpr
#     # add empty values to qc fields if qc not available
#     conn = sqlite3.connect(database)
#     cur = conn.cursor()
#     cur.execute('SELECT {0}.file_swid FROM {0} WHERE {0}.project_id = \"{1}\"'.format(table, project))
#     records = cur.fetchall()
#     records = [i[0] for i in records]
#     conn.close()
    
#     for file_swid in records:
#         if file_swid not in D[project]:
#             D[project][file_swid] = {'skip': '', 'user': '', 'date': '', 'status': '', 'reference': '', 'fresh': '', 'ticket': ''}
    
#     return D
    


# def match_donor_to_file_swid(fpr_data, projects):
#     '''
#     (str, dict) -> dict
    
#     Returns a dictionary matching each file swid for each valid project to its donor_id
    
#     Parameters
#     ----------
#     - fpr_data (dict): Dictionary with file information extracted from File Provenance Report file
#     - projects (list): List of projects 
#     '''
    
#     D = {}
#     for project in fpr_data:
#         for case_id in fpr_data[project]:
#             for file_swid in fpr_data[project][case_id]:
#                 if project not in D:
#                     D[project] = {}
#                 D[project][file_swid] = case_id
#     return D            
        
    


# def collect_qc_info(project, database, nabu_api):
#     '''
#     (str, str, str, str) -> dict

#     Returns a dictionary with QC status of all files in project    
        
#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - database (str): Path to the database file
#     - nabu_api (str): URL of the nabu API
#     '''

#     # track qc info for all files in project
#     D = {project: {}}

#     # make a list of workflows for project
#     workflows = get_project_workflows(project, database)
#     for workflow in workflows:
#         qc = get_QC_status_from_nabu(project, workflow, nabu_api)     
#         # update dict
#         D[project].update(qc[project])
#     # add project files that may be missing QC info in Nabu
#     D = add_missing_QC_status(D, project, database)

#     return D





# def donors_info_to_update(database, fpr_data, project, table = 'Checksums'):
#     '''
#     (str, str, str) -> dict

#     Returns a dictionary of donors, checksum for projhect of interest for which
#     the information extracted from FPR is different than the information stored in the database.
#     The information for these donors need to be updated
        
#     Parameters
#     ----------
#     - database (str): Path to the sqlite database
#     - fpr_data (dict): Dictionary with file information extracted from FPR and organized by donor
#     - project (str): Name of project of interest 
#     '''
        
#     # evaluate the checksum of extracted data
#     md5 = {}
#     for case_id in fpr_data[project]:
#         md5sum = get_md5(fpr_data[project][case_id])
#         assert case_id not in md5
#         md5[case_id] = md5sum
    
#     # connect to database, get recorded md5sums
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
#     tables = [i['name'] for i in data]
#     records = {}
#     if table in tables:
#         data = conn.execute('SELECT case_id, md5 FROM {0} WHERE project_id = \"{1}\"'.format(table, project)).fetchall()  
#         for i in data:
#             records[i['case_id']] = i['md5']
#     conn.close()
    
#     donors = {}
#     for case in md5:
#         # update if not already recorded
#         if case not in records:
#             donors[case] = md5[case]
#         # update if md5sums are different
#         else:
#             if records[case] != md5[case]:
#                 donors[case] = md5[case]
        
#     # delete records no longer in fpr
#     for i in records:
#         if i not in md5:
#             donors[i] = 'delete'
        
#     return donors



# def collect_file_info_from_fpr(fpr, projects):
#     '''
#     (str, str) -> dict

#     Returns a dictionary with file information for project_name extracted from File Provenance Report
        
#     Parameters
#     ----------
#     - fpr (str): Path to File Provenance Report file
#     - project_name (str): Name of project of interest 
#     '''

#     infile = open_fpr(fpr)
#     # skip header
#     infile.readline()
        
#     fpr_data = {}
    
#     for line in infile:
#         line = line.rstrip()
#         if line:
#             line = line.split('\t')
#             # get project and initiate dict
#             project = line[1]
#             if project in projects:
#                 if project not in fpr_data:
#                     fpr_data[project] = {}
                   
#                 # skip stale records
#                 stale = line[52]
#                 if stale.lower() != 'stale':
#                     case_id = line[7]
#                     # get creation date
#                     creation_date = line[0]
#                     # convert creation date into epoch time
#                     # remove milliseconds
#                     creation_date = creation_date.split('.')[0]
#                     pattern = '%Y-%m-%d %H:%M:%S'
#                     creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
#                     # get file path
#                     file = line[46]
#                     # get md5sums
#                     md5sum = line[47]
#                     # get file_swid
#                     file_swid = line[44]
            
#                     # get the library type
#                     library_type = line[17]
#                     if library_type:
#                         d = collect_info({k.split('=')[0]:k.split('=')[1] for k in line[17].split(';')},
#                                  ['geo_library_source_template_type'], ['library_type'])
#                         library_type = d['library_type']
#                     else:
#                         library_type = ''
                        
#                     # get file attributes
#                     attributes = line[45]
#                     if attributes:
#                         attributes = attributes.split(';')
#                         attributes = json.dumps({k.split('=')[0]: k.split('=')[1] for k in attributes})
#                     else:
#                         attributes = ''
            
#                     # get workflow attributes
#                     wf_attributes = line[37]
#                     if wf_attributes:
#                         wf_attributes = wf_attributes.split(';')
#                         wf_attributes = json.dumps({k.split('=')[0]: k.split('=')[1] for k in wf_attributes if k.split('=')[0] not in ['cromwell-workflow-id', 'major_olive_version']})
#                     else:
#                         wf_attributes = ''                
       
#                     # get workflow
#                     workflow = line[30]
#                     # get workflow version
#                     version = line[31]        
#                     # get workflow run accession
#                     workflow_run = line[36]
#                     # get limskey
#                     limskey = line[56]
                    
#                     skip = line[51]
#                     if skip.lower() == 'true':
#                         skip = 1
#                     else:
#                         skip = 0
                    
#                     if case_id not in fpr_data[project]:
#                         fpr_data[project][case_id] = {}
                
#                     # collect file info
#                     if file_swid not in fpr_data[project][case_id]:
#                         fpr_data[project][case_id][file_swid] = {'creation_date': creation_date, 'md5sum': md5sum,
#                                                  'workflow': workflow, 'version': version,
#                                                  'wfrun_id': workflow_run, 'file': file,
#                                                  'library_type': library_type, 'attributes': attributes,
#                                                  'wf_attributes': wf_attributes, 'limskey': [limskey],
#                                                  'skip': skip, 'stale': stale, 'case_id': case_id}
#                     else:
#                         fpr_data[project][case_id][file_swid]['limskey'].append(limskey)
    
#     infile.close()
#     return fpr_data



# def extract_workflow_info(fpr, projects):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with library input information for all workflows for project_name
    
#     Parameters
#     ----------
#     - fpr (str): Path to the File Provenance Report
#     - project_name (str): Name of project of interest
#     '''

#     # create a dict to store workflow input libraries
#     D = {}

#     infile = open_fpr(fpr)
#     # skip header
#     infile.readline()
    
#     for line in infile:
#         line = line.rstrip()
#         if line:
#             line = line.split('\t')
#             # get project name
#             project = line[1]
#             if project in projects:
#                 # do not include stale records
#                 stale = line[52]
#                 if stale.lower() != 'stale':
#                     # get sample name
#                     sample = line[7]
#                     # get workflow name and workflow run accession
#                     workflow, workflow_run = line[30], line[36]
                                  
#                     # get lane and run
#                     run, lane = line[18], line[24]            
#                     # get library and limskey
#                     library, limskey  = line[13], line[56]
            
#                     # get barcode and platform
#                     barcode, platform = line[27], line[22]
            
#                     d = {'run': run, 'lane': lane, 'library': library, 'limskey': limskey, 'barcode': barcode, 'platform': platform}
            
#                     if project not in D:
#                         D[project] = {}
#                     if workflow_run not in D[project]:
#                         D[project][workflow_run] = {}
#                     if sample not in D[project][workflow_run]:
#                         D[project][workflow_run][sample] = {'sample': sample, 'workflow': workflow, 'libraries': [d]}
#                     else:
#                         assert sample == D[project][workflow_run][sample]['sample']
#                         assert workflow == D[project][workflow_run][sample]['workflow']
#                         if d not in D[project][workflow_run][sample]['libraries']:
#                             D[project][workflow_run][sample]['libraries'].append(d)
    
#     return D            



# def extract_projects_from_fpr(fpr):
#     '''
#     (str) -> list

#     Returns a list of projects with data in FPR
        
#     Parameters
#     ----------
#     - fpr (str): Path to File Provenance Report file
#     '''

#     infile = open_fpr(fpr)
#     # skip header
#     infile.readline()
        
#     projects = set()
    
#     for line in infile:
#         line = line.rstrip()
#         if line:
#             line = line.split('\t')
#             projects.add(line[1])
#     infile.close()
#     projects = list(projects)
    
#     return projects



# def get_valid_projects(projects, fpr):
#     '''
#     (dict, str) -> list
    
#     Returns a list of projects defined in Pinery with data available in FPR
    
#     Parameters
#     ----------
#     - projects (dict): Dictionary with project information extracted from project provenance
#     - fpr (str): Path to the File Provenance Report
#     '''
    
#     # make a list of projects available in Pinery
#     P = list(projects.keys())
#     # list all projects with data in FPR
#     L = extract_projects_from_fpr(fpr)
#     # keep only projects listed in Pinery with data in FPR
#     projects = list(set(P).intersection(set(L)))
   
#     return projects


# def filter_completed_projects(projects):
#     '''
#     (dict) -> dict
   
#     Returns a dictionary in which only Active projects are kept (ie, Completed projects are removed)
   
#     Parameters
#     ----------
#     projects (dict): Dictionary with project information extracted from project provenance
#     '''

#     to_remove = [i for i in projects if projects[i]['active'].lower() != 'active'] 
#     for i in to_remove:
#         del projects[i]
#     return projects





# def match_file_worklow_ids(fpr_data, project):
#     '''
#     (dict, str) -> dict

#     Returns a dictionary of file swids matched to their workflow run id for a project of interest    
    
#     Parameters
#     ----------
#     - fpr_data (dict): Dictionary with file information parsed from FPR
#     - project (str): Name of the project of interest
#     '''
    
#     D = {}
#     D[project] = {}
    
#     for case_id in fpr_data[project]:
#         for file_swid in fpr_data[project][case_id]:
#             workflow_id = fpr_data[project][case_id][file_swid]['wfrun_id']
#             if file_swid in D[project]:
#                 assert D[project][file_swid] == workflow_id
#             else:
#                 D[project][file_swid] = workflow_id
#     return D        



# def get_workflow_information(fpr_data, donors_to_update, project, database, file_table='Files', workflow_input_table='Workflow_Inputs'):
#     '''
#     (dict, dict, str, str, str, str) -> dict    
    
#     Returns a dictinary with workflow information for the project of interest
        
#     Parameters
#     ----------
#     - fpr_data (dict): Dictionary with file information parsed from FPR
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - project (str): Name of the project of interest
#     - database (str): Path to the sqlite database
#     - file_table (str): Table in database storing file information
#     - workflow_input_table (str): Table in database storing the worflow input information
#     '''
    
#     D = {}
#     D[project] = {}
    
#     if donors_to_update:
#         for case_id in donors_to_update:
#             for file_swid in fpr_data[project][case_id]:
#                 workflow = fpr_data[project][case_id][file_swid]['workflow']
#                 version = fpr_data[project][case_id][file_swid]['version']
#                 workflow_id = fpr_data[project][case_id][file_swid]['wfrun_id']
#                 attributes = fpr_data[project][case_id][file_swid]['wf_attributes']
#                 skip = fpr_data[project][case_id][file_swid]['skip']
#                 stale = fpr_data[project][case_id][file_swid]['stale']
            
#                 D[project][workflow_id] = {'wfrun_id': workflow_id,
#                                            'wf': workflow,
#                                            'wfv': version,
#                                            'attributes': attributes,
#                                            'project_id': project,
#                                            'case_id': case_id,
#                                            'skip': skip,
#                                            'stale': stale}
    
#         # get the file count for each workflow
#         counts = count_files(project, database, file_table)
#         # update workflow information with file count
#         update_workflow_information(project, D, counts, 'file_count')
#         # get the amount of data for each workflow
#         limskeys = get_workflow_limskeys(project, database, workflow_input_table)
#         for i in limskeys:
#             limskeys[i] = len(limskeys[i])
#         # update workflow information with lane count
#         update_workflow_information(project, D, limskeys, 'lane_count')
    
#     return D


# def get_parent_file_ids(fpr, projects):
#     '''
#     (str, list) -> dict

#     Returns a dictionary with parent file ids for each workflow in each project 
    
#     Parameters
#     ----------
#     - fpr (str): Path to the File Provenance Report
#     - projects (list): List of valid projects
#     '''

#     P = {}

#     # open fpr
#     infile = open_fpr(fpr)
#     # skip header
#     infile.readline()
       
#     for line in infile:
#         line = line.rstrip()
#         if line:
#             line = line.split('\t')
#             # get project name
#             project = line[1]
#             if project in projects:
#                 stale = line[52]
#                 # do not include stale records
#                 if stale.lower() != 'stale':
#                     # get workflow, workflow version and workflow run accession
#                     workflow_run = line[36]
#                     # get input files          
#                     input_files = line[38]
#                     if input_files:
#                         input_files = sorted(input_files.split(','))
#                     else:
#                         input_files = []
                           
#                     if project not in P:
#                         P[project] = {}
#                     if workflow_run in P[project]:
#                         assert P[project][workflow_run] == input_files
#                     else:
#                         P[project][workflow_run] = input_files
     
#     infile.close()
    
#     return P         


    
# def get_provenance_data(provenance):
#     '''
#     (str) -> list
    
#     Returns a list of dictionary with lims information for each library
    
#     Parameters
#     ----------
#     - provenance (str): URL of the pinery provenance API 
#     '''
    
#     #response = requests.get(provenance, timeout = (10,120))
#     response = requests.get(provenance, timeout = (20,240))
#     if response.ok:
#         L = response.json()
#     else:
#         L = []
    
#     return L



# def get_sample_info(pinery, project):
#     '''
#     (str, str, str) -> dict
    
#     Return a dictionary with sample information, including samples not sequenced
#     for each project
    
#     Parameters
#     ----------
#     - pinery (str): Pinery API
#     - project (str): Name of project
#     '''
    
#     provenance = pinery + '/samples'
    
#     headers = {'accept': 'application/json',}
#     params = {'project': project,}
#     response = requests.get(provenance, params=params, headers=headers)

#     cases = []
#     if response.ok:
#         L = response.json()
#         cases = list(set(map(lambda x: '_'.join(x.split('_')[:2]), [i['name'] for i in L])))     
                
#     return cases   
    

# def get_parent_sample_info(pinery, project, samples):
#     '''
#     (str, str, dict) -> dict
    
#     Return a dictionary with sample information, including samples not sequenced
#     for each project
    
#     Parameters
#     ----------
#     - pinery (str): Pinery API
#     - project (str): Name of project
#     - samples (dict): Dictionary with samples information extracted from Pinery
#     '''
    
#     cases = get_sample_info(pinery, project)
#     L = [samples[i] for i in cases if i in samples]
        
#     return L   












# def remove_table(database, table):
#     '''
#     (str, str) -> None
    
#     Drop table in database
    
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - table (str): Table of interest
#     '''
    
#     conn = sqlite3.connect(database)
#     cur = conn.cursor()
#     cur.execute("DROP TABLE {0}".format(table))
#     conn.commit()
#     conn.close()





# def add_project_info_to_db(database, pinery, project, lims, table = 'Projects'):
#     '''
#     (str, str, str, dict, str) -> None
    
#     Add project information into Projects table of database
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - pinery (str): Pinery API, http://pinery.gsi.oicr.on.ca
#     - project (str): Name of project of interest
#     - lims (dict): Dictionary with lims information extracted from pinery
#     - table (str): Name of Table in database. Default is Projects
#     '''
    
#     # remove project info from table
#     conn = sqlite3.connect(database)
#     cur = conn.cursor()
#     cur.execute('DELETE FROM {0} WHERE project_id = \"{1}\"'.format(table, project))
#     conn.commit()
#     conn.close()

#     # add project info in database
    
#     # get project info
#     projects = extract_project_info(pinery)
#     projects = {project: projects[project]}
    
#     # get column names
#     column_names = define_column_names()[table]

#     # add time stamp when project data was updated
#     projects[project]['last_updated'] = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
    
#     # add number of expected cases for project
#     samples = get_sample_info(pinery, project)
#     projects[project]['expected_samples'] = len(set(samples))
    
#     # add number of sequenced cases for project
#     projects[project]['sequenced_samples'] = len(lims[project].keys())
    
#     # add library types
#     library_type = []
#     for i in lims[project]:
#         for j in lims[project][i]:
#             library_type.append(lims[project][i][j]['library_type'])
#     library_type = ','.join(sorted(list(set(library_type))))        
#     projects[project]['library_types'] = library_type
        
#     # connect to db
#     conn = sqlite3.connect(database,  timeout=30)
#     cur = conn.cursor()
    
#     # order values according to column names
#     vals = [projects[project][i] for i in column_names]
#     # insert project info
#     cur.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), tuple(vals)))
#     conn.commit()
    
#     conn.close()


# def parse_calculate_contamination_db(calcontaqc_db):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with information collected from the calculate contamination
#     QC-etl sqlite cache
    
#     Parameters
#     ----------
#     - calcontaqc_db (str): Path to the calculate contamination database
#     '''

#     # get tables in db
#     conn = sqlite3.connect(calcontaqc_db)
#     cur = conn.cursor()
#     cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
#     tables = cur.fetchall()
#     tables = [i[0] for i in tables]    
#     conn.close()
#     # connect to db and retrieve all data
#     conn = connect_to_db(calcontaqc_db)
#     data = conn.execute('SELECT * FROM {0}'.format(tables[0])).fetchall()  
#     conn.close()

#     D = {}
    
#     for i in data:
#         case_id = i['Donor']
#         group_id = i['Group ID']
#         library_type = i['Library Design']
#         tissue_origin = i['Tissue Origin']
#         tissue_type = i['Tissue Type']
#         contamination = i['contamination']
#         merged_limskey = i['Merged Pinery Lims ID']
#         merged_limskey = list(map(lambda x: x.strip(), merged_limskey.replace('[', '').replace(']', '').replace('\"', '').split(',')))
#         merged_limskey = ';'.join(sorted(merged_limskey))
#         sample_id = '_'.join([i['Donor'], i['Tissue Type'], i['Tissue Origin'],
#                               i['Library Design'], i['Group ID']])
        
#         if case_id not in D:
#             D[case_id] = {}
#         if sample_id not in D[case_id]:
#             D[case_id][sample_id] = []
#         d = {'case_id': case_id,
#              'group_id': group_id,
#              'library_type': library_type, 
#              'tissue_origin': tissue_origin,
#              'tissue_type': tissue_type,
#              'contamination': contamination,
#              'merged_limskey': merged_limskey,
#              'sample_id': sample_id
#              }
#         D[case_id][sample_id].append(d)
           
#     return D                         
                         

# def add_contamination_info(database, calcontaqc_db, donors_to_update, table = 'Calculate_Contamination'):
#     '''
#     (str, str, dict, str) -> None
    
#     Parse the calcontaqc_db, reformat data and add information to the Calculate_Contamination
#     table of the database
    
#     Parameters
#     ----------
#     - database (str): Path to the sqlite database
#     - calcontaqc_db (db): Path to the calculate contamination database
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): name of the table in the database. Default is Calculate_Contamination 
#     '''

#     # remove rows for donors to update
#     if donors_to_update:
#         delete_records(donors_to_update, database, table)
#         print('deleted records in Calculate_Contamination')

#         # collect ata from the calculate contamination cache
#         D = parse_calculate_contamination_db(calcontaqc_db)
        
#         if D:
#             # make a list of data to insert
#             newdata = []
#             # connect to db
#             conn = sqlite3.connect(database)
#             # get column names
#             column_names = define_column_names()[table]

#             # add data
#             for case_id in D:
#                 if case_id in donors_to_update and donors_to_update[case_id] != 'delete':
#                     for sample_id in D[case_id]:
#                         for i in D[case_id][sample_id]:
#                            L = [i[j] for j in  column_names]
#                            newdata.append(L)
        
#         if newdata:
#             vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#             conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#             conn.commit()
#             conn.close()
        


# def count_files(project_name, database, file_table = 'Files'):
#     '''
#     (str, str, str) -> dict
    
#     Returns a dictionary with the number of files for each workflow in project
    
#     Parameters
#     ----------
#     - project_name (str): Name of project of interest
#     - database (str): Path to the sqlite database
#     - file_table (str): Name of the table with File information
#     '''
    
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT DISTINCT {0}.file, {0}.wfrun_id FROM {0} WHERE {0}.project_id = '{1}'".format(file_table, project_name)).fetchall()
#     conn.close()
    
#     counts = {}
#     for i in data:
#         counts[i['wfrun_id']] = counts.get(i['wfrun_id'], 0) + 1
       
#     return counts


# def update_workflow_information(project, workflows, D, key):
#     '''
#     (str, dict, dict, str) -> None
    
#     Update in place the dictionary workflows for a given project using key to store information contains in dictionary D
    
#     Parameters
#     ----------
#     - project (str): Project of interest
#     - workflows (dict): Dictionary containing workflow information for a given project
#     - D (dictionary): Dictionary containing file count or amount of input sequencing lanes for each workflow of project
#     - key (str): Key use to store the information from D in workflows.
#                  Valid keys are file_count or lane_count
#     '''
    
#     for workflow_id in D:
#         # workflow id may not be in workflows when updating 
#         if workflow_id in workflows[project]:
#             workflows[project][workflow_id][key] = D[workflow_id]
    

            
# def add_fileQC_info_to_db(database, project, nabu_api, matched_ids, donors_to_update, table='FilesQC'):
#     '''
#     (str, str, str, dict, str, str) -> None
    
#     Inserts file QC information in database's FilesQC table
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - project (str): Name of project of interest
#     - nabu_api (str): URL of the nabu API
#     - matched_ids (dict): Dictionary of matched file swids and donor ids for each project in FPR
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): Table in database storing the file QC. Default is FilesQC
#     '''

#     # remove rows for donors to update
#     if donors_to_update:
#         delete_records(donors_to_update, database, table)
#         print('deleted records in FilesQC')

#         # collect QC info from nabu
#         D = collect_qc_info(project, database, nabu_api)
    
#         # check that data is recorded in nabu for project
#         if D:
#             # make a list of data to insert
#             newdata = []
#             # connect to db
#             conn = sqlite3.connect(database)
#             # get column names
#             column_names = define_column_names()[table]

#             # add data
#             for file_swid in D[project]:
#                 # check that file swid is recorded in FPR for the same project
#                 if file_swid in matched_ids[project]:
#                     if matched_ids[project][file_swid] in donors_to_update and donors_to_update[matched_ids[project][file_swid]] != 'delete':
#                         L = [D[project][file_swid][i] for i in column_names if i in D[project][file_swid]]
#                         L.insert(0, matched_ids[project][file_swid])
#                         L.insert(0, project)
#                         L.insert(0, file_swid)
#                         newdata.append(L)
#             vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#             conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#             conn.commit()
#             conn.close()
            



# def add_samples_info_to_db(database, project, pinery, table, donors_to_update, sample_info):
#     '''
#     (str, str, str, dict, dict, dict) -> None
    
#     Inserts samples data into Samples table of database    
    
#     Parameters
#     ----------
#     - database (str): Path to the databae file
#     - project (str): Name of project of interest
#     - pinery (str): Pinery API
#     - table (str): Name of table in database
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - sample_info (dict): Dictionary with sample information extracted from Pinery
#     '''
    
    
#     # remove rows for donors to update
#     if donors_to_update:
#         delete_records(donors_to_update, database, table)
#         print('deleted records in Samples')

#         # collect information about samples
#         samples = get_parent_sample_info(pinery, project, sample_info)
    
#         if samples:
#             # make a list of row data
#             newdata = []
#             # connect to db
#             conn = sqlite3.connect(database)
             
#             # get column names
#             data = conn.execute("SELECT * FROM {0} WHERE project_id = '{1}';".format(table, project))
#             column_names = [column[0] for column in data.description]

#             # add data into table
#             for i in samples:
#                 if i['case'] in donors_to_update and donors_to_update[i['case']] != 'delete':
#                    L = [i['case'], i['donor_id'], i['species'], i['sex'], i['miso'],
#                          i['created_date'], i['modified_date'], project, i['project']]          
#                    newdata.append(L)
            
#             vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#             conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#             conn.commit()
#             conn.close()


# def add_blocks_to_db(database, project, expected_workflows, table, pipeline, donors_to_update):
#     '''
#     (str, str, str, str) -> None
    
#     Inserts WGS blocks data into WGS_blocks table of database    
    
#     Parameters
#     ----------
#     - database (str): Path to the databae file
#     - project (str): Name of project of interest
#     - expected_workflows (list): List of expected workflow names to define a complete block
#     - table (str): Name of table in database
#     - pipeline (str): Indicate WT or WGS pipeline
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     '''
    
#     # remove rows for donors to update
#     if donors_to_update:
#         delete_records(donors_to_update, database, table)
    
#         if pipeline == 'WGS':
#             # get the WGS blocks for donors in project
#             blocks = find_WGS_blocks(project, database, expected_workflows, donors_to_update)
#         elif pipeline == 'WT':
#             blocks = find_WT_blocks(project, database, expected_workflows, donors_to_update)

#         if blocks:
#             # make a list of row data
#             newdata = []
            
#             # connect to db
#             conn = sqlite3.connect(database)
#             # get column names
#             data = conn.execute("SELECT * FROM {0} WHERE project_id = '{1}';".format(table, project))
#             column_names = [column[0] for column in data.description]

#             # add data into table
#             for d in blocks:
#                 # loop over samples and blocks
#                 for samples in d:
#                     for block in d[samples]:
#                         if d[samples][block]['case_id'] in donors_to_update and donors_to_update[d[samples][block]['case_id']] != 'delete':
#                             L = [d[samples][block]['project_id'],
#                                  d[samples][block]['case_id'],
#                                  samples,
#                                  block,
#                                  ';'.join(d[samples][block]['workflows']),
#                                  d[samples][block]['name'],
#                                  d[samples][block]['date'],
#                                  d[samples][block]['release'],
#                                  d[samples][block]['complete'],
#                                  d[samples][block]['clean'],
#                                  d[samples][block]['network']]
#                             newdata.append(L)
#             vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#             conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#             conn.commit()
#             conn.close()





# def collect_lims_info(pinery):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with lims information extracted from sample-provenance in Pinery
        
#     Parameters
#     ----------
#     - pinery (str): Pinery API
#     ''' 

#     sample_provenance = pinery + '/provenance/v9/sample-provenance'
#     # get sample info from pinery
#     L = get_provenance_data(sample_provenance)
    
#     # store lims information for each library of each project
#     # {project: {sample: {library_info}}}    
#     D = {}

#     for i in L:
#         project = i['studyTitle']
#         sample = i['rootSampleName']  
#         if project not in D:
#             D[project] = {}
#         if sample not in D[project]:
#             D[project][sample] = {}
    
#         # collect sample information
#         d = collect_info(i['sampleAttributes'], ['geo_tissue_type', 'geo_external_name', 'geo_tissue_origin',
#                  'geo_library_source_template_type', 'geo_prep_kit',
#                  'geo_tissue_preparation', 'geo_receive_date', 'geo_group_id', 'geo_group_id_description'], ['tissue_type', 'ext_id', 'tissue_origin', 'library_type', 
#                          'prep', 'tissue_prep', 'sample_received_date', 'group_id', 'group_id_description'])
#         # add library name
#         library = i['sampleName']
                        
#         # store sample information
#         if library not in D[project][sample]:
#             D[project][sample][library] = d
            
#         #update sample information for each library if some fields are missing
#         for k in d:
#             if D[project][sample][library][k] == '':
#                 D[project][sample][library][k] = d[k]    
    
#     return D




# def collect_parent_sample_info(pinery):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with sample information extracted from sample-provenance in Pinery
        
#     Parameters
#     ----------
#     - pinery (str): Pinery API
#     ''' 

#     provenance = pinery + '/samples'
#     # get sample info from pinery
#     L = get_provenance_data(provenance)
    
#     D = {}

#     for i in L:
#         if i['name'].count('_') == 1:
#             case = i['name']
#             sex, species, donor_id = '', '', ''
#             for j in i['attributes']:
#                 if j['name'] == 'Sex':
#                     sex = j['value']
#                 if j['name'] == 'Organism':
#                     species= j['value']
#                 if j['name'] == 'External Name':
#                     donor_id = j['value']
#             sample_id = i['id']
#             miso = 'https://miso.oicr.on.ca/miso/sample/{0}'.format(sample_id.replace('SAM', ''))    
#             project = i['project_name']
#             modified_date = i['modified_date']
#             created_date = i['created_date']
#             D[case] = {'project': project, 'sex': sex, 'species': species,
#                        'miso': miso, 'donor_id': donor_id, 'case': case, 'modified_date': modified_date,
#                        'created_date': created_date}
        
#     return D


# def generate_database(args):
#     '''
#     (str, str, str, str, str, str)
    
#     Pull data from resources to generate the waterzooi database
        
#     Parameters:
#     -----------
#     - fpr (str): Path to Path to the File Provenance Report
#     - nabu (str): URL of the Nabu API
#     - pinery (str): Pinery API
#     - database (str): Path to the database file
#     - lims_info (str): Path to the json file with lims information
#     - samples_info (str): Path to the json file with sample information
#     '''

#     # collect lims info
#     lims_info = collect_lims_info(args.pinery)
#     print('collected lins info')
#     # collect sample info
#     samples_info = collect_parent_sample_info(args.pinery)
#     print('samples info')
           
#     # collect file info from FPR
#     # get the list of valid projects
#     # extract project information from project provenance
#     projects = extract_project_info(args.pinery)
#     # filter out completed projects
#     projects = filter_completed_projects(projects)
#     # make a list of projects with data in Prinery and FPR
#     projects = get_valid_projects(projects, args.fpr)
#     projects.sort()
    
    
    
    
    
#     fpr_data = collect_file_info_from_fpr(args.fpr, projects)
    
#     # create database if file doesn't exist
#     if os.path.isfile(args.database) == False:
#         initiate_db(args.database)
    
#     # collect workflow inputlibraries
#     wf_input = extract_workflow_info(args.fpr, projects)
#     # match file swids to donor ids
#     matched_ids = match_donor_to_file_swid(fpr_data, projects)
#     # get the parent file ids of each workflow
#     parents = get_parent_file_ids(args.fpr, projects)
    
#     # loop over projects
#     for project in projects:
#         print(project)
#         try:
#             # get the donors for which records need to be updated
#             start = time.time()
#             donors_to_update = donors_info_to_update(args.database, fpr_data, project, 'Checksums')
#             print('donors to update:', len(donors_to_update))
#             if donors_to_update:
#                 


                  
#                 
#                 # add sample information
#                 add_samples_info_to_db(args.database, project, args.pinery, 'Samples', donors_to_update, samples_info)
#                 print('added samples') 
#                 

                  

#                 # add calculate contamination info
#                 add_contamination_info(args.database, args.calcontaqc_db, donors_to_update, table = 'Calculate_Contamination')
#                 print('added call-ready contamination')
#                 # add file QC info
#                 add_fileQC_info_to_db(args.database, project, args.nabu, matched_ids, donors_to_update, 'FilesQC')
#                 print('added filesqc')
#                 
                  
#                 print('added workflow relationships')
#                 # add WGS blocks
#                 expected_WGS_workflows = sorted(['mutect2', 'variantEffectPredictor', 'delly', 'varscan', 'sequenza', 'mavis']) 
#                 add_blocks_to_db(args.database, project, expected_WGS_workflows, 'WGS_blocks', 'WGS', donors_to_update)
#                 print('added wgs blocks')  
#                 # add WT blocks
#                 expected_WT_workflows = sorted(['arriba', 'rsem', 'starfusion', 'mavis'])
#                 add_blocks_to_db(args.database, project, expected_WT_workflows, 'WT_blocks', 'WT', donors_to_update)
#                 print('added WT blocks')  
#                 
                  
#         except:
#             print('could not add data for {0}'.format(project))
#             print(print(traceback.format_exc()))
#             print('----')
    
    
    
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(prog = 'waterzooiDataCollector.py', description='Script to add data to the waterzooi database')
#     parser.add_argument('-f', '--fpr', dest='fpr', default = '/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to the File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
#     parser.add_argument('-n', '--nabu', dest='nabu', default='https://nabu-prod.gsi.oicr.on.ca', help='URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca')
#     parser.add_argument('-p', '--pinery', dest = 'pinery', default = 'http://pinery.gsi.oicr.on.ca', help = 'Pinery API. Default is http://pinery.gsi.oicr.on.ca')
#     parser.add_argument('-cq', '--calcontaqc', dest = 'calcontaqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest', help = 'Path to the merged rnaseq calculateContamination database. Default is /scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest')
#     parser.add_argument('-db', '--database', dest='database', help='Path to the waterzooi database', required=True)    
#     parser.set_defaults(func=generate_database)
 
#     # get arguments from the command line
#     args = parser.parse_args()
#     args.func(args)
    

