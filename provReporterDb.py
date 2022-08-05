# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 14:30:23 2022

@author: rjovelin
"""


import sqlite3
import json
import requests
import gzip
import os
import argparse
import time


def extract_project_info(project_provenance = 'http://pinery.gsi.oicr.on.ca/sample/projects'):
    '''
    (str) -> dict
    
    Returns a dictionary with project information pulled down from the
    project_provenance Pinary API

    Parameters
    ----------
    - project_provenance (str): Pinery API, http://pinery.gsi.oicr.on.ca/sample/projects
    '''
    
    response = requests.get(project_provenance)
    if response.ok:
        L = response.json()
    else:
        L = []
    
    D = {}
    
    if L:
        for i in L:
            name = i['name']
            assert name not in D
            D[name] = {'project_id': name}
            for j in ['pipeline', 'description', 'active', 'contact_name', 'contact_email']:
                if j in i:
                    D[name][j] = i[j] 
                else:
                    D[name][j] = ''
                if j == 'active' and j in i:
                    if i[j]:
                        D[name][j] = 'Active'
                    else:
                        D[name][j] = 'Completed'
    return D                            


def collect_info(data, names, keys):
    '''
    (dict, list, list) -> dict

    Returns a dictionary with keys by extracting specific information using names from data
    Note, the value of any name in data may be a string or single-element-list

    Parameters
    ----------
    - data (dict): Dictionary with data to be extracted
    - names (list): List of keys of interest in data
    - keys (list): List of renamed keys in output dictionary
    '''
    
    d = {}
    for j in range(len(names)):
        if names[j] in data:
            if type(data[names[j]]) == list:
                d[keys[j]] = data[names[j]][0]         
            else:
                d[keys[j]] = data[names[j]]
        else:
            d[keys[j]] = ''

    return d



def extract_qc_status_from_nabu(project, database, file_table = 'Files', nabu_api = 'http://gsi-dcc.oicr.on.ca:3000'):
    '''
    (str, str, str, str) -> dict

    Returns a dictionary with QC status of all files in project    
        
    Parameters
    ----------
    - project (str): Name of project of interest
    - database (str): Path to the database file
    - file_table (str): Table in database storing file information. Default is Files
    - nabu_api (str): URL of the nabu API
    '''
        
    response = requests.get(nabu_api + '/fileqcs?project={0}'.format(project))
    
    if response.ok:
        L = response.json()['fileqcs']
    else:
        L = []
    
    D = {project: {}}
    
    if L:
        for i in L:
            qc = collect_info(i, ['skip', 'user', 'date', 'qcstatus', 'ref', 'stalestatus'], ['skip', 'user', 'date', 'status', 'reference', 'fresh']) 
            file_swid = int(i['fileswid'])
            D[project][file_swid] = qc 
    
    # QC may not be available for all files
    # retrieve file swids from database instead of prasing again fpr
    # add empty values to qc fields if qc not available
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    cur.execute('SELECT {0}.file_swid FROM {0} WHERE {0}.project_id = \"{1}\"'.format(file_table, project))
    records = cur.fetchall()
    
    if records:
        print('extract_qc_status', records[0])
    
    conn.close()
    
    for project in D:
        for file_swid in records:
            if file_swid not in D[project]:
                D[project][file_swid] = {'skip': '', 'user': '', 'date': '', 'status': '', 'ref': '', 'fresh': ''}
    return D


def collect_qc_info(database, file_table = 'Files', project_provenance = 'http://pinery.gsi.oicr.on.ca/sample/projects', nabu_api = 'http://gsi-dcc.oicr.on.ca:3000'):
    '''
    (str, str, str, str) -> dict
    
    Returns a dictionary with QC information for all files in FPR and for all projects defined in Pinery
             
    Parameters
    ----------
    - database (str): Path to the database file
    - file_table (str): Table in database storing file information. Default is Files 
    - project_provenance (str): Pinery API
    - nabu_api (str): URL of the nabu API
    '''
        
    # make a list of projects defined in Pinery
    projects = list(extract_project_info(project_provenance).keys())
        
    # collect QC information for all files in each project
    D = {}
    for project in projects:
        qc = extract_qc_status_from_nabu(project, database, file_table, nabu_api)
        assert project not in D
        D[project] = qc[project]
    
    return D


def collect_file_info_from_fpr(fpr):
    '''
    (str) -> dict

    Returns a dictionary with file information extracted from File Provenance Report
        
    Parameters
    ----------
    - fpr (str): Path to File Provenance Report file
    '''

    infile = open_fpr(fpr)
    
    D = {}
    
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            # get project and initiate dict
            project = line[1]
            if project not in D:
                D[project] = {}
            # get file path
            file = line[46]
            # get md5sums
            md5sum = line[47]
            # get file_swid
            file_swid = int(line[44])
            # get the library type
            d = collect_info({k.split('=')[0]:k.split('=')[1] for k in line[17].split(';')},
                         ['geo_library_source_template_type'], ['library_type'])
            
            # get file attributes
            attributes = line[45]
            if attributes:
                attributes = attributes.split(';')
                attributes = json.dumps({k.split('=')[0]: k.split('=')[1] for k in attributes})
            else:
                attributes = ''
            
            # get workflow
            workflow = line[30]
            # get workflow version
            version = line[31]        
            # get workflow run accession
            workflow_run = int(line[36])
            # collect file info
            D[project][file_swid] = {'md5sum': md5sum, 'workflow': workflow, 'version': version,
                   'wfrun_id': workflow_run, 'file': file, 'library_type': d['library_type'], 'attributes': attributes}
    infile.close()
    return D



def is_gzipped(file):
    '''
    (str) -> bool

    Return True if file is gzipped

    Parameters
    ----------
    - file (str): Path to file
    '''
    
    # open file in rb mode
    infile = open(file, 'rb')
    header = infile.readline()
    infile.close()
    if header.startswith(b'\x1f\x8b\x08'):
        return True
    else:
        return False


def open_fpr(fpr):
    '''
    (str) -> _io.TextIOWrapper
    
    Returns a file open for reading
    
    Parameters
    ----------
    - fpr (str): Path to File Provenance Report file
    '''

    # open provenance for reading. allow gzipped file or not
    if is_gzipped(fpr):
        infile = gzip.open(fpr, 'rt', errors='ignore')
    else:
        infile = open(fpr)
    return infile


def get_workflow_relationships(fpr):
    '''
    (str, str) -> (dict, dict, dict)

    Returns a tuple with dictionaries with worklow information, parent file ids for each workflow,
    and workflow id for each file and project of interest
    
    Parameters
    ----------
    - fpr (str): Path to the File Provenance Report
    '''

    F, P, W = {}, {}, {}


    start = time.time()
    

    # open fpr
    infile = open_fpr(fpr)
    # skip header
    infile.readline()
    
    
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            # get project name
            project = line[1]
            # get workflow, workflow version and workflow run accession
            workflow, workflow_version, workflow_run = line[30], line[31], int(line[36])
            # get file swid
            file_swid = int(line[44])
            input_files = line[38]
            if input_files:
                input_files = sorted(input_files.split(';'))
            else:
                input_files = []
                   
            # get workflow attributes
            attributes = line[37]
            if attributes:
               attributes = attributes.split(';')
               attributes = json.dumps({k.split('=')[0]: k.split('=')[1] for k in attributes if k.split('=')[0] not in ['cromwell-workflow-id', 'major_olive_version']})
            else:
               attributes = ''                
       
            if project not in P:
                P[project] = {}
      
            if workflow_run in P[project]:
                assert P[project][workflow_run] == input_files
            else:
                P[project][workflow_run] = input_files
                
            if project not in W:
                W[project] = {}
            if workflow_run in W[project]:
                assert W[project][workflow_run] == {'wfrun_id': workflow_run, 'wfv': workflow_version, 'wf': workflow, 'attributes': attributes}
            else:
                W[project][workflow_run] = {'wfrun_id': workflow_run, 'wfv': workflow_version, 'wf': workflow, 'attributes': attributes}
        
            if project not in F:
                F[project] = {}
            if file_swid in F[project]:
                assert F[project][file_swid] == workflow_run
            else:
                F[project][file_swid] = workflow_run
    
    infile.close()
    
    
    end = time.time()
    print('parsing workflow relationships', end - start)

    
    return W, P, F        



def identify_parent_children_workflows(P, F):
    '''
    (dict, dict, dict) -> (dict, dict)     
    
    Returns a tuple with dictionaries of children: parents workflows and parents: children workflow  
    relationsips for a given project
        
    Parameters
    ----------
    - P (dict): Input file ids for each workflow run id
    - F (dict): Map of file id and workflow id
    '''
    
    start = time.time()
    
    
    # parents record parent-child workflow relationships
    # children record child-parent workflow relationships
    parents, children = {}, {}
       
    for project in P:
        if project not in children:
            children[project] = {}
        for workflow in P[project]:
            parent_workflows = sorted(list(set([F[project][i] for i in P[project][workflow] if i in F[project]])))
            children[project][workflow] = parent_workflows
    
    for project in children:
        if project not in parents:
            parents[project] = {}
        for workflow in children[project]:
            for parent_workflow in children[project][workflow]:
                if parent_workflow not in parents[project]:
                    parents[project][parent_workflow] = [workflow]
                else:
                    parents[project][parent_workflow].append(workflow)
    
    for project in parents:
        for workflow in parents[project]:
            parents[project][workflow] = sorted(list(set(parents[project][workflow])))
    
    end = time.time()
    
    print('identify parent_children', end - start)
    
    return parents, children




def get_provenance_data(provenance):
    '''
    (str) -> list
    
    Returns a list of dictionary with lims information for each library
    
    Parameters
    ----------
    - provenance (str): URL of the pinery provenance API 
    '''
    
    response = requests.get(provenance)
    if response.ok:
        L = response.json()
    else:
        L = []
    
    return L



def collect_lims_info(sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
    '''
    (str) -> dict
    
    Returns a dictionary with library information for each sample of each project

    Parameters
    ----------
    - sample_provenance (str): URL of the pinery sample_provenance API
    ''' 

    # get sample info from pinery
    L = get_provenance_data(sample_provenance)
    
    # store lims information for each library of each project
    # {project: {sample: {library_info}}}    
    D = {}

    for i in L:
        project = i['studyTitle']
        sample = i['rootSampleName']  
        if project not in D:
            D[project] = {}
        if sample not in D[project]:
            D[project][sample] = {}
    
        # collect sample information
        d = collect_info(i['sampleAttributes'], ['geo_tissue_type', 'geo_external_name', 'geo_tissue_origin',
                 'geo_library_source_template_type', 'geo_prep_kit',
                 'geo_tissue_preparation', 'geo_receive_date', 'geo_group_id', 'geo_group_id_description'], ['tissue_type', 'ext_id', 'tissue_origin', 'library_type', 
                         'prep', 'tissue_prep', 'sample_received_date', 'group_id', 'group_id_description'])
        # add library name
        library = i['sampleName']
                        
        # store sample information
        if library not in D[project][sample]:
            D[project][sample][library] = d
            
        #update sample information for each library if some fields are missing
        for k in d:
            if D[project][sample][library][k] == '':
                D[project][sample][library][k] = d[k]    
    return D    




def extract_workflow_info(fpr):
    '''
    (str) -> dict
    
    Returns a dictionary with library input information for all workflows for each project
    
    Parameters
    ----------
    - fpr (str): Path to the File Provenance Report
    '''

    # create a dict to store workflow input libraries
    D = {}

    infile = open_fpr(fpr)
    for line in infile:
        line = line.rstrip()
        if line:
            line = line.split('\t')
            # get project name
            project = line[1]
            if project not in D:
                D[project] = {}
            
            # get sample name
            sample = line[7]
            # get workflow name and workflow run accession
            workflow, workflow_run = line[30], int(line[36])
                                  
            # get lane and run
            run, lane = line[18], line[24]            
            # get library and limskey
            library, limskey  = line[13], line[56]
            
            # get barcode and platform
            barcode, platform = line[27], line[22]
            
            d = {'run': run, 'lane': lane, 'library': library, 'limskey': limskey, 'barcode': barcode, 'platform': platform}
            
            if project not in D:
                D[project] = {}
            if workflow_run not in D[project]:
                D[project][workflow_run] = {'sample': sample, 'workflow': workflow, 'libraries': [d]}
            else:
                assert sample == D[project][workflow_run]['sample']
                assert workflow == D[project][workflow_run]['workflow']
                if d not in D[project][workflow_run]['libraries']:
                    D[project][workflow_run]['libraries'].append(d)
    return D            



def define_column_names():
    '''
    (None) -> dict

    Returns a dictionary with column names for each table in database
    '''

    # create dict to store column names for each table {table: [column names]}
    column_names = {'Workflows': ['wfrun_id', 'wf', 'wfv', 'project_id', 'attributes'],
                    'Parents': ['children_id', 'parents_id', 'project_id'],
                    'Children': ['parents_id', 'children_id', 'project_id'],
                    'Projects': ['project_id', 'pipeline', 'description', 'active', 'contact_name', 'contact_email'],
                    'Files': ['file_swid', 'project_id', 'md5sum', 'workflow', 'version', 'wfrun_id', 'file', 'library_type', 'attributes'],
                    'FilesQC': ['file_swid', 'project_id', 'skip', 'user', 'date', 'status', 'ref', 'fresh'],
                    'Libraries': ['library', 'sample', 'tissue_type', 'ext_id', 'tissue_origin',
                                  'library_type', 'prep', 'tissue_prep', 'sample_received_date', 'group_id', 'group_id_description', 'project_id'],
                    'Workflow_Inputs': ['library', 'run', 'lane', 'wfrun_id', 'limskey', 'barcode', 'platform', 'project_id']}
        
    return column_names


def define_column_types():
    '''
    (None) -> dict

    Returns a dictionary with column types for each table in database
    '''

    # create dict to store column names for each table {table: [column names]}
    column_types = {'Workflows': ['INTEGER', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'TEXT'],
                    'Parents': ['INTEGER', 'INTEGER', 'VARCHAR(128)'],
                    'Children': ['INTEGER', 'INTEGER', 'VARCHAR(128)'],
                    'Projects': ['VARCHAR(128) PRIMARY KEY NOT NULL UNIQUE', 'VARCHAR(128)',
                                  'TEXT', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(256)'],
                    'Files': ['INTEGER PRIMARY KEY NOT NULL UNIQUE', 'VARCHAR(128)',
                              'VARCHAR(256)', 'VARCHAR(128)', 'VARCHAR(128)',
                              'INTEGER', 'TEXT', 'VARCHAR(128)', 'TEXT'],
                    'FilesQC': ['INTEGER PRIMARY KEY NOT NULL UNIQUE', 'VARCHAR(128)',
                                'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)'],
                    'Libraries': ['VARCHAR(256) PRIMARY KEY NOT NULL', 'VARCHAR(128)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 
                                  'VARCHAR(256)', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(128)'],
                    'Workflow_Inputs': ['VARCHAR(128)', 'VARCHAR(256)', 'INTEGER', 'INTEGER', 
                                        'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)']}
    
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


    if table  in ['Workflows', 'Parents', 'Children', 'Files', 'FilesQC', 'Libraries', 'Workflow_Inputs']:
        constraints = '''FOREIGN KEY (project_id)
            REFERENCES Projects (project_id)
            ON DELETE CASCADE ON UPDATE CASCADE'''
        table_format = table_format + ', ' + constraints 
    
    if table in ['Parents', 'Children']:
        constraints = '''FOREIGN KEY (parents_id)
          REFERENCES Workflows (wfrun_id)
          ON DELETE CASCADE ON UPDATE CASCADE,
          FOREIGN KEY (children_id)
              REFERENCES Workflows (wfrun_id)
              ON DELETE CASCADE ON UPDATE CASCADE''' 
        table_format = table_format + ', ' + constraints
    
    if table == 'Parents':
        table_format = table_format + ', PRIMARY KEY (children_id, parents_id, project_id)'
    
    if table == 'Children':
        table_format = table_format + ', PRIMARY KEY (parents_id, children_id, project_id)'
    
    if table == 'Worklows':
        table_format = table_format + ', PRIMARY KEY (wfrun_id, project_id)'
        
    if table == 'Files':
        constraints = '''FOREIGN KEY (wfrun_id)
            REFERENCES Workflows (wfrun_id)
            ON DELETE CASCADE ON UPDATE CASCADE,
            FOREIGN KEY (file_swid)
               REFERENCES FilesQC (file_swid)
               ON DELETE CASCADE ON UPDATE CASCADE'''
        table_format = table_format + ', ' + constraints
    
    if table == 'Workflow_Inputs':
        constraints = '''FOREIGN KEY (wfrun_id)
            REFERENCES Workflows (wfrun_id)
            ON DELETE CASCADE ON UPDATE CASCADE,
            FOREIGN KEY (library)
              REFERENCES Libraries (library)
              ON DELETE CASCADE ON UPDATE CASCADE'''
        table_format = table_format + ', ' + constraints
    
    # connect to database
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    #cur.execute('DROP TABLE IF EXISTS {0};'.format(table))
    #conn.commit()
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
    for i in ['Projects', 'Workflows', 'Parents', 'Children', 'Files', 'FilesQC', 'Libraries', 'Workflow_Inputs']:
        if i not in tables:
            create_table(database, i)



def add_project_info_to_db(database, table = 'Projects', project_provenance = 'http://pinery.gsi.oicr.on.ca/sample/projects'):
    '''
    (str, str) -> None
    
    Add project information into Projects table of database
       
    Parameters
    ----------
    - database (str): Path to the database file
    - project_provenance (str): Pinery API, http://pinery.gsi.oicr.on.ca/sample/projects
    '''

    # connect to db
    conn = sqlite3.connect(database)
    cur = conn.cursor()
            
    # get project info
    projects = extract_project_info(project_provenance)
    
    to_remove = []
    for i in projects:
        if i not in ['HCCCFD', 'TGL01MOH', 'KLCS', 'BARON', 'SIMONE', 'HLCS', 'ARCH1']:
            to_remove.append(i)
    for i in to_remove:
        del projects[i]
    
    
    
    
    
    
    # get existing records
    cur.execute('SELECT {0}.project_id FROM {0}'.format(table))
    records = cur.fetchall()
    records = [i[0] for i in records]

    # get column names
    column_names = define_column_names()[table]

    # add or update data
    for project in projects:
        # order values according to column names
        vals = [projects[project][i] for i in column_names]
        if project in records:
            # update project info
            for i in projects[project]:
                cur.execute('UPDATE {0} SET {1} = \"{2}\" WHERE project_id=\"{3}\"'.format(table, i, projects[project][i], project))  
                conn.commit()
        else:
            # insert project info
            cur.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), tuple(vals)))
            conn.commit()
    # remove any projects in database not anymore defined in pinery    
    for project in records:
        if project not in projects:
            cur.execute('DELETE FROM {0} WHERE project_id=\"{1}\"'.format(table, projects))
            conn.commit()
    conn.close()




def add_workflows(workflows, database, table = 'Workflows'):
    '''
    (dict, str, str) -> None
    
    Inserts or updates workflow information to table Workflows in database
           
    Parameters
    ----------
    - workflows (dict): Dictionary with workflow information
    - database (str): Path to the database file
    - workflow_table (str): Name of the table storing workflow information. Default is Workflows
    '''
       
    # connect to db
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    
    #get existing records
    cur.execute('SELECT {0}.wfrun_id, {0}.project_id FROM {0}'.format(table))
    records = cur.fetchall()
    
    # get column names
    column_names = define_column_names()[table]
    
    # make a list of fields to update
    to_update = [i for i in column_names]
    to_update.remove('project_id')
    to_update.remove('wfrun_id')
    
    for project in workflows:
        for workflow_run in workflows[project]:
            if (workflow_run, project) in records:
                # update workflows info
                for i in to_update:
                    cur.execute("UPDATE {0} SET {1} = '{2}' WHERE wfrun_id='{3}' AND project_id ='{4}'".format(table, i, workflows[project][workflow_run][i], workflow_run, project))
                    conn.commit()
            else:
                # insert data into table
                values = [workflows[project][workflow_run][i] for i in column_names if i in workflows[project][workflow_run]]
                values.insert(-1, project)
                cur.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), tuple(values)))
                conn.commit()
        
    # remove any workflows not defined anymore in FPR    
    for (workflow_run, project) in records:
        if project not in workflows:
            cur.execute('DELETE FROM {0} WHERE {0}.project_id=\"{1}\"'.format(table, project))
            conn.commit()
        elif workflow_run not in workflows[project]:
            cur.execute('DELETE FROM {0} WHERE wfrun_id = \"{1}\" AND project_id=\"{2}\"'.format(table, workflow_run, project))
            conn.commit()
    conn.close()


    
def add_workflow_relationships(D, database, table):    
    '''
    (dict, str, str)
    
    Inserts or updates parent-children workflow relatiionships to table Children in database
    
    Parameters
    ----------    
    - D (dict): Dictionary with children-parents (Parents table) or parent-children (Children table) workflow relationships
    - database (str): Path to the database file
    - table (str): Name of the table storing children-parents workflow relationships
    '''
    
    # connect to db
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    
    #get existing records
    if table == 'Parents':
        cur.execute('SELECT {0}.children_id, {0}.parents_id, {0}.project_id FROM {0}'.format(table))
    elif table == 'Children':
        cur.execute('SELECT {0}.parents_id, {0}.children_id, {0}.project_id FROM {0}'.format(table))
    records = cur.fetchall()
    
    # get column names
    column_names = define_column_names()[table]
    
    for project in D:
        for i in D[project]:
            for j in D[project][i]:
                if (i, j, project) not in records:
                    # insert data into table
                    cur.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), (i, j, project)))
                    conn.commit()
                            
    # remove any workflow relationships not defined anymore in FPR    
    for (i, j, project) in records:
        if project not in D:
            cur.execute('DELETE FROM {0} WHERE {0}.project_id=\"{1}\"'.format(table, project))
            conn.commit()
        elif i not in D[project]:
            if table == 'Parents':
                cmd1 = 'DELETE FROM {0} WHERE {0}.children_id = \"{1}\" AND {0}.project_id=\"{2}\"'.format(table, i, project)
            elif table == 'Children':
                cmd1 = 'DELETE FROM {0} WHERE {0}.parents_id = \"{1}\" AND {0}.project_id=\"{2}\"'.format(table, i, project)
            cur.execute(cmd1)
            conn.commit()
        elif j not in D[project][i]:
            if table == 'Parents':
                cmd2 = 'DELETE FROM {0} WHERE {0}.children_id = \"{1}\" AND {0}.parents_id = \"{2}\" AND {0}.project_id=\"{3}\"'.format(table, i, j, project)
            elif table == 'Children':
                cmd2 = 'DELETE FROM {0} WHERE {0}.parents_id = \"{1}\" AND {0}.children_id = \"{2}\" AND {0}.project_id=\"{3}\"'.format(table, i, j, project)
            cur.execute(cmd2)            
            conn.commit()
    conn.close()


def add_workflows_info_to_db(fpr, database, workflow_table = 'Workflows', parent_table = 'Parents', children_table = 'Children'):
    '''
    (str, str, str, str, str) -> None
    
    Inserts or updates workflow information and parent-children workflow relationships
 
    Parameters
    ----------    
    - fpr (str): Path to the File Provenance Report
    - database (str): Path to the database file
    - workflow_table (str): Name of the table storing workflow information. Default is Workflows
    - parent_table (str): Name of the table storing parents-children workflow relationships. Default is Parents
    - children_table (str): Name of the table storing children-parents workflow relationships. Default is Children
    '''

    # get the workflow inputs and workflow info
    workflows, parents, files = get_workflow_relationships(fpr)
    
    to_remove = []
    for i in workflows:
        if i not in ['HCCCFD', 'TGL01MOH', 'KLCS', 'BARON', 'SIMONE', 'HLCS', 'ARCH1']:
            to_remove.append(i)
    for i in to_remove:
        del workflows[i]
        del parents[i]
        del files[i]
    
    
    # identify parent-children workflow relationships
    parent_workflows, children_workflows = identify_parent_children_workflows(parents, files)

    # add workflow info
    add_workflows(workflows, database, workflow_table)
    print('added workflows')
    
    
    # add parents-children workflow relationships
    add_workflow_relationships(parent_workflows, database, parent_table)    
    print('added parents')
    
    # add children-parents workflow relationships    
    add_workflow_relationships(children_workflows, database, children_table)    
    print('added children')
   
    
    
def add_file_info_to_db(database, table, fpr, file_table = 'Files', project_provenance = 'http://pinery.gsi.oicr.on.ca/sample/projects', nabu_api = 'http://gsi-dcc.oicr.on.ca:3000'):
    '''
    (str, str, str, str, str, str) -> None
    
    Inserts or updates file QC information in database's FilesQC table
       
    Parameters
    ----------
    - database (str): Path to the database file
    - table (str): Table in database storing the QC or file information.
                   Accepted values are FilesQC or Files
    - fpr (str): Path to the File Provenance Report
    - file_table (str): Table in database storing file information. Default is Files 
    - project_provenance (str): Pinery API
    - nabu_api (str): URL of the nabu API
    '''

    assert table in ['Files', 'FilesQC']        

    # check if adding QC info or file info
    if table == 'FilesQC':
        # adding file QC information
        D = collect_qc_info(database, file_table, project_provenance, nabu_api)
    elif table == 'Files':
        # collect file info from FPR
        D = collect_file_info_from_fpr(fpr)
       
    # connect to db
    conn = sqlite3.connect(database)
    cur = conn.cursor()
        
    # get existing records
    cur.execute('SELECT {0}.file_swid FROM {0}'.format(table))
    records = cur.fetchall()
    records = [i[0] for i in records]
    
    if records:
        print('add_file_info', table, records[0])
    else:
        print('add_file_info', table, records)

    # get column names
    column_names = define_column_names()[table]

    # make a list of file_swids
    file_swids = []

    # add or update data
    for project in D:
        for file_swid in D[project]:
            file_swids.append(file_swid)
            L = [D[project][file_swid][i] for i in column_names if i in D[project][file_swid]]
            L.insert(0, project)
            L.insert(0, file_swid)
            if file_swid in records:
                # update QC info
                for i in range(1, len(column_names)):
                    cur.execute('UPDATE {0} SET {1} = \"{2}\" WHERE file_swid=\"{3}\"'.format(table, column_names[i], L[i], file_swid))  
                    conn.commit()
            else:
                # insert project info
                cur.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), tuple(L)))
                conn.commit()
    
    # remove any file in database not anymore defined in Nabu    
    for file_swid in records:
        if file_swid not in file_swids:
            cur.execute('DELETE FROM {0} WHERE {0}.file_swid=\"{1}\"'.format(table, file_swid))
            conn.commit()
    conn.close()



def add_library_info_to_db(database, table = 'Libraries', sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
    '''
    (str, str, str) -> None
    
    Inserts or updates library information in Libraries table of database    
    
    Parameters
    ----------
    - database (str): Path to the databae file
    - table (str): Table storing library in database. Default is Libraries
    - sample_provenance (str): URL of the sample provenance API
    '''
    
    # collect lims information
    lims = collect_lims_info(sample_provenance)
    
    # connect to db
    conn = sqlite3.connect(database)
    cur = conn.cursor()
       
    # get existing records
    cur.execute('SELECT {0}.library FROM {0}'.format(table))
    records = cur.fetchall()

    print('add_library_info', table, records[0])

    # get column names
    column_names = define_column_names()[table]

    # make a list of libraries
    libraries = []

    # add or update data
    for project in lims:
        for sample in lims[project]:
            for library in lims[project][sample]:
                libraries.append(library)
                L = [lims[project][library][i] for i in column_names if i in lims[project][library]]
                L.insert(0, sample)
                L.insert(0, library)
                L.append(project)
                if library in records:
                    # update QC info
                    for i in range(1, len(column_names)):
                        cur.execute('UPDATE {0} SET {0}.{1} = \"{2}\" WHERE {0}.library=\"{3}\"'.format(table, column_names[i], L[i], library))  
                        conn.commit()
                else:
                    # insert project info
                    cur.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), tuple(L)))
                    conn.commit()
   
    # remove any library in database not anymore defined in Pinery    
    for library in records:
       if library not in libraries:
           cur.execute('DELETE FROM {0} WHERE {0}.library=\"{1}\"'.format(table, library))
           conn.commit()
    conn.close()






def add_workflow_inputs_to_db(database, fpr, table = 'Workflow_Inputs'):
    '''
    (str, str, str) -> None
    
    Inserts or updates workflow input library information in table Workflow_Inputs of database    
    
    Parameters
    ----------
    - database (str): Path to the databae file
    - fpr (str): Path to the File Provenance Report
    - table (str): Table storing library in database. Default is Libraries
    '''

    # collect information about library inputs
    libraries = extract_workflow_info(fpr)
        
    # connect to db
    conn = sqlite3.connect(database)
    cur = conn.cursor()
       
    # get existing records
    cur.execute('SELECT {0}.library, {0}.run, {0}.lane FROM {0}'.format(table))
    records = cur.fetchall()

    print('add_wkf_input', table, records[0])


    # get column names
    column_names = define_column_names()[table]

    # make a list of (library, run, lane)
    # library, run and lane defines an input to workflow
    inputs = []

    # add or update data
    for project in libraries:
        for workflow_run in libraries[project]:
            for i in libraries[project][workflow_run]['libraries']:
                library, run, lane = i['library'], i['run'], i['lane']
                key = tuple(library, run, lane)
                inputs.append(key)
                L = [i[j] for j in column_names if j in i]
                L.append(project)
                L.insert(2, workflow_run)
                if key in records:
                    # update inputs info
                    for i in range(3, len(column_names)):
                        cur.execute('UPDATE {0} SET {0}.{1} = \"{2}\" WHERE {0}.library=\"{3}\" AND {0}.run=\"{4}\" AND {0}.lane=\"{5}\"'.format(table, column_names[i], L[i], library, run, lane))  
                        conn.commit()
                else:
                    # insert project info
                    cur.execute('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), tuple(L)))
                    conn.commit()
   
    # remove any inputs in database not anymore defined in FPR    
    for k in records:
       if library not in inputs:
           cur.execute('DELETE FROM {0} WHERE {0}.library=\"{1}\" AND {0}.run=\"{2}\" AND {0}.lane=\"{3}\"'.format(table, library, run, lane))
           conn.commit()
    conn.close()



if __name__ == '__main__':

    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'provReporterDb.py', description='Script to add data from FPR, Nabu and Pinery to Provenance Reporter Db')
    parser.add_argument('-f', '--fpr', dest='fpr', default = '/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz', help='Path to the File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    parser.add_argument('-n', '--nabu', dest='nabu', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
    parser.add_argument('-sp', '--sample_provenance', dest='sample_provenance', default = 'http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance', help = 'URL of the Sample Provenance in Pinery. Default is http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance')
    parser.add_argument('-pp', '--project_provenance', dest = 'project_provenance', default = 'http://pinery.gsi.oicr.on.ca/sample/projects', help = 'URL of the Project Provenance in Pinery. Default is http://pinery.gsi.oicr.on.ca/sample/projects')
    parser.add_argument('-d', '--database', dest='database', help='Path to the database file', required=True)
    
    # get arguments from the command line
    args = parser.parse_args()
    
    # initiate database
    initiate_db(args.database)
    # add or update information in tables
    add_project_info_to_db(args.database, 'Projects', args.project_provenance)
    print('added data into Projects')
    start = time.time()
    add_workflows_info_to_db(args.fpr, args.database, 'Workflows', 'Parents', 'Children')
    end = time.time()
    print('added data into Workflows', end - start)
    start = time.time()
    add_file_info_to_db(args.database, 'FilesQC', args.fpr, 'Files', args.project_provenance, args.nabu)
    end = time.time()
    print('added file info into FilesQC', end - start)
    # add_file_info_to_db(args.database, 'Files', args.fpr, 'Files', args.project_provenance, args.nabu)
    # add_library_info_to_db(args.database, 'Libraries', args.sample_provenance)
    # add_workflow_inputs_to_db(args.database, args.fpr, 'Workflow_Inputs')
