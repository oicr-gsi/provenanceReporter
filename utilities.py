# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:04:02 2023

@author: rjovelin
"""

import sqlite3
import time



def connect_to_db(database='merged.db'):
    '''
    (None) -> sqlite3.Connection
    
    Returns a connection to SqLite database prov_report.db.
    This database contains information extracted from FPR
    '''
    
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    return conn


def get_children_workflows(project_name):
    '''
    (str) -> list
    
    Returns a dictionary with workflow name, list of workflow_ids that are all children of 
    workflow_id (i.e immediate downstream workflow) for a given project_name
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT DISTINCT Workflows.wf, Parents.parents_id, \
                        Parents.children_id FROM Parents JOIN Workflows \
                        WHERE Parents.project_id = '{0}' \
                        AND Workflows.project_id = '{0}' AND \
                        Workflows.wfrun_id = Parents.children_id;".format(project_name)).fetchall()
    data= list(set(data))
    conn.close()
    
    D = {}
    for i in data:
        if i['parents_id'] not in D:
            D[i['parents_id']] = []
        D[i['parents_id']].append({'wf': i['wf'], 'children_id': i['children_id']})
    
    return D


def get_workflow_names(project_name):
    '''
    (str) -> dict
    
    Returns a dictionary with workflow_id and workflow name key, value pairs
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''

    conn = connect_to_db()
    data = conn.execute("SELECT DISTINCT Workflows.wfrun_id, Workflows.wf FROM \
                        Workflows WHERE Workflows.project_id = '{0}';".format(project_name)).fetchall()
    data= list(set(data))
    conn.close()
    
    D = {}
    for i in data:
        D[i['wfrun_id']] = i['wf']
       
    return D












# def remove_non_analysis_workflows(data):
#     '''
#     (list) -> list
    
    
#     '''
    
    
#     non_analysis_workflows = ('wgsmetrics', 'insertsizemetrics', 'bamqc', 'calculatecontamination',
#                     'calculatecontamination_lane_level', 'callability', 'fastqc',
#                     'crosscheckfingerprintscollector_bam', 'crosscheckfingerprintscollector',
#                     'fingerprintcollector', 'bamqc_lane_level', 'bamqc_call_ready', 'bwamem', 
#                     'bammergepreprocessing', 'ichorcna_lane_level', 'ichorcna', 'tmbanalysis', 
#                     'casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 
#                     'import_fastq', 'dnaseqqc', 'hotspotfingerprintcollector', 
#                     'wgsmetrics_call_ready', 'rnaseqqc_lane_level', 'rnaseqqc_call_ready')
    
#     to_remove = [i for i in data if i['wf'].lower() in non_analysis_workflows]
#     for i in to_remove:
#         data.remove(i)
    
#     return data
        



def remove_non_analysis_workflows(L):
    '''
    (list) -> list
    
    Returns a list L of dictionaries with workflows, removing any non-analysis workflows
    
    Parameters
    ----------
    - L (list): List of dictionaries with workflow name and workflow ids
    '''
    
    non_analysis_workflows = ('wgsmetrics', 'insertsizemetrics', 'bamqc', 'calculatecontamination',
                              'callability', 'fastqc', 'crosscheckfingerprintscollector',
                              'fingerprintcollector', 'bwamem', 'bammergepreprocessing',
                              'ichorcna', 'tmbanalysis', 'casava', 'bcl2fastq',
                              'fileimportforanalysis', 'fileimport', 'import_fastq',
                              'dnaseqqc', 'hotspotfingerprintcollector', 'rnaseqqc')
    
    to_remove = [i for i in L if i['wf'].split('_')[0].lower() in non_analysis_workflows]
    for i in to_remove:
        L.remove(i)
    
    return L











def convert_epoch_time(epoch):
    '''
    (str) -> str
    
    Returns epoch time in readable format
    
    Parameters
    ----------
    - epoch (str)
    '''
    
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(epoch)))


def get_miso_sample_link(project_name, case):
    '''
    (str, str) -> str
    
    Returns a link to the sample MISO page
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - case (str): Sample name
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT miso FROM Samples WHERE project_id = '{0}' AND case_id = '{1}';".format(project_name, case)).fetchall()
    data = list(set(data))
    miso_link = data[0]['miso']
    
    return miso_link



def get_library_design(library_source):
    '''
    (str) -> str
    
    Returns the description of library_source as defined in MISO
    
    Parameters
    ----------
    - library_source (str): Code of the library source as defined in MISO
    '''

    library_design = {'WT': 'Whole Transcriptome', 'WG': 'Whole Genome', 'TS': 'Targeted Sequencing',
                      'TR': 'Total RNA', 'SW': 'Shallow Whole Genome', 'SM': 'smRNA', 'SC': 'Single Cell',
                      'NN': 'Unknown', 'MR': 'mRNA', 'EX': 'Exome', 'CT': 'ctDNA', 'CM': 'cfMEDIP',
                      'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq'}

    if library_source in library_design:
        return library_design[library_source]
    else:
        return None





def get_pipelines(project_name):
    '''
    (str) -> list
    
    Returns a list of pipeline names based on the library codes extracted from database for project_name
    
    Parameters
    ----------
    - project_name (str) Name of the project of interest
    '''    
    
    # connect to db
    conn = connect_to_db()
    # extract library source
    library_source = conn.execute("SELECT DISTINCT library_type FROM Files WHERE project_id = '{0}';".format(project_name)).fetchall()
    library_source = list(set([i['library_type'] for i in  list(set(library_source))]))
    # get the library definitions
    pipelines = [get_library_design(j) for j in library_source if get_library_design(j)]
    conn.close()
    
    return pipelines
