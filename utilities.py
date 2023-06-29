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


def get_children_workflows(project_name, workflow_id):
    '''
    (str, str) -> list
    
    Returns a dictionary with workflow name, list of workflow_ids that are all children of 
    workflow_id (i.e immediate downstream workflow) for a given project_name
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    conn = connect_to_db()
    
    data = conn.execute("SELECT Workflows.wf, Parents.children_id FROM Parents JOIN Workflows \
                        WHERE Parents.project_id = '{0}' AND Workflows.project_id = '{0}' \
                        AND Parents.parents_id = '{1}' AND Workflows.wfrun_id = Parents.children_id;".format(project_name, workflow_id)).fetchall()
    data= list(set(data))
    
    D = {}
    for i in data:
        if i['wf'] in D:
            D[i['wf']].append(i['children_id'])
            D[i['wf']] = list(set(D[i['wf']]))
        else:
            D[i['wf']] = [i['children_id']]
    conn.close()
    
    return D


def filter_out_QC_workflows(project_name, workflows):
    '''
    (str, dict) -> dict
    
    Returns a dictionary of workflow name, list of workflow ids removing any QC workflows
        
    Parameters
    ----------
    - project_name (str): name of project of interest
    - workflows (dict): Dictionary of workflow, list of workflow ids that are either parent or children of an other workflow
    '''

    to_remove = [i for i in workflows if i.lower().split('_')[0] in ['wgsmetrics', 'insertsizemetrics', 
                         'bamqc', 'calculatecontamination', 'callability', 'rnaseqqc']]
    for i in to_remove:
        del workflows[i]
    return workflows

def remove_non_analysis_workflows(data):
    '''
    (list) -> list
    
    
    '''
    
    
    non_analysis_workflows = ('wgsmetrics', 'insertsizemetrics', 'bamqc', 'calculatecontamination',
                    'calculatecontamination_lane_level', 'callability', 'fastqc',
                    'crosscheckfingerprintscollector_bam', 'crosscheckfingerprintscollector',
                    'fingerprintcollector', 'bamqc_lane_level', 'bamqc_call_ready', 'bwamem', 
                    'bammergepreprocessing', 'ichorcna_lane_level', 'ichorcna', 'tmbanalysis', 
                    'casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 
                    'import_fastq', 'dnaseqqc', 'hotspotfingerprintcollector', 
                    'wgsmetrics_call_ready', 'rnaseqqc_lane_level', 'rnaseqqc_call_ready')
    
    to_remove = [i for i in data if i['wf'].lower() in non_analysis_workflows]
    for i in to_remove:
        data.remove(i)
    
    return data
        








def convert_epoch_time(epoch):
    '''
    (str) -> str
    
    Returns epoch time in readable format
    
    Parameters
    ----------
    - epoch (str)
    '''
    
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(epoch)))



def get_workflow_name(wfrun_id):
    '''
    (str) -> str

    Returns the name of the workflow defined by workflow run id    
    
    Parameters
    ----------
    - wfrun_id (str): Workflow run identifier
    '''
        
    conn = connect_to_db()    
    data = conn.execute("SELECT Workflows.wf FROM Workflows WHERE Workflows.wfrun_id='{0}'".format(wfrun_id)).fetchall()
    conn.close()   
    data = list(set(data))

    return data[0]['wf']


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

