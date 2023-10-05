# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:35:40 2023

@author: rjovelin
"""

from utilities import connect_to_db, get_children_workflows, remove_non_analysis_workflows
from whole_genome import map_analysis_workflows_to_sample



def get_WT_call_ready_cases(project_name, platform, database, library_type='WT'):
    '''
    (str, str, str, str) -> dict

    Returns a dictionary with samples and libraries and bmpp and downstream workflow ids for each case in a project,
    restricting data to specified platform and library type

    Parameters
    ----------
    - project_name (str): Name of the project
    - platform (str): Name of sequencing platform.
                      Accepted values: novaseq, nextseq, hiseq, miseq
    - database (str): Path to the sqlite database
    - library_type (str): 2 letters-code indicating the type of library                   
    '''

    # get all the samples for project name 
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Libraries.library, Libraries.sample, Libraries.project_id, \
                          Libraries.ext_id, Libraries.group_id, Libraries.library_type, \
                          Libraries.tissue_type, Libraries.tissue_origin, \
                          Workflows.wf, Workflows.wfrun_id, Workflow_Inputs.platform \
                          from Workflow_Inputs JOIN Libraries JOIN Workflows \
                          WHERE Libraries.project_id = '{0}' AND Workflow_Inputs.project_id = '{0}' \
                          AND Workflows.project_id = '{0}' AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id \
                          AND Workflow_Inputs.library = Libraries.library AND Libraries.library_type = '{1}';".format(project_name, library_type)).fetchall()
    conn.close()

    cases = {}
    for i in data:
        # select star data sequenced on novaseq
        if platform in i['platform'].lower():
            if 'star_call_ready' in i['wf'].lower():
                if i['sample'] not in cases:
                    cases[i['sample']] = {'project': i['project_id'], 'samples': [], 'libraries': [], 'star': []}
                cases[i['sample']]['star'].append(i['wfrun_id'])
                sample = '_'.join([i['sample'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
                cases[i['sample']]['samples'].append(sample)
                cases[i['sample']]['libraries'].append(i['library'])
            
    # get parent-children workflow relationships
    parents = get_children_workflows(project_name, database)

    # find the bmpp downstream workflows
    for sample in cases:
        downstream = []
        for star in cases[sample]['star']:
            if star in parents:
                # get the star downstream workflows
                children = parents[star]
                # removed any non-analysis workflow
                children = remove_non_analysis_workflows(children)
                # list all downtream workflows
                downstream.extend([i['children_id'] for i in children])
                # get the downstream workflows of downstream workflows
                # remove non-analysis workflows
                for workflow in downstream:
                    if workflow in parents:
                        L = remove_non_analysis_workflows(parents[workflow])
                        downstream.extend([i['children_id'] for i in L])
        cases[sample]['downstream'] = list(set(downstream)) 
        
    return cases




def get_star_case(project_name, case, platform, library_type, database):
    '''
    (str, str, str, str, database) -> list
    
    Returns a list of star call ready workflow Ids corresponding to the specific case from project 
    with input sequences from platform and library_type
    
    Parameters
    ----------   
    - project_name (str): Name of project of interest
    - case (str): Donor id 
    - platform (str): Sequencing platform. Values are novaseq, miseq, nextseq and hiseq
    - library_type (str): 2-letters code describing the type of the library (eg, WG, WT,..)
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Libraries.library_type, Workflow_Inputs.platform, \
                        Workflow_Inputs.wfrun_id FROM Libraries JOIN Workflow_Inputs \
                        JOIN Workflows WHERE Libraries.project_id = '{0}' \
                        AND Workflow_Inputs.project_id = '{0}' AND Workflows.project_id = '{0}' \
                        AND Workflows.wfrun_id = Workflow_Inputs.wfrun_id \
                        AND Workflow_Inputs.library = Libraries.library \
                        AND LOWER(Workflows.wf) = 'star_call_ready' \
                        AND Libraries.sample ='{1}'".format(project_name, case)).fetchall()
    conn.close()

    stars = list(set([i['wfrun_id'] for i in data if platform in i['platform'].lower() and library_type == i['library_type']]))
    
    return stars





def get_WT_call_ready_samples(project_name, star_id, database):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with normal and tumour samples from project_name processed through bamMergePreprcessing
    workflow with bmpp_run_id 
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - star_id (str): star_call_ready workflow run identifiers
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Libraries.sample, Libraries.group_id, Libraries.library, Libraries.tissue_type, \
                        Libraries.tissue_origin, Libraries.library_type \
                        FROM Libraries JOIN Workflow_Inputs WHERE Workflow_Inputs.library = Libraries.library \
                        AND Workflow_Inputs.wfrun_id = '{0}' AND Libraries.project_id = '{1}' \
                        AND Workflow_Inputs.project_id = '{1}'".format(star_id, project_name)).fetchall()
    conn.close()

    data = list(set(data))
    
    # only considering tumour samples for WT data     
    samples = {'tumour': []}
    for i in data:
        if i['tissue_type'] != 'R':
            tissue = 'tumour'
        sample = '_'.join([i['sample'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
        if sample not in samples[tissue]:
            samples[tissue].append(sample)

    return samples



def get_WT_case_call_ready_samples(project_name, star_samples):
    '''
    (str, dict) -> dict
    
    Returns a dictionary with tumor samples processed for all star_call_ready
    workflow run ids for a specific case and project
        
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - star_samples (dict): Dictionary with tumor samples for each star run id
    '''
     
    D = {'tumour': []}
    for i in star_samples:
        D['tumour'].extend(star_samples[i]['tumour'])
    return D    
 


def map_samples_to_star_runs(project_name, star_ids, database):
    '''
    (str, list, str) -> dict
    
    Returns a dictionary with normal, tumor samples for each bmpp run id
      
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - star_ids (list): List of star_call_ready workflow run identifiers for a single case
    - database (str): Path to the sqlite database
    '''

    D = {}
    for i in star_ids:
        # initiate dictionary
        samples = get_WT_call_ready_samples(project_name, i, database)
        D[i] = samples
    return D



def map_workflows_to_samples(project_name, platform, samples, database):
    '''
    (str, str, dict, str) -> dict
    
    Returns a dictionary with workflows processed by the WT pipeline for tumor samples
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - platform (str): Sequencing platform. Values are: novaseq, miseq, hiseq and nextseq
    - samples (dict): Dictionary with tumor samples
    - database (str): Path to the sqlite database
    '''
    
    D = {}
    for i in samples['tumour']:
        L =  map_analysis_workflows_to_sample(project_name, i, platform, database)
        if L:
            D[i] = L 
    
    return D



def find_WT_analysis_blocks(D, parents, parent_workflows, star):
    '''
    Returns a dictionary with analysis workflows grouped by sample and blocks 
    (ie, originating from common call-ready workflows)
   
    Parameters
    ----------
    - D (dict): Dictionary with workflows for tumor samples
    - parents (dict): Dictionary with parent workflows of each workflow in a given project
    - parent_worfklows (dict): Dictionary with parent workflows of each workflow for the tumour samples
    - star (list): List of star workflow run id for a given case     
    '''
    
    # track blocks
    blocks = {}

    # record all parent workflows
    L = []

    # get the bmpp sort bmpp-dowsntream workflows by block and sample     
    for samples in D:
        for j in D[samples]:
            parent = parents[j['wfrun_id']]
            if len(parent.keys()) == 1 and parent[list(parent.keys())[0]][0] in star:
                assert 'star_call_ready' in list(parent.keys())[0]
                parent_workflow = '.'.join(sorted(parent[list(parent.keys())[0]]))
                if samples not in blocks:
                    blocks[samples] = {}
                if parent_workflow not in blocks[samples]:
                    blocks[samples][parent_workflow] = []
                wfrunid = j['wfrun_id']
                d = {wfrunid: {'parent': j, 'children': []}}
                blocks[samples][parent_workflow].append(d)
                L.append(wfrunid)
    
    # sort workflows downstream of callers by block and sample, 
    # keeping track of workflow aprent-child relationshsips    
    for samples in blocks:
        for parent_workflow in blocks[samples]:
            for j in D[samples]:
                if j['wfrun_id'] not in L:
                    upstream = parent_workflows[j['wfrun_id']]
                    for k in blocks[samples][parent_workflow]:
                        for m in upstream:
                            if m in k:
                                k[m]['children'].append(j)
    return blocks                


def name_WT_blocks(ordered_blocks):
    '''
    (dict) -> dict  
    
    Returns a dictionary with sub-blocks names for each sample (ie, block)
    
    Parameters
    ----------
    - ordered_blocks (dict): Dictionary with star parent worflows ordered by amount of data for each sample 
    '''
    
    
    names = {}
    for block in ordered_blocks:
        counter = 1
        names[block] = {}
        for i in ordered_blocks[block]:
            k = 'WT Analysis Block {0}'.format(counter)
            names[block][i] = k
            counter += 1
    return names
