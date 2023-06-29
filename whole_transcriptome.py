# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:35:40 2023

@author: rjovelin
"""

from utilities import connect_to_db, get_children_workflows, filter_out_QC_workflows
from whole_genome import map_analysis_workflows_to_sample, get_parent_workflows



def get_WT_call_ready_cases(project_name, platform, library_type='WT'):
    '''
    (str, str, str) -> dict

    Returns a dictionary with samples and libraries and bmpp and downstream workflow ids for each case in a project,
    restricting data to specified platform and library type

    Parameters
    ----------
    - project_name (str): Name of the project
    - platform (str): Name of sequencing platform.
                      Accepted values: novaseq, nextseq, hiseq, miseq
    - library_type (str): 2 letters-code indicating the type of library                   
    '''

    # get all the samples for project name 
    conn = connect_to_db()
    data = conn.execute("SELECT Libraries.library, Libraries.sample, Libraries.project_id, \
                          Libraries.ext_id, Libraries.group_id, Libraries.library_type, \
                          Libraries.tissue_type, Libraries.tissue_origin, \
                          Workflows.wf, Workflows.wfrun_id, Workflow_Inputs.platform \
                          from Workflow_Inputs JOIN Libraries JOIN Workflows \
                          WHERE Libraries.project_id = '{0}' AND Workflow_Inputs.project_id = '{0}' \
                          AND Workflows.project_id = '{0}' AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id \
                          AND Workflow_Inputs.library = Libraries.library AND Libraries.library_type = '{1}';".format(project_name, library_type)).fetchall()
    conn.close()

    data = list(set(data))


    cases = {}
    for i in data:
        # select star data sequenced on novaseq
        if platform in i['platform'].lower():
            if i['sample'] not in cases:
                cases[i['sample']] = {'project': i['project_id'], 'samples': [], 'libraries': [], 'star': []}
            
            sample = '_'.join([i['sample'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
            cases[i['sample']]['samples'].append(sample)
                       
            cases[i['sample']]['libraries'].append(i['library'])
            if 'star_call_ready' in i['wf'].lower():
                cases[i['sample']]['star'].append(i['wfrun_id'])
    
    for i in cases:
        cases[i]['samples'] = list(set(cases[i]['samples']))
        cases[i]['libraries'] = list(set(cases[i]['libraries']))
        cases[i]['star'] = list(set(cases[i]['star']))

    # find the bmpp downstream workflows
    for i in cases:
        downstream = []
        for j in cases[i]['star']:
            d = get_children_workflows(cases[i]['project'], j)
            # filter out QC workflows
            d = filter_out_QC_workflows(cases[i]['project'], d)
            # list all downstream workflows
            downstream_wf = [k for m in d.values() for k in m]
            downstream.extend(downstream_wf)
            # get downstream workflows of dowmstream workflows
            for k in downstream_wf:
                d = get_children_workflows(cases[i]['project'], k)
                # filter out QC workflows
                d = filter_out_QC_workflows(cases[i]['project'], d)
                df = [n for m in d.values() for n in m]
                downstream.extend(df)
        cases[i]['downstream'] = list(set(downstream)) 
        
    return cases



def get_star_case(project_name, case, platform, library_type):
    '''
    
    
    
    '''
    conn = connect_to_db()
    data = conn.execute("SELECT Libraries.sample, Libraries.library, Libraries.library_type, Workflow_Inputs.lane, \
                         Workflow_Inputs.platform, Workflow_Inputs.wfrun_id, Workflows.wf, \
                         Workflows.wfrun_id FROM Libraries JOIN Workflow_Inputs JOIN Workflows \
                         WHERE Libraries.project_id = '{0}' AND Workflow_Inputs.project_id = '{0}' \
                         AND Workflows.project_id = '{0}' AND Workflows.wfrun_id = Workflow_Inputs.wfrun_id \
                         AND Workflow_Inputs.library = Libraries.library \
                         AND LOWER(Workflows.wf) = 'star_call_ready' \
                         AND Libraries.sample ='{1}'".format(project_name, case)).fetchall()
    conn.close()

    stars = list(set([i['wfrun_id'] for i in data if platform in i['platform'].lower() and library_type == i['library_type']]))
    
    return stars


###############


def get_WT_call_ready_samples(project_name, star_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary with normal and tumour samples from project_name processed through bamMergePreprcessing
    workflow with bmpp_run_id 
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - star_id (str): star_call_ready workflow run identifiers
    '''
    
    conn = connect_to_db()
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



def get_WT_case_call_ready_samples(project_name, star_ids):
    '''
    (str, list) -> list
    
    Returns a list of tumour samples processed for all star_call_ready workflow
    run ids for a specific case
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - star_ids (list): List of star_call_ready workflow run identifiers for a single case
    '''
    
    L = []
    for i in star_ids:
        # initiate dictionary
        samples = get_WT_call_ready_samples(project_name, i)
        if samples not in L:
            L.append(samples)    
    samples = []
    for d in L:
        samples.extend(d['tumour'])
    return samples


def map_workflows_to_samples(project_name, platform, samples):
    '''
    -> dict
    
    
    
    '''
    
    D = {}
    for i in samples:
        L =  map_analysis_workflows_to_sample(project_name, i, platform)
        D[i] = L 
    to_remove = [i for i in D if len(D[i]) == 0]    
    for i in to_remove:
        del D[i]
    
    return D



def find_WT_analysis_blocks(project_name, D, parent_workflows, star):
    '''
    
    
    '''
    
    # track blocks
    blocks = {}

    # record all parent workflows
    L = []

    # get the bmpp sort bmpp-dowsntream workflows by block and sample     
    for samples in D:
        for j in D[samples]:
            parents = get_parent_workflows(project_name, j['wfrun_id'])
            if len(parents.keys()) == 1 and parents[list(parents.keys())[0]][0] in star:
                assert 'star_call_ready' in list(parents.keys())[0]
                parent_workflow = '.'.join(sorted(parents[list(parents.keys())[0]]))
                if samples not in blocks:
                    blocks[samples] = {}
                if parent_workflow not in blocks[samples]:
                    blocks[samples][parent_workflow] = []
                wfrunid = j['wfrun_id']
                d = {wfrunid: {'parent': j, 'children': []}}
                #blocks[parent_workflow][samples].append(j)
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


