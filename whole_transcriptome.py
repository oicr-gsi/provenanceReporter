# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:35:40 2023

@author: rjovelin
"""

from utilities import connect_to_db, get_children_workflows, filter_out_QC_workflows




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




