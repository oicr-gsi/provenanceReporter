# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:25:46 2023

@author: rjovelin
"""


from utilities import connect_to_db



def get_call_ready_samples(project_name, bmpp_run_id):
    '''
    
    
    
    '''
    conn = connect_to_db()
    data = conn.execute("SELECT Libraries.sample, Libraries.group_id, Libraries.library, Libraries.tissue_type, \
                        Libraries.tissue_origin, Libraries.library_type \
                        FROM Libraries JOIN Workflow_Inputs WHERE Workflow_Inputs.library = Libraries.library \
                        AND Workflow_Inputs.wfrun_id = '{0}' AND Libraries.project_id = '{1}' \
                        AND Workflow_Inputs.project_id = '{1}'".format(bmpp_run_id, project_name)).fetchall()
    conn.close()

    data = list(set(data))
    
    samples = {'normal': [], 'tumour': []}
    for i in data:
        if i['tissue_type'] == 'R':
            tissue = 'normal'
        else:
            tissue = 'tumour'
        sample = '_'.join([i['sample'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
        if sample not in samples[tissue]:
            samples[tissue].append(sample)

    return samples


def get_workflow_name(wfrun_id):
    '''
    
    
    
    '''
    
    
    
    conn = connect_to_db()    
    data = conn.execute("SELECT Workflows.wf FROM Workflows WHERE Workflows.wfrun_id='{0}'".format(wfrun_id)).fetchall()
    conn.close()   
    data = list(set(data))

    assert len(data) == 1
    return data[0]['wf']




def sort_call_ready_samples(project_name, blocks):
    
    D = {}

    for block in blocks:
        bmpp_ids = block.split('.')
        for i in bmpp_ids:
            if block not in D:
                D[block] = {}
            D[block][i] = {'samples': get_call_ready_samples(project_name, i), 'name': get_workflow_name(i)}
            
    
    return D




def get_case_call_ready_samples(project_name, bmpp_ids):
    '''
    (str, list)
    
    
    '''
    
    
    L = []
    for i in bmpp_ids:
        samples = get_call_ready_samples(project_name, i)
        if samples not in L:
            L.append(samples)    
    D = {'normal': [], 'tumour': []}
    for d in L:
        D['normal'].extend(d['normal'])
        D['tumour'].extend(d['tumour'])
    return D


