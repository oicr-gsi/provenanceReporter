# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:39:27 2023

@author: rjovelin
"""


from utilities import connect_to_db
from whole_genome import get_parent_workflows



def get_shallow_wg(project_name, database, workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries'):
    '''
    (str, str, str, str, str) -> 
    
    
    
    '''
    
    conn = connect_to_db(database)
    data = conn.execute('SELECT {0}.wfrun_id, {0}.wf, {1}.library, {1}.limskey, \
                        {1}.platform, {2}.sample, {2}.tissue_type, {2}.tissue_origin, \
                        {2}.library_type, {2}.group_id FROM {0} JOIN {1} JOIN {2} WHERE \
                        {0}.wfrun_id = {1}.wfrun_id AND {0}.project_id = {1}.project_id AND \
                        {0}.project_id = {2}.project_id AND {0}.project_id="{3}" AND \
                        {1}.library = {2}.library'.format(workflow_table, wf_input_table,
                        library_table, project_name)).fetchall()
    conn.close()

    # get all the
    L = [i for i in data if 'ichorcna' in i['wf'].lower()]
    
    D = {}
    
    for i in L:
        donor = i['sample']
        sample =  '_'.join([donor, i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
        workflow_id = i['wfrun_id']
        library = i['library']
        wf = i['wf']
        platform = i['platform']
        limskey = i['limskey']
        
        if donor not in D:
            D[donor] = {}
        if sample not in D[donor]:
            D[donor][sample] = {}
        if workflow_id in D[donor][sample]:
            D[donor][sample][workflow_id]['library'].add(library)
            D[donor][sample][workflow_id]['platform'].add(platform)
            D[donor][sample][workflow_id]['limskey'].add(limskey)
            
        else:
            D[donor][sample][workflow_id] = {'donor': donor,
                                             'sample': sample,
                                             'workflow_id': workflow_id,
                                             'library': {library},
                                             'workflow': wf,
                                             'platform': {platform},
                                             'limskey': {limskey}}
    return D


def get_input_release_status(swg, release_status):
    '''
    
    
    
    

    Parameters
    ----------
    swg : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    D = {}
    
    for donor in swg:
        if donor not in D:
            D[donor] = {}
        for sample in swg[donor]:
            if sample not in D[donor]:
                D[donor][sample] = {}
            for workflow_id in swg[donor][sample]:
                status = []
                for limskey in swg[donor][sample][workflow_id]['limskey']:
                    for i in release_status[limskey]:
                        status.append(i[1])
                if all(map(lambda x: x.lower() == 'pass', status)):
                    D[donor][sample][workflow_id] = 1
                else:
                    D[donor][sample][workflow_id] = 0
    return D

  
    

def review_swg(swg, selected_workflows, input_release_status):
    '''
    (dict, dict, dict) -> dict 
    
    Returns a dictionary with status for analysis blocks for each case in project
                  
    Parameters
    ----------
    - blocks (dict): 
    - selected_workflows (dict): 
    '''
    
    D = {}
    
    for donor in swg:
        if donor not in D:
            D[donor] = {}
        for sample in swg[donor]:
            if sample not in D[donor]:
                D[donor][sample] = {}
            for workflow_id in swg[donor][sample]:
                if selected_workflows[workflow_id]:
                    D[donor][sample][workflow_id] = workflow_id 
                else:
                    if input_release_status[donor][sample][workflow_id]:
                        D[donor][sample][workflow_id] = 'ready'
                    else:
                        D[donor][sample][workflow_id] = 'review'
                
    return D



