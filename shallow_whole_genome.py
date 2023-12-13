# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:39:27 2023

@author: rjovelin
"""


from utilities import connect_to_db
from whole_genome import get_parent_workflows



def get_shallow_wg(project_name, database, workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries'):
    '''
    
    
    
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
        d = {'donor': donor,
            'sample': sample,
            'library': i['library'],
            'wf': i['wf'],
            'workflow_id': i['wfrun_id']}
            
        if donor not in D:
            D[donor] = {}
        if sample not in D[donor]:
            D[donor][sample] = []
        D[donor][sample].append(d)
    
    return D




