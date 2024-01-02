# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:39:27 2023

@author: rjovelin
"""


from utilities import connect_to_db
from whole_genome import get_parent_workflows, map_limskey_to_library, map_library_to_sample, \
    get_workflow_output    


def get_shallow_wg(project_name, database, workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries'):
    '''
    (str, str, str, str, str) -> dict
    
    Returns a dictionary with shallow whole genome data for all samples in project
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - database (str): Path to the sqlite database
    - workflow_table (str): Name of the table storing workflow information
    - wf_input_table (str): Name of the table storing information about the input data to the workflows 
    - library_table (str): Name of the table storing information about the libraries
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
    (dict, dict) -> dict
        
    Returns a dictionary with the release status of the input sequences of the ichorcna workflows
    
    Parameters
    ----------
    - swg (dict): Dictionary storing the shallow whole genome data for a given project
    - release_status (dict): Release status of individual files for a given project
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
    
    Returns a dictionary with review status of shallow whole genome data for each samples
    of each case in project
                  
    Parameters
    ----------
    - swg (dict): Dictionary storing the shallow whole genome data for a given project
    - selected_workflows (dict): Dictionary with workflow selection status
    - input_release_status (status): Dictionary with the release status of the input
                                     sequences of the ichorcna workflows
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
                    D[donor][sample] = workflow_id
                    break
                else:
                    if input_release_status[donor][sample][workflow_id]:
                        D[donor][sample] = 'ready'
                    else:
                        D[donor][sample] = 'review'
                
    return D



def get_SWG_standard_deliverables():
    '''
    (None) -> dict
    
    Returns a dictionary with the file extensions for ichorCNA workflow for which
    output files are released as part of the standard WT package
    
    Parameters
    ----------
     None
    '''
    
    deliverables = {'ichorcna': ['.params.txt']}
           
    return deliverables





def create_swg_sample_json(database, project_name, swg, case, sample, workflow_id, workflow_names, selected_workflows, selection):
    '''
    (str, str, dict, str, str, dict, dict, str)
    
    Returns a dictionary with ichorcna workflow information for a given sample 
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - project_name (str): Name of project of interest
    - case (str): Donor identifier 
    - swg (dict): Dictionary storing the shallow whole genome data for a given project
    - sample (str): Sample
    - workflow_id (str): ichorCNA workflow identifier
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - selection (str): Include files from all selected workflows or files from the standard deliverables
                       Values: standard or all
    '''
    
    # get the deliverables
    if selection == 'standard':
        deliverables = get_SWG_standard_deliverables()
    elif selection == 'all':
        deliverables = {}
    
    # organize the data
    D = {}
    
    if selected_workflows[workflow_id]:
        # get workflow name and version
        workflow_name = workflow_names[workflow_id][0]
        workflow_version = workflow_names[workflow_id][1]
                
        # get workflow output files
        libraries = map_limskey_to_library(project_name, workflow_id, database, 'Workflow_Inputs')
        sample_names = map_library_to_sample(project_name, case, database, 'Libraries')
        outputfiles = get_workflow_output(project_name, case, workflow_id, database, libraries, sample_names, 'Files')
                
        D[case] = {}
        assert sample not in D[case]
        D[case][sample] = {}
        assert workflow_name not in D[case][sample]
        D[case][sample][workflow_name] = []

                
        # check that only workflows in standard WGS deliverables are used
        if deliverables:
            # keep track of the files to be released                                            
            L = []
            key = workflow_name.split('_')[0].lower()
            if key in deliverables:
                assert sample in outputfiles
                for i in outputfiles[sample]:
                    file = i[0]
                    for file_ending in deliverables[key]:
                        if file_ending in file and file[file.rindex(file_ending):] == file_ending:
                            L.append(file)
                
                if L:
                    D[case][sample][workflow_name].append({'workflow_id': workflow_id,
                                                           'workflow_version': workflow_version,
                                                           'files': L})
        else:
            D[case][sample][workflow_name].append({'workflow_id': workflow_id, 'workflow_version': workflow_version})
    
    return D                




def create_swg_project_json(database, project_name, swg, workflow_names, selected_workflows, deliverables=None):
    '''
    (str, str, dict, dict, dict, str)
    
    Returns a dictionary with ichorcna workflow information for a given sample 
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - project_name (str): Name of project of interest
    - swg (dict): Dictionary storing the shallow whole genome data for a given project
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - deliverables (None | dict): None or dictionary with file extensions of standard WGS deliverables
    '''
    
    # organize the data
    D = {}
    
    for case in swg:
        for sample in swg[case]:
            for workflow_id in swg[case][sample]:
                if selected_workflows[workflow_id]:
                    # get workflow name and version
                    workflow_name = workflow_names[workflow_id][0]
                    workflow_version = workflow_names[workflow_id][1]
                
                    # get workflow output files
                    libraries = map_limskey_to_library(project_name, workflow_id, database, 'Workflow_Inputs')
                    sample_names = map_library_to_sample(project_name, case, database, 'Libraries')
                    outputfiles = get_workflow_output(project_name, case, workflow_id, database, libraries, sample_names, 'Files')
                
                    if case not in D:
                        D[case] = {}
                    if sample not in D[case]:
                        D[case][sample] = {}
                    assert workflow_name not in D[case][sample]
                    D[case][sample][workflow_name] = []

                
                    # check that only workflows in standard WGS deliverables are used
                    if deliverables:
                        # keep track of the files to be released                                            
                        L = []
                        key = workflow_name.split('_')[0].lower()
                        if key in deliverables:
                            assert sample in outputfiles
                            for i in outputfiles[sample]:
                                file = i[0]
                                for file_ending in deliverables[key]:
                                    if file_ending in file and file[file.rindex(file_ending):] == file_ending:
                                        L.append(file)
                
                            if L:
                                D[case][sample][workflow_name].append({'workflow_id': workflow_id,
                                                                       'workflow_version': workflow_version,
                                                                       'files': L})
                    else:
                        D[case][sample][workflow_name].append({'workflow_id': workflow_id, 'workflow_version': workflow_version})
    
    return D                
