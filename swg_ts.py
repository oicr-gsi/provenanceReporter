# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:39:27 2023

@author: rjovelin
"""


from utilities import connect_to_db
from whole_genome import map_limskey_to_library, map_library_to_sample, \
    get_workflow_output    


def get_swg_ts(project_name, database, workflow, workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries'):
    '''
    (str, str, str, str, str) -> dict
    
    Returns a dictionary with ichorcna (shallow whole genome) or consensus cruncher
    (targeted sequencing) data for all samples in project
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - database (str): Path to the sqlite database
    - workflow_table (str): Name of the table storing workflow information
    - workflow (str): Workflow of interest
                      Valid workflows: ichorcna, consensuscruncher
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
    L = [i for i in data if workflow in i['wf'].lower()]
    
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


def get_input_release_status(data, release_status):
    '''
    (dict, dict) -> dict
        
    Returns a dictionary with the release status of the input sequences of the ichorcna workflows
    
    Parameters
    ----------
    - data (dict): Dictionary storing the shallow whole genome or targeted sequencing
                   data for a given project
    - release_status (dict): Release status of individual files for a given project
    '''
    
    D = {}
    
    for donor in data:
        if donor not in D:
            D[donor] = {}
        for sample in data[donor]:
            if sample not in D[donor]:
                D[donor][sample] = {}
            for workflow_id in data[donor][sample]:
                status = []
                for limskey in data[donor][sample][workflow_id]['limskey']:
                    for i in release_status[limskey]:
                        status.append(i[1])
                if all(map(lambda x: x.lower() == 'pass', status)):
                    D[donor][sample][workflow_id] = 1
                else:
                    D[donor][sample][workflow_id] = 0
    return D

   

def review_data(data, selected_workflows, input_release_status):
    '''
    (dict, dict, dict) -> dict 
    
    Returns a dictionary with review status of shallow whole genome data for each samples
    of each case in project
                  
    Parameters
    ----------
    - data (dict): Dictionary storing the shallow whole genome or targeted sequencing data for a given project
    - selected_workflows (dict): Dictionary with workflow selection status
    - input_release_status (status): Dictionary with the release status of the input
                                     sequences of the ichorcna workflows
    '''
    
    D = {}
    
    for donor in data:
        if donor not in D:
            D[donor] = {}
        for sample in data[donor]:
            if sample not in D[donor]:
                D[donor][sample] = {}
            for workflow_id in data[donor][sample]:
                if selected_workflows[workflow_id]:
                    D[donor][sample] = workflow_id
                    break
                else:
                    if input_release_status[donor][sample][workflow_id]:
                        D[donor][sample] = 'ready'
                    else:
                        D[donor][sample] = 'review'
                
    return D



def get_swg_ts_standard_deliverables(datatype):
    '''
    (str) -> dict
    
    Returns a dictionary with the file extensions for ichorCNA or consensusCruncher
    workflow for which output files are released as part of the standard WT package
    
    Parameters
    ----------
    - datatype (str): The data of interest. Valid options are: swg and ts
    '''
    
    if datatype == 'swg':
        deliverables = {'ichorcna': ['.params.txt']}
    elif datatype == 'ts':
        deliverables = {'consensuscruncher': ['.bai', '.bam', '.maf.gz',
                                              '.vcf.gz', '.tbi', '.tar.gz',
                                              '.stats.txt', '.read_families.txt']}
    return deliverables




def create_swg_ts_sample_json(datatype, database, project_name, case, sample, workflow_id, workflow_names, selected_workflows, selection):
    '''
    (str, str, str, dict, str, str, dict, dict, str)
    
    Returns a dictionary with ichorcna workflow information for a given sample 
    
    Parameters
    ----------
    - datatype (str): The data of interest. Valid options are: swg and ts
    - database (str): Path to the sqlite database
    - project_name (str): Name of project of interest
    - case (str): Donor identifier
    - sample (str): Sample
    - workflow_id (str): ichorCNA workflow identifier
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - selection (str): Include files from all selected workflows or files from the standard deliverables
                       Values: standard or all
    '''
    
    # get the deliverables
    if selection == 'standard':
        deliverables = get_swg_ts_standard_deliverables(datatype)
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
            assert sample in outputfiles
            d = {'workflow_id': workflow_id,
                 'workflow_version': workflow_version,
                 'files': [i[0] for i in outputfiles[sample]]}
            D[case][sample][workflow_name].append(d)
            
    return D                




def create_swg_ts_project_json(database, project_name, data, workflow_names, selected_workflows, deliverables=None):
    '''
    (str, str, dict, dict, dict, str)
    
    Returns a dictionary with ichorcna or consensus workflow information for a given sample 
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - project_name (str): Name of project of interest
    - data (dict): Dictionary storing the shallow whole genome or targeted sequencing data for a given project
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - deliverables (None | dict): None or dictionary with file extensions of standard WGS deliverables
    '''
    
    # organize the data
    D = {}
    
    for case in data:
        for sample in data[case]:
            for workflow_id in data[case][sample]:
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
                        assert sample in outputfiles
                        d = {'workflow_id': workflow_id,
                             'workflow_version': workflow_version,
                             'files': [i[0] for i in outputfiles[sample]]}
                        D[case][sample][workflow_name].append(d)
    
    return D                




def score_workflows(data, amount_data, release_status, dates):
    '''
    (dict, dict, dict, dict) -> dict
    
    Returns a dictionary with scores (based on amount of data, release status and creation dates)
    for each workflow for each sample
       
    Parameters
    ----------
    - data (dict): Dictionary storing the shallow whole genome or targeted sequencing data for a given project
    - amount_data (dict): Dictionary of amount of data for each workflow
    - release_status (dict): Dictionary with release status for each sub-blocks
    - dates (dict): Dictionary with workflow creation dates
    '''
    
    # rank blocks according to lane count and creation date
    ranks = {}
    for case in data:
        ranks[case] = {}
        for sample in data[case]:
            ranks[case][sample] = {}
            d = {}
            e = {}
            for workflow_id in data[case][sample]:
                d[workflow_id] = amount_data[workflow_id]
                e[workflow_id] = dates[workflow_id]
            L = sorted(list(set(d.values())))
            C = sorted(list(set(e.values())))
            # get the indices (ranks) for each workflow
            # weighted by the number of values
            for workflow_id in d:
                for i in range(len(L)):
                    if L[i] == d[workflow_id]:
                        ranks[case][sample][workflow_id] = (i + 1) / len(L)
                for i in range(len(C)):
                    if C[i] == e[workflow_id]:
                        ranks[case][sample][workflow_id] = (i + 1) / len(C)
     
    # score workflows for each sample
    D = {}
    for case in data:
        for sample in data[case]:
            for workflow_id in data[case][sample]:
                score = 0
                score += ranks[case][sample][workflow_id]
                score += release_status[case][sample][workflow_id]
                if case not in D:
                    D[case] = {}
                if sample not in D[case]:
                    D[case][sample] = {}
                D[case][sample][workflow_id] = score
   
    return D



def order_workflows(data, amount_data, release_status, dates):
    '''
    (dict, dict, dict, dict) -> dict
    
    Returns a dictionary with ichorcna or consensusCruncher workflows ordered by amount of lane of data,
    creation dates and release status for each sample
    
    Parameters
    ----------
    - data (dict): Dictionary storing the shallow whole genome or targeted sequencing data for a given project
    - amount_data (dict): Dictionary of amount of data for each workflow
    - release_status (dict): Dictionary with release status for each sub-blocks
    - dates (dict): Dictionary with workflow creation dates
    '''
    
    # score the workflows
    scores = score_workflows(data, amount_data, release_status, dates)
        
    D = {}
    for case in data:
        for sample in data[case]:
            L = []
            for workflow_id in data[case][sample]:
                L.append([scores[case][sample][workflow_id], workflow_id])
            # sort workflows according to scores
            L.sort(key=lambda x: x[0])
            L.reverse()
            if case not in D:
                D[case] = {}
            assert sample not in D[case]
            D[case][sample] = [i[1] for i in L]
    return D

