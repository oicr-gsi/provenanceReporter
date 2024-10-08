# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:42:36 2023

@author: rjovelin
"""

import os
import itertools
from utilities import connect_to_db, convert_epoch_time, remove_non_analysis_workflows,\
    get_children_workflows, get_workflow_names, get_donors



    




def get_parent_workflows(project_name, database):
    '''
    (str, str) -> dict
    
    Returns a dictionary with workflow name, list of workflow_ids that are parent
    to all each workflow (i.e immediate upstream workflow) for a given project
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Workflows.wf, Parents.parents_id, Parents.children_id \
                        FROM Parents JOIN Workflows WHERE Parents.project_id = ? \
                        AND Workflows.project_id = ? AND Workflows.wfrun_id = Parents.parents_id;", (project_name, project_name)).fetchall()
    data= list(set(data))
    conn.close()
    
    D = {}
    for i in data:
        if i['children_id'] not in D:
            D[i['children_id']] = {}
        if i['wf'] not in D[i['children_id']]:
            D[i['children_id']][i['wf']] = []
        D[i['children_id']][i['wf']].append(i['parents_id'])
    return D
      




def get_workflows_analysis_date(project_name, database):
    '''
    (str, str) -> dict
    
    Returns the creation date of any file for each workflow id for the project of interest
       
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
        
    # connect to db
    conn = connect_to_db(database)
    # extract project info
    data = conn.execute("SELECT DISTINCT creation_date, wfrun_id FROM Files WHERE project_id= ?;", (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        D[i['wfrun_id']] = i['creation_date']
        
    return D



def get_workflow_file_count(project_name, database, workflow_table='Workflows'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with the number of files for each workflow in project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - workflow_table (str): Name of the table containing the workflow information in database
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT {0}.file_count, {0}.wfrun_id FROM {0} WHERE {0}.project_id = ?;".format(workflow_table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        counts[i['wfrun_id']] = i['file_count']
    
    return counts


def get_workflow_limskeys(project_name, database, workflow_input_table='Workflow_Inputs'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with list of limskeys for each workflow id in project
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - workflow_input_table (str): Name of the table with worklow input information in the database
    '''
        
    conn = connect_to_db(database)
    query = "SELECT {0}.limskey, {0}.wfrun_id FROM {0} WHERE {0}.project_id = ?;".format(workflow_input_table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['wfrun_id'] not in D:
            D[i['wfrun_id']] = []
        D[i['wfrun_id']].append(i['limskey'])
            
    return D


   
def get_amount_data(project_name, database, workflow_table='Workflows'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with the amount of data (ie, lane count) for each workflow in project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - workflow_table (str): Name of the table containing the workflow information in database
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT {0}.lane_count, {0}.wfrun_id FROM {0} WHERE {0}.project_id = ?;".format(workflow_table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        counts[i['wfrun_id']] = i['lane_count']
    
    return counts



def create_WG_block_json(database, project_name, case, blocks, block, anchor_workflow, workflow_names, selected_workflows, selection):
    '''
    (str, str, dict, str, str, dict, dict, str)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor parent workflow (bmpp or star)
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - project_name (str): Name of project of interest
    - case (str): Donor identifier 
    - blocks (dict): Dictionary with block information
    - block (str): Sample pair in blocks
    - anchor_workflow (str): bamMergePreprocessing parent workflow(s) or star_call_ready parent workflow
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - selection (str): Include files from all selected workflows or files from the standard deliverables
                       Values: standard or all
    '''
    
    
    libraries = map_limskey_to_library(project_name, database, table='Workflow_Inputs')
    sample_names = map_library_to_sample(project_name, database, table = 'Libraries')
    donors = map_library_to_case(project_name, database, table = 'Libraries')
    workflow_outputfiles = get_workflow_output(project_name, database, libraries, sample_names, donors, 'Files')
    
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
            
    # get the deliverables
    if selection == 'standard':
        deliverables = get_WGS_standard_deliverables()
    elif selection == 'all':
        deliverables = {}
    
    # organize the workflows by block and samples
    D = {}
    # re-organize the sample pair
    sample_id = '.'.join(list(map(lambda x: x.strip(), block.split('|'))))
    # get the workflow ids for that block
    for i in blocks[block]:
        if i['anchor_wf'] == anchor_workflow:
            D[sample_id] = map(lambda x: x.strip(), i['workflows'].split(';'))
    
    block_data = {}
    for sample in D:
        for workflow_id in D[sample]:
            # check if workflow is selected
            if selected_workflows[workflow_id]:
                # get workflow name and version
                workflow_name = workflow_names[workflow_id][0]
                workflow_version = workflow_names[workflow_id][1]
                # needed to sort outputs by sample pairs or by sample for call-ready workflows
                # even if all files are recorded
                #outputfiles = get_workflow_output(project_name, case, workflow_id, database, libraries, sample_names, 'Files')
                outputfiles = workflow_outputfiles[workflow_id]
                
                # check that only workflows in standard WGS deliverables are used
                if deliverables:
                    key = workflow_name.split('_')[0].lower()
                    if key in deliverables:
                        
                        for j in outputfiles:
                            # list all deliverable files
                            L = []
                            # gather all file paths for workflow and sample(s)
                            files = [i[0] for i in outputfiles[j]]
                            # map all file endings of deliverables with files
                            groups = list(itertools.product(files, deliverables[key]))
                            # determine which files are part of the deliverables
                            F = list(map(G, groups))
                            L = [groups[k][0] for k in range(len(F)) if F[k]]
                            
                            if L:
                                sample_id = j.replace(';', '.')
                                if case not in block_data:
                                    block_data[case] = {}
                                if sample_id not in block_data[case]:
                                    block_data[case][sample_id] = {}
                                if workflow_name not in block_data[case][sample_id]:
                                    block_data[case][sample_id][workflow_name] = []
                                
                                d = {'workflow_id': workflow_id,
                                     'workflow_version': workflow_version,
                                     'files': L}
                                if d not in block_data[case][sample_id][workflow_name]:
                                    block_data[case][sample_id][workflow_name].append(d)
                                    
                else:
                    for j in outputfiles:
                        sample_id = j.replace(';', '.')
                        d =  {'workflow_id': workflow_id, 
                              'workflow_version': workflow_version,
                              'files': [i[0] for i in outputfiles[j]]}
                        if case not in block_data:
                            block_data[case] = {}
                        if sample_id not in block_data[case]:
                            block_data[case][sample_id] = {}
                        if workflow_name not in block_data[case][sample_id]:
                            block_data[case][sample_id][workflow_name] = []
                        if d not in block_data[case][sample_id][workflow_name]:
                            block_data[case][sample_id][workflow_name].append(d)
                    
    return block_data                





def create_WGS_project_block_json(project_name, database, blocks, block_status, selected_workflows, workflow_names, deliverables=None):
    '''
    (str, str, dict, dict, dict, dict, None | dict)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor bmpp parent workflow
    
    Parameters
    ----------
    - project_name (None | str): None or name of project of interest
    - database (None | str): None or path to the sqlite database
    - blocks (dict): Dictionary with block information
    - block_status (dict): Dictionary with review status of each block
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - deliverables (None | dict): None or dictionary with file extensions of standard WGS deliverables
    '''
    
    libraries = map_limskey_to_library(project_name, database, table='Workflow_Inputs')
    sample_names = map_library_to_sample(project_name, database, table = 'Libraries')
    donors = map_library_to_case(project_name, database, table = 'Libraries')
    workflow_outputfiles = get_workflow_output(project_name, database, libraries, sample_names, donors, 'Files')
      
    
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
    
    D = {}
    for case in blocks:
        for samples in blocks[case]:
            # check the selection status of the block
            if block_status[case][samples] not in ['ready', 'review']:
                # block already reviewed and workflows selected
                anchor_wf = block_status[case][samples]
                for workflow in blocks[case][samples][anchor_wf]['workflows']:
                    workflow = os.path.basename(workflow)
                    
                    # get workflow name and version
                    workflow_name = workflow_names[workflow][0]
                    workflow_version = workflow_names[workflow][1]
                               
                    # check workflow status
                    if selected_workflows[workflow]:
                        # get workflow output files
                        # needed to sort outputs by sample pairs or by sample for call-ready workflows
                        # even if all files are recorded
                        #outputfiles = get_workflow_output(project_name, case, workflow, database, libraries, sample_names, 'Files')
                        outputfiles = workflow_outputfiles[workflow]                        
                        
                        # check that only workflows in standard WGS deliverables are used
                        if deliverables:
                            key = workflow_names[workflow][0].split('_')[0].lower()
                                                       
                            if key in deliverables:
                                for j in outputfiles:
                                    # list all deliverable files
                                    L = []
                                    # gather all file paths for workflow and sample(s)
                                    files = [i[0] for i in outputfiles[j]]
                                    # map all file endings of deliverables with files
                                    groups = list(itertools.product(files, deliverables[key]))
                                    # determine which files are part of the deliverables
                                    F = list(map(G, groups))
                                    L = [groups[k][0] for k in range(len(F)) if F[k]]
                                    
                                    if L:
                                        sample_id = j.replace(';', '.')
                                        if case not in D:
                                            D[case] = {}
                                        if sample_id not in D[case]:
                                            D[case][sample_id] = {}
                                        if workflow_name not in D[case][sample_id]:
                                            D[case][sample_id][workflow_name] = []
                                        
                                        d = {'workflow_id': workflow,
                                             'workflow_version': workflow_version,
                                             'files': L}    
                                        if d not in D[case][sample_id][workflow_name]: 
                                            D[case][sample_id][workflow_name].append(d)
                        
                        else:
                            for j in outputfiles:
                                sample_id = j.replace(';', '.')
                                d =  {'workflow_id': workflow,
                                      'workflow_version': workflow_version,
                                      'files': [i[0] for i in outputfiles[j]]}
                                if case not in D:
                                    D[case] = {}
                                if sample_id not in D[case]:
                                    D[case][sample_id] = {}
                                if workflow_name not in D[case][sample_id]:
                                    D[case][sample_id][workflow_name] = []
                                if d not in D[case][sample_id][workflow_name]:
                                    D[case][sample_id][workflow_name].append(d)
                                        
    
    return D



def get_call_ready_cases(project_name, platform, library_type, database):
    '''
    (str, str, str, str) -> dict

    Returns a dictionary with samples and libraries and bmpp and downstream workflow ids for each case in a project,
    restricting data to specified platform and library type

    Parameters
    ----------
    - project_name (str): Name of the project
    - platform (str): Name of sequencing platform.
                      Accepted values: novaseq, nextseq, hiseq, miseq
    - library_type (str): 2 letters-code indicating the type of library                   
    - database (str): Path to the sqlite database
    '''

    # get all the samples for project name 
    conn = connect_to_db(database)
    query = "SELECT Libraries.library, Libraries.case_id, Libraries.project_id, \
             Libraries.ext_id, Libraries.group_id, Libraries.library_type, \
             Libraries.tissue_type, Libraries.tissue_origin, \
             Workflows.wf, Workflows.wfrun_id, Workflow_Inputs.platform \
             from Workflow_Inputs JOIN Libraries JOIN Workflows \
             WHERE Libraries.project_id = ? AND Workflow_Inputs.project_id = ? \
             AND Workflows.project_id = ? AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id \
             AND Workflow_Inputs.library = Libraries.library AND Libraries.library_type = ?;"
    data = conn.execute(query, (project_name, project_name, project_name, library_type)).fetchall()
    conn.close()
    
    cases = {}
    for i in data:
        # select bmpp data sequenced on novaseq
        if platform in i['platform'].lower():
            if 'bammergepreprocessing' in i['wf'].lower():
                if i['case_id'] not in cases:
                    cases[i['case_id']] = {'project': i['project_id'], 'samples': set(), 'libraries': set(), 'bmpp': set()}
                cases[i['case_id']]['bmpp'].add(i['wfrun_id'])
                sample = '_'.join([i['case_id'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
                cases[i['case_id']]['samples'].add(sample)
                cases[i['case_id']]['libraries'].add(i['library'])
            
    # get parent-children workflow relationships
    parents = get_children_workflows(project_name, database)
    
    # find the bmpp downstream workflows
    for sample in cases:
        downstream = []
        for bmpp in cases[sample]['bmpp']:
            if bmpp in parents:
                # get the bmpp downstream workflows
                children = parents[bmpp]
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




def get_call_ready_samples(project_name, bmpp_run_id, database):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with normal and tumour samples from project_name processed through bamMergePreprcessing
    workflow with bmpp_run_id 
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - bmpp_run_id (str): BamMergePreprocessing workflow run identifier
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    query = "SELECT Libraries.case_id, Libraries.group_id, Libraries.library, Libraries.tissue_type, \
                        Libraries.tissue_origin, Libraries.library_type \
                        FROM Libraries JOIN Workflow_Inputs WHERE Workflow_Inputs.library = Libraries.library \
                        AND Workflow_Inputs.wfrun_id = ? AND Libraries.project_id = ? \
                        AND Workflow_Inputs.project_id = ?"
    data = conn.execute(query, (os.path.basename(bmpp_run_id), project_name, project_name)).fetchall()
    conn.close()

    data = list(set(data))
    
    samples = {'normal': [], 'tumour': []}
    for i in data:
        if i['tissue_type'] == 'R':
            tissue = 'normal'
        else:
            tissue = 'tumour'
        sample = '_'.join([i['case_id'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
        if sample not in samples[tissue]:
            samples[tissue].append(sample)

    return samples



def map_samples_to_bmpp_runs(project_name, bmpp_ids, database):
    '''
    (str, list, str) -> dict
    
    Returns a dictionary with normal, tumor samples for each bmpp run id
      
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - bmpp_ids (list): List of BamMergePreprocessing workflow run identifiers for a single case
    - database (str): Path to the sqlite database
    '''

    D = {}
    for i in bmpp_ids:
        # initiate dictionary
        samples = get_call_ready_samples(project_name, i, database)
        D[i] = samples
    return D



def get_WGTS_blocks_info(project_name, case, database, table):
    '''
    (str, str, str, str) -> list 
    
    Returns a list of dictionaries containing WGS or WT block information for a given project and case
                
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - case (str): Donor id 
    - database (str): Path to the sqlite database
    - table (str): Table with block information: WGS_blocks, WT_blocks or EX_blocks
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT samples, anchor_wf, workflows, name, date, \
             complete, clean, network from {0} WHERE project_id = ? AND \
             case_id = ?;".format(table)
    data = conn.execute(query, (project_name, case)).fetchall() 
    conn.close()

    L = [dict(i) for i in data]

    D = {}
    # group by samples
    for i in L:
        samples = i['samples']
        if samples not in D:
            D[samples] = []
        # add call ready workflows
        call_ready = list(map(lambda x: x.strip(), i['anchor_wf'].split('.')))
        i['call_ready'] = call_ready    
        workflows = list(map(lambda x: x.strip(), i['workflows'].split(';')))
        # add caller workflows
        callers = set(workflows).difference(set(call_ready))
        i['callers'] = callers
        # map each sample to the
        bmpp_samples = map_samples_to_bmpp_runs(project_name, call_ready, database)
        i['pairs'] = bmpp_samples
        D[samples].append(i)
        # sort according to sub-block name
        D[samples].sort(key = lambda x: x['name'])

    return D    


def get_sequencing_platform(project_name, database, table = 'Workflow_Inputs'):
    '''
    (str, str, str) -> list 
    
    Returns a dictionary with the sequencing platform of input raw sequences
    for each workflow for project
                
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table with workflow input information
    '''
                
    conn = connect_to_db(database)
    query = "SELECT DISTINCT wfrun_id, platform FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall() 
    conn.close()

    D = {}
    for i in data:
        D[i['wfrun_id']] = i['platform']
    
    return D




def get_selected_workflows(project_name, database, table = 'Workflows'):
    '''
    (str, str, str) -> dict 
    
    Returns a dictionary with the selected status of each workflow for the given project
                  
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table with workflow information
    '''
                
    conn = connect_to_db(database)
    query = "SELECT DISTINCT wfrun_id, selected FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall() 
    conn.close()

    D = {}
    for i in data:
        D[i['wfrun_id']] = int(i['selected'])
    
    return D

    
def get_case_workflows(case, database, table = 'WGS_blocks'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary of all workflows in each block and sample pair for a given case
    
    Parameters
    ----------
    - case (str): Donor identifier
    - database (str): Path to the sqlite database
    - table (str): Name of the table storing analysis blocks
    '''
        
    conn = connect_to_db(database)    
    query = "SELECT samples, anchor_wf, workflows FROM {0} WHERE case_id = ?;".format(table)
    data = conn.execute(query, (case,)).fetchall()
    conn.close()
    
    D = {}
        
    for i in data:
        samples = i['samples']
        block = i['anchor_wf']
        workflows = i['workflows'].split(';')
        if samples not in D:
            D[samples] = {}
        if block not in D[samples]:
            D[samples][block] = []
        D[samples][block].extend(workflows)
        D[samples][block] = list(set(workflows))
    
    return D


def update_wf_selection(workflows, selected_workflows, selection_status, database, table='Workflows'):
    '''
    (list, list, dict, str, str)
    
    Update the selection status of workflows 
    
    Parameters
    ----------
    - workflows (list): List of workflows from a single analysis block for which 
                        status needs to be updated
    - selected_workflows (list): List of selected workflows from the application form for a given case
    - selection_status (dict): Selection status of all workflows for a given project
    - database (str): Path to the sqlite database
    - table (str): Table storing workflows information
    '''
    
    # update selected status
    conn = connect_to_db(database)
    for i in workflows:
        if i in selected_workflows:
            status = 1
        else:
            status = 0
        
        # update only if status has changed
        if selection_status[os.path.basename(i)] != status:
            query = 'UPDATE Workflows SET selected = ? WHERE wfrun_id = ?;'
            conn.execute(query, (status, i))
            conn.commit()
    conn.close()



def get_wgs_blocks(project, database, table = 'WGS_blocks'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with analysis block names for each case in project
    
    Parameters
    ----------
    - project (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table with analysis blocks
    '''

    conn = connect_to_db(database)    
    query = "SELECT DISTINCT * FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project,)).fetchall()
    conn.close()
    
    D = {}
    
    for i in data:
        case, samples, anchor = i['case_id'], i['samples'], i['anchor_wf']
        workflows = i['workflows'].split(';')
        clean, complete = int(i['clean']), int(i['complete'])
        if case not in D:
            D[case] = {}
        if samples not in D[case]:
            D[case][samples] = {}
        assert anchor not in D[case][samples]
        D[case][samples][anchor] = {'workflows': workflows, 'clean': clean,
                                    'complete': complete}
            
    return D
    
             
def get_block_counts(analysis_blocks):
    '''
    (dict) -> dict
    
    Returns a dictionary with block counts for each case and sample pairs in given project
    
    Parameters
    ----------
    - analysis_blocks (dict): Dictionary with analysis blocks for a given project
    '''
    
    D = {}
    
    for i in analysis_blocks:
        for j in analysis_blocks[i]:
            if i not in D:
                D[i] = {}
            assert j not in D[i]
            D[i][j] = len(analysis_blocks[i][j])
    return D  



def review_wgs_blocks(blocks, selected_workflows):
    '''
    (dict, dict) -> dict 
    
    Returns a dictionary with status for analysis blocks for each case in project
                  
    Parameters
    ----------
    - blocks (dict): 
    - selected_workflows (dict): 
    '''
    
    D = {}
        
    for case in blocks:
        if case not in D:
            D[case] = {}
        for samples in blocks[case]:
            for anchor in blocks[case][samples]:
                # do not include call-ready workflows to determine selection/review status
                # these may be shared across multiple blocks
                L = [selected_workflows[os.path.basename(i)] for i in blocks[case][samples][anchor]['workflows'] if i not in anchor] 
                
                if any(L):
                    D[case][samples] = anchor
                    break
                else:
                    if blocks[case][samples][anchor]['clean'] and \
                      blocks[case][samples][anchor]['complete']:
                          D[case][samples] = 'ready'
                          break
                    else:
                        D[case][samples] = 'review'
    return D



def map_limskey_to_library(project_name, database, table='Workflow_Inputs'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary mapping limskey ids to library ids for each workflow in project
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - database (str): Path to the sqlite database
    - table (str): Table storing the workflow input information
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT library, limskey, wfrun_id FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        workflow = i['wfrun_id']
        if workflow not in D:
            D[workflow] = {}
        assert i['limskey'] not in D[workflow]    
        D[workflow][i['limskey']] = i['library']  
    
    return D




def map_library_to_sample(project_name, database, table = 'Libraries'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary mapping sample ids to library ids
        
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - database (str): Path to the sqlite database
    - table (str): Table storing the libraries information
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT library, case_id, tissue_type, tissue_origin, \
             library_type, group_id FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    
    D = {}
    for i in data:
        donor = i['case_id']
        library = i['library']
        sample = [i['case_id'], i['tissue_type'], i['tissue_origin'],
                           i['library_type'], i['group_id']]
        if not i['group_id']:
            sample = sample[:-1]
        sample = '_'.join(sample)    
        
        if donor not in D:
            D[donor] = {}
        if library in D[donor]:
            assert D[donor][library] == sample
        else:
            D[donor][library] = sample
        
    return D



def map_library_to_case(project_name, database, table = 'Libraries'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary mapping each library to its donor identifier
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table in database storing library information.
                   Default is Libraries
    '''
    
    # get the workflow output files sorted by sample
    conn = connect_to_db(database)
    query = "SELECT DISTINCT library, case_id FROM {0} WHERE project_id = ?".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        assert i['library'] not in D
        D[i['library']] = i['case_id']
          
    return D



def get_workflow_output(project_name, database, libraries, samples, donors, table = 'Files'):
    '''
    (str, str, dict, dict, dict, str) -> dict
    
    Returns a dictionary with workflow output files sorted by sample
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - libraries (dict): Dictionary mapping libraries to limskeys
    - samples (dict): Dictionary mapping libraries to samples
    - donors (dict): Dictionary mapping libraries to donors
    - table (str): Table in database storing File information
    '''

    # get the workflow output files sorted by sample
    conn = connect_to_db(database)
    query = "SELECT DISTINCT file, limskey, file_swid, wfrun_id FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    D = {}

    for i in data:
        file = i['file']
        limskeys = i['limskey'].split(';')
        fileswid = i['file_swid']
        workflow_id = i['wfrun_id']
        libs = list(set([libraries[workflow_id][j] for j in limskeys]))
        
        #sample_names = ';'.join(sorted(list(set([samples[case][j] for j in libs]))))
        
        sample_names = ';'.join(sorted(list(set([samples[donors[j]][j] for j in libs]))))
        
        if workflow_id not in D:
            D[workflow_id] = {}
        if sample_names in D[workflow_id]:
            D[workflow_id][sample_names].append([file, fileswid])
        else:
            D[workflow_id][sample_names] = [[file, fileswid]]
    return D



def map_fileswid_to_filename(project_name, database, table='Files'):
   '''


   '''

   # get the workflow output files sorted by sample
   conn = connect_to_db(database)
   query = "SELECT DISTINCT file_swid, file FROM {0} WHERE project_id = ?;".format(table)
   data = conn.execute(query, (project_name,)).fetchall()
   conn.close()

   D = {}
   for i in data:
       D[i['file_swid']] = i['file']
   
   return D



def get_WGS_standard_deliverables():
    '''
    (None) -> dict
    
    Returns a dictionary with the file extensions or file endings for each workflow
    for which output files are released as part of the standard WGS package
    
    Parameters
    ----------
     None
    
    '''
    
    deliverables = {'bammergepreprocessing': ['.bai', '.bam'],
                    'varianteffectpredictor': ['.mutect2.filtered.vep.vcf.gz',
                                               '.mutect2.filtered.vep.vcf.gz.tbi',
                                               '.mutect2.filtered.maf.gz'],
                    'delly': ['.somatic_filtered.delly.merged.vcf.gz',
                      '.somatic_filtered.delly.merged.vcf.gz.tbi'],
                    'sequenza': ['results.zip', 'summary.pdf', 'alternative_solutions.json'],
                    'mavis': ['.tab', '.zip']}
    
    return deliverables

    
    
def get_contamination(sample_id, database, table = 'Calculate_Contamination'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with call-ready contamination and merged limskey for sample_id     
    
    Parameters
    ----------
    - sample_id (str): Sample identifier
    - database (str): Path to the sqlite database
    - table (str): Table in database storing the call-ready contamination. Default is Calculate_Contamination
    '''    
   
    conn = connect_to_db(database)
    query = "SELECT DISTINCT contamination, merged_limskey FROM {0} WHERE sample_id = ?;".format(table)
    data = conn.execute(query, (sample_id,)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['merged_limskey'] in D:
            D[i['merged_limskey']].append(i['contamination'])
        else:
            D[i['merged_limskey']] = [i['contamination']]
    for i in D:
        D[i] = max(D[i])    
    
    return D



def group_limskeys(block_limskeys):
    '''
    (list) -> list
    
    Sort the limskeys of an analysis block by sample
    
    Parameters
    ----------
    - block_limskeys (list): List of limskeys for a given block
        
    Examples
    --------
    >>> group_limskeys(['4991_1_LDI51430', '5073_4_LDI57812', '5073_3_LDI57812', '5073_2_LDI57812'])
    ['4991_1_LDI51430', '5073_4_LDI57812;5073_3_LDI57812;5073_2_LDI57812']
    '''
    
    D = {}
    for i in block_limskeys:
        j = i.split('_')[-1]
        if j in D:
            D[j].append(i)
        else:
            D[j] = [i]
        D[j].sort()
    
    L = [';'.join(D[j]) for j in D]

    return L


def get_block_level_contamination(project_name, database, blocks, sample_pair):
    '''
    (str, str,, dict, str) -> dict
    
    Returns a dictionary with contamination for each anchor workflow in the 
    analysis block for the given sample pair.
    The contamination is the maximum contamination among samples and among 
    bmpp workflows for each anhor workflow (sub-block)
    
    Parameters
    ----------
    - project_name (str): 
    - database (str): Path to the sqlite database
    - blocks (dict): Dictionary with analysis block information
    - sample_pair (str): Normal/tumor sample pair
    '''
       
    # get the block limskeys
    limskeys = get_workflow_limskeys(project_name, database, 'Workflow_Inputs')
    
    # list the bmpp anchor workflow ids
    bmpp_workflows = [blocks[sample_pair][i]['anchor_wf'] for i in range(len(blocks[sample_pair]))]
    bmpp_workflows = list(set(bmpp_workflows))
    
    # get contamination for each sample in sample pair
    contamination = {}
    for sample in sample_pair.split('|'):
        sample = sample.strip()
        contamination.update(get_contamination(sample, database, 'Calculate_Contamination'))
    
    # map each contamination to the workflow anchor id
    D = {}    
    for workflow in bmpp_workflows:
        block_conta = []
        for workflow_id in workflow.split('.'):
            # get the limskeys for that workflow
            workflow_limskeys = limskeys[os.path.basename(workflow_id)]
            # group the limskeys by sample
            workflow_limskeys = group_limskeys(workflow_limskeys)
            # use the maximum contamination of each sample 
            conta = [contamination[i] for i in workflow_limskeys if i in contamination]
            if conta:
                conta = round(max(conta) * 100, 3)
            else:
                conta = 'NA'
            block_conta.append(conta)
        # use the maximum contamination of each bmpp workflow
        if 'NA' in block_conta:
            D[workflow] = 'NA'
        else:
            D[workflow] = max(block_conta)
    
    return D       
    


def map_workflow_to_platform(project_name, database, table = 'Workflow_Inputs'):
    '''
    (str, str, str) -> dict
    
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table storing workflow_input information
    '''
    
    # connect to db
    conn = connect_to_db(database)
    # extract library source
    data = conn.execute("SELECT DISTINCT wfrun_id, limskey, platform FROM {} WHERE project_id = ?;".format(table), (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        workflow = i['wfrun_id']
        limskey = i['limskey']
        platform = i['platform']
        
        if workflow in D:
            D[workflow]['limskey'].add(limskey)
            D[workflow]['platform'].add(platform)
        
        else:
            D[workflow] = {'limskey': {limskey}, 'platform' : {platform}}
    
    return D
    

def get_sample_sequencing_amount(project_name, case, samples, database, workflow_table = 'Workflows', workflow_input = 'Workflow_Inputs', library_table = 'Libraries'):
    '''
    (str, str, str, str, str, str, str) -> dict

    Returns a dictionary with the lane counts for sequencing workflows
    for each sequencing platform for each sample in samples 

    Parameters
    ----------
    - project_name (str): Name of project of interest 
    - case (str): Donor identifier
    - samples (str): Sample or sample pair
    - database (str): Path to the sqlite database
    - workflow_table (str): Table storing the workflow information
    - workflow_input (str): Table storing the workflow input information
    - library_table (str): Table storing the library information
    '''    
        
    libraries = map_limskey_to_library(project_name, database, table = workflow_input)
    sample_names = map_library_to_sample(project_name, database, table = library_table)
    workflow_names = get_workflow_names(project_name, database)
    amount_data = get_amount_data(project_name, database, workflow_table)
    platforms = map_workflow_to_platform(project_name, database, table = workflow_input)

    sequencing_workflows = ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq')
    workflows = [i for i in workflow_names if workflow_names[i][0].lower() in sequencing_workflows]

    # find libraries for each sample
    L = {}
    for sample in samples.split('|'):
        sample = sample.strip()
        L[sample] = []
        for i in sample_names[case]:
            if sample_names[case][i] == sample:
                L[sample].append(i)
    # find sequencing workflows for each sample
    wfs = {}
    for sample in L:
        for workflow in libraries:
            for limskey in libraries[workflow]:
                if workflow in workflows and libraries[workflow][limskey] in L[sample]:
                    if sample not in wfs:
                        wfs[sample] = []
                    wfs[sample].append(workflow)
          
    # sort lane counts, workflow and limskey by sample and platform    
    lanes = {}
    for sample in wfs:
        for workflow in wfs[sample]:
            platform = platforms[workflow]['platform']
            assert len(platform) == 1
            platform = list(platform)[0]
            if sample not in lanes:
                lanes[sample] = {}
            if platform not in lanes[sample]:
                lanes[sample][platform] = {'count': 0, 'workflows': [], 'released': [], 'limskeys': []}
            lanes[sample][platform]['count'] += amount_data[workflow]    
            lanes[sample][platform]['workflows'].append(workflow)  
    
    return lanes    
       





def get_input_sequences(project_name, database):
    '''
    (str, str) -> dict
    
    Returns a dictionary with file swid for each limskey
    corresponding to fastq-generating workflows casava and bcl2fastq for the project of interest
       
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database 
    '''
    
    # ignore fastq-import workflows because files from these workflows are not released    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Files.file_swid, Files.limskey FROM Files WHERE \
                        Files.project_id = ? AND LOWER(Files.workflow) \
                        IN ('casava', 'bcl2fastq', 'fileimportforanalysis', \
                        'fileimport', 'import_fastq');", (project_name,)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['limskey'] not in D:
            D[i['limskey']] = []
        D[i['limskey']].append(i['file_swid'])
    return D    



