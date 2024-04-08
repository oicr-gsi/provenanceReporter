# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:35:40 2023

@author: rjovelin
"""

import itertools
from utilities import connect_to_db, get_children_workflows, remove_non_analysis_workflows, \
     get_workflow_names, get_donors
from whole_genome import map_analysis_workflows_to_sample, get_parent_workflows, \
    map_workflows_to_parent, list_block_workflows, get_workflows_analysis_date, \
    get_block_analysis_date, sort_call_ready_samples, is_block_complete, \
    get_workflow_limskeys, get_file_release_status, get_block_release_status, is_block_clean, \
    get_amount_data, order_blocks, map_limskey_to_library, map_library_to_sample, \
    map_library_to_case, get_workflow_output    
from networks import get_node_labels, make_adjacency_matrix, plot_workflow_network


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
        tissue = ''
        if i['tissue_type'] != 'R':
            tissue = 'tumour'
        if tissue == 'tumour':
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
        if samples['tumour']:
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


def find_case_WT_blocks(project_name, case, database, expected_workflows):
    '''
    (str, str, str, list) -> dict
    
    Returns a dictionary with the WT blocks for case in project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - case (str): Case in project
    - database (str): Path to the sqlite database
    - expected_workflows (list): List of expected workflow names to define a complete block
    '''
       
    WT_blocks = {}
            
    # build the somatic calling block
    # identify all call ready star runs for novaseq
    star = get_star_case(project_name, case, 'novaseq', 'WT', database)
    
    # proceed only in star ids exist
    if star:
        # get the tumor samples for each star id
        star_samples = map_samples_to_star_runs(project_name, star, database)
        
        # identify all the samples processed
        samples = get_WT_case_call_ready_samples(project_name, star_samples)
        if samples['tumour']:
            # remove samples without analysis workflows
            D = map_workflows_to_samples(project_name, 'novaseq', samples, database)
            # find the parents of each workflow
            parents = get_parent_workflows(project_name, database)
            # get the parent workflows for each block
            parent_workflows = map_workflows_to_parent(D, parents)
            # find the blocks by mapping the analysis workflows to their parent workflows    
            blocks = find_WT_analysis_blocks(D, parents, parent_workflows, star)
            # list all workflows for each block
            block_workflows = list_block_workflows(blocks)
            # get the workflow creation date for all the workflows in project
            creation_dates = get_workflows_analysis_date(project_name, database)
            # assign date to each block. most recent file creation date from all workflows within block 
            block_date = get_block_analysis_date(block_workflows, creation_dates)
            # map each workflow run id to its workflow name
            workflow_names = get_workflow_names(project_name, database)
            # get the workflow names
            block_workflow_names = get_node_labels(block_workflows, workflow_names)
            # convert workflow relationships to adjacency matrix for each block
            matrix = make_adjacency_matrix(block_workflows, parent_workflows)
            # create figures
            figures = plot_workflow_network(matrix, block_workflow_names)
            # get the samples for each star id
            samples_star = sort_call_ready_samples(project_name, blocks, star_samples, workflow_names)
            # get release status of input sequences for each block
            # get the input limskeys for each workflow in project
            limskeys = get_workflow_limskeys(project_name, database)
            # get the file swid and release status for each limskey for fastq-generating workflows
            # excluding fastq-import workflows
            status = get_file_release_status(project_name, database)
            release_status = get_block_release_status(block_workflows, limskeys, status)
            # check if blocks are complete
            complete = is_block_complete(block_workflows, expected_workflows, workflow_names)
            # check if blocks have extra workflows
            clean = is_block_clean(block_workflows, expected_workflows)
            # get the amount of data for each block
            amount_data = get_amount_data(project_name, database)
            # order blocks based on scores
            ordered_blocks = order_blocks(blocks, amount_data, release_status, complete, clean)
            # name each block according to the selected block order
            names = name_WT_blocks(ordered_blocks)
                        
            for samples in blocks:
                WT_blocks[samples] = {}
                for block in blocks[samples]:
                    WT_blocks[samples][block] = {}
                    # record network image
                    WT_blocks[samples][block]['network'] = figures[samples][block]
                    # record all workflow ids
                    WT_blocks[samples][block]['workflows'] = block_workflows[samples][block]
                    # record release status
                    WT_blocks[samples][block]['release'] = release_status[samples][block]
                    # record block date
                    WT_blocks[samples][block]['date'] = block_date[samples][block]
                    # record complete status
                    WT_blocks[samples][block]['complete'] = complete[samples][block]
                    # record extra workflow status
                    WT_blocks[samples][block]['clean'] = clean[samples][block]
                    # record block name
                    WT_blocks[samples][block]['name'] = names[samples][block]
                    # add project and case ids
                    WT_blocks[samples][block]['project_id'] = project_name
                    WT_blocks[samples][block]['case_id'] = case
        
        return WT_blocks



def find_WT_blocks(project_name, database, expected_workflows):
    '''
    (str, str, list) -> list
    
    Returns a list of WGS blocks for donors in project in which WGS blocks exist
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - expected_workflows (list): List of expected workflow names to define a complete block
    '''
    
    # make a list of donors:
    donors = get_donors(project_name, database)
    L = []
    for case in donors:
        blocks = find_case_WT_blocks(project_name, case, database, expected_workflows)
        if blocks:
            L.append(blocks)
       
    return L


def get_WT_standard_deliverables():
    '''
    (None) -> dict
    
    Returns a dictionary with the file extensions or file endings for each workflow
    for which output files are released as part of the standard WT package
    
    Parameters
    ----------
     None
    '''
    
    deliverables = {'star': ['.bai', '.bam'],
                    'mavis': ['.zip', '.tab'],
                    'rsem': ['.genes.results', '.isoforms.results', '.transcript.bam'],
                    'starfusion': ['.tsv'],
                    'arriba': ['.tsv']}
       
    return deliverables


def create_WT_project_block_json(project_name, database, blocks, block_status, selected_workflows, workflow_names, deliverables=None):
    '''
    (str, str, dict, dict, dict, dict, None | dict)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor star parent workflow
    
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
        for sample in blocks[case]:
            # check the selection status of the block
            if block_status[case][sample] not in ['ready', 'review']:
                # block already reviewed and workflows selected
                anchor_wf = block_status[case][sample]
                                
                for workflow in blocks[case][sample][anchor_wf]['workflows']:
                    # get workflow name and version
                    workflow_name = workflow_names[workflow][0]
                    workflow_version = workflow_names[workflow][1]
                    
                    # check workflow status
                    if selected_workflows[workflow]:
                                                
                        # get workflow output files
                        # needed to sort outputs by sample pairs or by sample for call-ready workflows
                        # even if all files are recorded
                        outputfiles = workflow_outputfiles[workflow]
                        
                        # check that only workflows in standard WGS deliverables are used
                        if deliverables:
                            # keep track of the files to be released                                            
                            L = []
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
                                        sample_id = j.replace('.', ';')
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



def create_WT_block_json(database, project_name, case, blocks, sample, anchor_workflow, workflow_names, selected_workflows, selection):
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
    - sample (str): Tumor sample
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
        deliverables = get_WT_standard_deliverables()
    elif selection == 'all':
        deliverables = {}
    
    # organize the workflows by block and samples
    D = {}
    # get the workflow ids for that block
    for i in blocks[sample]:
        if i['anchor_wf'] == anchor_workflow:
            D[sample] = map(lambda x: x.strip(), i['workflows'].split(';'))
    
    block_data = {}
    for sample in D:
        for workflow_id in D[sample]:
            # check if workflow is selected
            if selected_workflows[workflow_id]:
                # get workflow name and version
                workflow_name = workflow_names[workflow_id][0]
                workflow_version = workflow_names[workflow_id][1]
                
                # get workflow output files
                # needed to sort outputs by sample pairs or by sample for call-ready workflows
                # even if all files are recorded
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
                                sample_id = '.'.join(j.split(';'))
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
                        sample_id = '.'.join(j.split(';'))
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
