# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:42:36 2023

@author: rjovelin
"""


from utilities import connect_to_db, convert_epoch_time, remove_non_analysis_workflows,\
    get_children_workflows, get_workflow_names, get_donors
from networks import get_node_labels, make_adjacency_matrix, plot_workflow_network



def get_bmpp_case(project_name, case, platform, library_type, database):
    '''
    (str, str, str, str, str) -> list
    
    Returns a list of bmpp workflow Ids corresponding to the specific case from project 
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
                        Workflow_Inputs.wfrun_id FROM Libraries JOIN Workflow_Inputs JOIN Workflows \
                        WHERE Libraries.project_id = '{0}' AND Workflow_Inputs.project_id = '{0}' \
                        AND Workflows.project_id = '{0}' AND Workflows.wfrun_id = Workflow_Inputs.wfrun_id \
                        AND Workflow_Inputs.library = Libraries.library \
                        AND LOWER(SUBSTR(Workflows.wf, 1, 21)) = 'bammergepreprocessing' \
                        AND Libraries.sample ='{1}'".format(project_name, case)).fetchall()
    conn.close()

    bmpps = list(set([i['wfrun_id'] for i in data if platform in i['platform'].lower() and library_type == i['library_type']]))
    
    return bmpps



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



def get_case_call_ready_samples(project_name, bmpp_samples):
    '''
    (str, dict)
    
    Returns a dictionary with all normal and tumor samples for project_name
    processed for all bamMergePreprocessing workflow run ids for a specific case
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - bmpp_samples (dict): Dictionary with normal, tumor samples for each bmpp run id   '''
    
    
    D = {'normal': [], 'tumour': []} 
    
    for i in bmpp_samples:
        D['normal'].extend(bmpp_samples[i]['normal'])
        D['tumour'].extend(bmpp_samples[i]['tumour'])
    return D
    

def map_analysis_workflows_to_sample(project_name, sample, platform, database):
    '''
    (str, str, str, str) -> list
    
    Returns a list with workflow information for a sample from project_name with input 
    sequences sequenced on a given platform
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - sample (str): Specific sample from project
    - platform (str): Sequencing platform: novaseq, miseq, nextseq or hiseq
    - database (str): Path to the sqlite database
    '''
        
    sample = sample.split('_')
    case = '_'.join(sample[0:2])
    tissue_type = sample[2]
    tissue_origin = sample[3]
    library_type = sample[4]
    group_id = '_'.join(sample[5:])
    
    conn = connect_to_db(database)    
    data = conn.execute("SELECT Workflow_Inputs.wfrun_id, Workflow_Inputs.platform, Workflows.wf FROM \
                          Workflow_Inputs JOIN Workflows JOIN Libraries WHERE Workflow_Inputs.library = Libraries.library \
                          AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id AND \
                          Workflow_Inputs.project_id = '{0}' AND Libraries.project_id = '{0}' \
                          AND Workflows.project_id = '{0}' AND Libraries.sample = '{1}' AND \
                          Libraries.library_type = '{2}' AND Libraries.tissue_origin = '{3}' AND \
                          Libraries.tissue_type = '{4}' AND Libraries.group_id = '{5}'".format(project_name, case, library_type, tissue_origin, tissue_type, group_id)).fetchall()
    conn.close()   
    data = list(set(data))
    
    # remove non analysis worfkflows
    data = remove_non_analysis_workflows(data)
        
    to_remove = [i for i in data if platform not in i['platform'].lower()]
    for i in to_remove:
        data.remove(i)
    return data



def find_common_workflows(project_name, platform, samples, database):
    '''
    (str, str, list, str) -> list
    
    Returns a list of workflows processed with the normal and tumor sample pair in samples
    for data from a given sequencing platform and project of interest
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - platform (str): Sequencing platform: novaseq, miseq, hiseq or nextseq
    - samples (list): list with a normal and a tumor sample
    - database (str): Path to the sqlite database
    '''
    
    # get the analysis workflows for the normal and tumor sample    
    L1 = map_analysis_workflows_to_sample(project_name, samples[0], platform, database)
    L2 = map_analysis_workflows_to_sample(project_name, samples[1], platform, database)
    # get workflows in common, processed with the tumor, normal sample pair
    merged = list(set(L1).intersection(L2))
    
    return merged
    


def map_workflows_to_sample_pairs(project_name, platform, pairs, database):
    '''
    (str, str, list, str) -> dict
    
    Returns a dictionary with workflow information for each normal, tumor sample pair
    in pairs for a given project and for sequences done on a given platform
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - platform (str): Sequencing platform: novaseq, miseq, hiseq or nextseq
    - pairs (list): List of normal, tumor sample pairs
    - database (str): Path to the sqlite database
    '''
    
    D = {}
    for i in pairs:
        j = ' | '.join(sorted(i))
        L = find_common_workflows(project_name, platform, i, database)
        if len(L) != 0:
            D[j] = L 
    return D


def group_normal_tumor_pairs(samples):
    '''
    (dict) -> list
    
    Returns a list with all possible normal, tumor pairs from the dictionary
    of normal and tumor samples for a given case
    
    Parameters
    ----------
    - samples (dict): Dictionary with normal and tumour samples for a given case
    '''
    
    pairs = []
    if samples['normal'] and samples['tumour']:
        for i in samples['normal']:
            for j in samples['tumour']:
                pairs.append(sorted([i, j]))
    elif samples['normal']:
        for i in samples['normal']:
            pairs.append([i])
    elif samples['tumour']:
        for i in samples['tumour']:
            pairs.append([i])
       
    return pairs



def get_sample_bmpp(project_name, sample, platform, database):
    '''
    (str, str, str, str) -> list
    
    Returns a list of bmpp workflow ids with input data from sample from project sequenced on platform
       
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - sample (str): Normal or tumour sample identifier
    - platform (str): Sequencing platform
    - database (str): Path to the sqlite database
    '''
    
    sample = sample.split('_')
    
    case = '_'.join(sample[0:2])
    tissue_type = sample[2]
    tissue_origin = sample[3]
    library_type = sample[4]
    group_id = sample[5]
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Libraries.library, Workflow_Inputs.platform, Workflow_Inputs.wfrun_id \
                        FROM Libraries JOIN Workflow_Inputs JOIN Workflows WHERE \
                        Libraries.project_id = '{0}' AND Workflow_Inputs.project_id = '{0}' \
                        AND Workflows.project_id = '{0}' AND Libraries.sample = '{1}' AND \
                        Libraries.library_type = '{2}' AND Libraries.tissue_origin = '{3}' \
                        AND Libraries.tissue_type = '{4}' AND Libraries.group_id = '{5}' \
                        AND Workflows.wfrun_id = Workflow_Inputs.wfrun_id \
                        AND Workflow_Inputs.library = Libraries.library \
                        AND LOWER(SUBSTR(Workflows.wf, 1, 21)) = 'bammergepreprocessing'".format(project_name,
                        case, library_type, tissue_origin, tissue_type, group_id)).fetchall()
    conn.close()

    data = list(set(data))

    bmpps = list(set([i['wfrun_id'] for i in data if platform in i['platform'].lower()]))
    
    return bmpps


def map_workflows_to_parent(D, parents):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with parent workflows of each workflow 
    of the sample pairs in D
    
    Parameters
    ----------
    - D (dict): Dictionary with workflows of each normal, tumour sample pair
    - parents (dict): Dictionary with parent workflows of each workflow in a given project
    '''
    
    parent_workflows = {}
    
    for samples in D:
        for j in D[samples]:
            parent = parents[j['wfrun_id']]
            if j['wfrun_id'] not in parent_workflows:
                parent_workflows[j['wfrun_id']] = []
            for k in parent:
                parent_workflows[j['wfrun_id']].extend(parent[k])

    return parent_workflows



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
                        FROM Parents JOIN Workflows WHERE Parents.project_id = '{0}' \
                        AND Workflows.project_id = '{0}' AND Workflows.wfrun_id = Parents.parents_id;".format(project_name)).fetchall()
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
      



def find_analysis_blocks(D, parents, parent_workflows, bmpp):
    '''
    (dict, dict, dict, list) -> list
    
    Returns a dictionary with analysis workflows grouped by sample pairs and blocks 
    (ie, originating from common call-ready workflows)
    
    Parameters
    ----------
    - D (dict): Dictionary with workflows of each normal, tumour sample pair
    - parents (dict): Dictionary with parent workflows of each workflow in a given project
    - parent_worfklows (dict): Dictionary with parent workflows of each workflow for the normal, tumour sample pairs
    - bmpp (list): List of bmpp workflow run id for a given case                           
    '''
    
    # track blocks
    blocks = {}

    # record all parent workflows
    L = []

    # get the bmpp sort bmpp-dowsntream workflows by block and sample     
    for samples in D:
        for j in D[samples]:
            parent = parents[j['wfrun_id']]
            if len(parent.keys()) == 1 and parent[list(parent.keys())[0]][0] in bmpp:
                assert 'bamMergePreprocessing' in list(parent.keys())[0]
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



def list_block_workflows(blocks):
    '''
    (dict) -> dict
    
    Returns a dictionary with list of workflow ids grouped by sample pair and analysis block
    
    Parameters
    ----------
    - blocks (dict): Dictionary with analysis workflows grouped by sample pairs and blocks 
    (ie, originating from common call-ready workflows)
    '''
        
    # list all workflows downstream of each bmpp anchors
    W = {}
    for block in blocks:
        W[block] = {}
        for bmpp in blocks[block]:
            L = []
            for d in blocks[block][bmpp]:
                for workflow in d:
                    L.append(d[workflow]['parent']['wfrun_id'])
                    for k in d[workflow]['children']:
                        L.append(k['wfrun_id'])
            L.extend(bmpp.split('.'))
            L = sorted(list(set(L)))            
            W[block][bmpp] = L
    
    return W



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
    data = conn.execute("SELECT DISTINCT creation_date, wfrun_id FROM Files WHERE project_id='{0}'".format(project_name)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        D[i['wfrun_id']] = i['creation_date']
        
    return D


def get_block_analysis_date(block_workflows, creation_dates):
    '''
    (dict) -> dict, dict
    
    Returns a dictionary with the most recent analysis date of any workflow downstream
    of bmpp parent workflows for each sample pair
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
    - creation_dates (dict): Dictionary with creation date of each worklow in project
    '''

    block_date = {}
    for block in block_workflows:
        block_date[block] = {}
        for bmpp in block_workflows[block]:
            most_recent = 0
            # get the analysis date for each workflow
            for wf in block_workflows[block][bmpp]:
                if wf in creation_dates:
                    wf_date = creation_dates[wf]
                else:
                    wf_date = 'NA'
                if wf_date != 'NA':
                    if wf_date  > most_recent:
                        most_recent = wf_date
                
            if most_recent:
                block_date[block][bmpp] = convert_epoch_time(most_recent)
            else:
                block_date[block][bmpp] = 'NA'
    return block_date



def sort_call_ready_samples(project_name, blocks, bmpp_samples, workflow_names):
    '''
    (str, dict, dict, dict) -> dict
    
    Returns a dictionary with samples sorted by block and bmpp parent workflow
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - blocks (dict): Dictionary with workflows sorted by block and bmpp parent workflows
    - bmpp_samples (dict): Dictionary with normal and tumor samples for each bmpp run id
    - workflow_names (dict): Dictionary with workflow name for each workflow run id in project
    '''
    
    D = {}

    for block in blocks:
        D[block] = {}
        for bmpp in blocks[block]:
            D[block][bmpp] = {}
            bmpp_ids = bmpp.split('.')
            for i in bmpp_ids:
                D[block][bmpp][i] = {'samples': bmpp_samples[i], 'name': workflow_names[(i)]}
       
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
    data = conn.execute("SELECT DISTINCT {0}.file_count, {0}.wfrun_id FROM {0} WHERE {0}.project_id = '{1}'".format(workflow_table, project_name)).fetchall()
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
    data = conn.execute("SELECT {0}.limskey, {0}.wfrun_id FROM {0} WHERE {0}.project_id = '{1}'".format(workflow_input_table, project_name)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['wfrun_id'] not in D:
            D[i['wfrun_id']] = []
        D[i['wfrun_id']].append(i['limskey'])
            
    return D


def get_file_release_status(project_name, database):
    '''
    (str, str) -> dict
    
    Returns a dictionary with file swid and file QC status from Nabu for each limskey
    corresponding to fastq-generating workflows casava and bcl2fastq for the project of interest
       
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database 
    '''
    
    # ignore fastq-import workflows because files from these workflows are not released    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Files.file_swid, Files.limskey, FilesQC.status FROM Files JOIN \
                        FilesQC WHERE FilesQC.file_swid = Files.file_swid AND Files.project_id = '{0}' \
                        AND FilesQC.project_id = '{0}' AND LOWER(Files.workflow) IN ('casava', 'bcl2fastq');".format(project_name)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['limskey'] not in D:
            D[i['limskey']] = []
        D[i['limskey']].append([i['file_swid'], i['status']])
    return D    



def get_workflow_release_status(workflow_id, limskeys, release_status):
    '''
    (str, dict, dict) -> bool
    
    Returns True if ALL input files of workflow with workflow_id have been released
    
    Parameters
    ----------
    - workflow_id (str): Workflow run identifier
    - limskeys (dict): Dictionary with workflow run id, list of limskeys
    - release_status (dict): Dictionary with file, swid, file release status for each limskey
    '''

    # get workflow limskeys
    workflow_limskeys = limskeys[workflow_id]

    # get release status of all files for the given workflow
    status = []
    for i in workflow_limskeys:
        if i not in release_status:
            return False
        for j in release_status[i]:
            status.append(j[1])
    status = all(map(lambda x: x.lower() == 'pass', status))
    
    return status
    
    

def get_block_release_status(block_workflows, limskeys, release_status):
    '''
    (dict, dict, dict) -> dict
    
    Returns a dictionary with the release status of the input fastqs to each bmpp parent workflow
    of each sample pair
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
    - limskeys (dict): Dictionary with workflow run id, list of limskeys
    - release_status (dict): Dictionary with file, swid, file release status for each limskey
    '''
        
    D = {}
    for block in block_workflows:
        D[block] = {}
        for bmpp in block_workflows[block]:
            workflows = bmpp.split('.')
            status = all([get_workflow_release_status(j, limskeys, release_status) for j in workflows])
            # record boolean as 0 or 1
            D[block][bmpp] = int(status)
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
    data = conn.execute("SELECT DISTINCT {0}.lane_count, {0}.wfrun_id FROM {0} WHERE {0}.project_id = '{1}'".format(workflow_table, project_name)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        counts[i['wfrun_id']] = i['lane_count']
    
    return counts


# def is_block_complete(blocks, expected_workflows):
#     '''
#     (dict, list) -> dict
    
#     Returns a dictionary indicating if the downstream workflows of each parent bmpp workflows
#     are complete (ie, contains all the expected workflows) for each block (ie,sample pair)
    
#     Parameters
#     ----------
#     - blocks (dict): Dictionary of workflow information organized by sample pair and parent bmpp workflows
#     - expected_workflows (list): List of expected generic workflows names
#     '''
    
#     D = {}
    
#     for block in blocks:
#         D[block] = {}
#         for bmpp in blocks[block]:
#             if len(blocks[block][bmpp]) == 0:
#                 complete = False
#             else:
#                 complete = True
#                 c = []
        
#                 workflows = []
#                 for d in blocks[block][bmpp]:
#                     for workflow in d:
#                         workflows.append(d[workflow]['parent']['wf'])
#                         if d[workflow]['children']:
#                             for k in d[workflow]['children']:
#                                 workflows.append(k['wf'])
#                 # homogeneize workflow names by removing the matched suffix
#                 for i in range(len(workflows)):
#                     if '_' in workflows[i]:
#                         workflows[i] = workflows[i].split('_')[0]
#                 # check that all workflows are present
#                 if sorted(list(set(workflows))) != sorted(list(set(expected_workflows))):
#                     complete = False
#                 c.append(complete)
        
#                 complete = all(c)
#             # record boolean as 0 or 1
#             D[block][bmpp] = int(complete)
            
#     return D



def is_block_complete(block_workflows, expected_workflows, workflow_names):
    '''
    (dict, list) -> dict
    
    Returns a dictionary indicating if the downstream workflows of each parent bmpp workflows
    are complete (ie, contains all the expected workflows) for each block (ie,sample pair)
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary wirh list of workflows for each sample pair and parent bmpp
    - expected_workflows (list): List of expected generic workflows names
    '''
    
    D = {}
    
    for block in block_workflows:
        D[block] = {}
        for bmpp in block_workflows[block]:
            if len(block_workflows[block][bmpp]) == 0:
                complete = False
            else:
                # make a list of caller workflows (ie, remove bmpp)
                call_ready = list(map(lambda x: x.strip(), bmpp.split('.')))
                callers = set(block_workflows[block][bmpp]).difference(set(call_ready))
                workflows = [workflow_names[i][0] for i in callers]
                # homogeneize workflow names by removing the matched suffix
                for i in range(len(workflows)):
                    workflows[i] = workflows[i].split('_')[0]
                # check that all workflows are present
                complete = set(expected_workflows).issubset(set(workflows))
            # record boolean as 0 or 1
            D[block][bmpp] = int(complete)
            
    return D


def extra_workflows(block_workflows, expected_workflows):
    '''
    (dict, list) -> dict
    
    Returns a dictionary indicating if each parent bmpp workflow has extra (more than expected)
    downstream workflows for each block (ie,sample pair)
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary wirh list of workflows for each sample pair and parent bmpp
    - expected_workflows (list): List of expected generic workflows names
    '''
    
    D = {}
    
    for block in block_workflows:
        D[block] = {}
        for bmpp in block_workflows[block]:
            L = block_workflows[block][bmpp]
            if len(L) == 0:
                extra = False
            else:
                # remove call ready workflows
                call_ready = list(map(lambda x: x.strip(), bmpp.split('.')))
                callers = set(call_ready).difference(set(call_ready))
                if len(callers) > len(expected_workflows):
                    extra = True
                else:
                    extra = False
            # record boolean as 0 or 1
            D[block][bmpp] = int(extra)
             
            
    return D
    

def order_blocks(blocks, amount_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with bmpp parent workflows ordered by amount of lane data
    for each sample pair
    
    Parameters
    ----------
    - blocks (dict): Dictionary of workflow information organized by sample pair and parent bmpp workflows
    - amount_data (dict): Dictionary of amount of data for each workflow
    '''
        
    # order sub-blocks (anchored by bmpp parent workflows) for each block (ie, sample pair)
    D = {}
        
    for block in blocks:
        L = []
        for bmpp in blocks[block]:
            total = 0
            # sum all lanes for all call-ready workflows within each block
            workflows = list(map(lambda x: x.strip(), bmpp.split('.')))
            for workflow_id in workflows:
                total += amount_data[workflow_id]
            L.append([total, bmpp])
        L = sorted(L, key=lambda x: x[0], reverse=True)
        D[block] = [i[1] for i in L]
    
    return D    



def name_WGS_blocks(ordered_blocks):
    '''
    (dict) -> dict  
    
    Returns a dictionary with sub-blocks names for each sample pair (ie, block)
    
    Parameters
    ----------
    - ordered_blocks (dict): Dictionary with bmpp parent worflows ordered by amount of data for each sample pair
    '''
    
    
    names = {}
    for block in ordered_blocks:
        counter = 1
        names[block] = {}
        for i in ordered_blocks[block]:
            k = 'WGS Analysis Block {0}'.format(counter)
            names[block][i] = k
            counter += 1
    return names



def create_block_json(project_name, blocks, block, bmpp_parent, workflow_names):
    '''
    (str, dict, str, str, dict)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and bmpp parent workflow
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - blocks (dict): Dictionary with block information
    - block (str): Sample pair in blocks
    - bmpp_parent (str): bamMergePreprocessing parent workflow(s)
    - workflow_names (dict): Dictionaru with workflow name and version for each workflow in project
    '''
    
    # organize the workflows by block and samples
    D = {}
    # re-organize the sample pair
    sample_id = '.'.join(list(map(lambda x: x.strip(), block.split('|'))))
    # get the workflow ids for that block
    D[sample_id] = list_block_workflows(blocks)[block][bmpp_parent]
    
    block_data = {}
    
    for sample in D:
        if sample not in block_data:
            block_data[sample] = {}
            
        for workflow_id in D[sample]:
            workflow_name = workflow_names[workflow_id][0]
            workflow_version = workflow_names[workflow_id][1]
            block_data[sample][workflow_name] = {'workflow_id': workflow_id, 'workflow_version': workflow_version}
    
    return block_data                





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
    data = conn.execute("SELECT Libraries.library, Libraries.sample, Libraries.project_id, \
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
        # select bmpp data sequenced on novaseq
        if platform in i['platform'].lower():
            if 'bammergepreprocessing' in i['wf'].lower():
                if i['sample'] not in cases:
                    cases[i['sample']] = {'project': i['project_id'], 'samples': set(), 'libraries': set(), 'bmpp': set()}
                cases[i['sample']]['bmpp'].add(i['wfrun_id'])
                sample = '_'.join([i['sample'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
                cases[i['sample']]['samples'].add(sample)
                cases[i['sample']]['libraries'].add(i['library'])
            
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



def find_WGS_blocks(project_name, database, expected_workflows):
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
        blocks = find_case_WGS_blocks(project_name, case, database, expected_workflows)
        if blocks:
            L.append(blocks)
       
    return L



def find_case_WGS_blocks(project_name, case, database, expected_workflows):
    '''
    (str, str, str, list) -> dict
    
    Returns a dictionary with the WGS blocks for case in project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - case (str): Case in project
    - database (str): Path to the sqlite database
    - expected_workflows (list): List of expected workflow names to define a complete block
    '''
    
    # organize the WGS blocks as a dictionary
    WGS_blocks = {}
    
    # get all the bmpp runs for WG library type and Novaseq platform
    bmpp = get_bmpp_case(project_name, case, 'novaseq', 'WG', database)
    
    # proceed only in bmpp ids exist
    if bmpp:
        # get the normal and tumor samples for each bmpp id
        bmpp_samples = map_samples_to_bmpp_runs(project_name, bmpp, database)
        # identify all the samples processed
        samples = get_case_call_ready_samples(project_name, bmpp_samples)
        # proceed only if tumor/normal samples exist
        if samples['normal'] and samples['tumour']:
            # get all pairs N/T samples
            pairs = group_normal_tumor_pairs(samples)
            # find analysis workflows for each N/T pairs
            # remove sample pairs without analysis workflows
            D = map_workflows_to_sample_pairs(project_name, 'novaseq', pairs, database)
            # find the parents of each workflow
            parents = get_parent_workflows(project_name, database)
            # get the parent workflows for each block
            parent_workflows = map_workflows_to_parent(D, parents)
            # find the blocks by mapping the analysis workflows to their parent workflows    
            blocks = find_analysis_blocks(D, parents, parent_workflows, bmpp)
            # list all workflows for each block
            block_workflows = list_block_workflows(blocks)
            # get the workflow creation date for all the workflows in project
            creation_dates = get_workflows_analysis_date(project_name, database)
            # assign date to each block. most recent file creation date from all workflows within block 
            # get the date of each workflow within block
            block_date = get_block_analysis_date(block_workflows, creation_dates)
            # map each workflow run id to its workflow name
            workflow_names = get_workflow_names(project_name, database)
            # get the workflow names
            block_workflow_names = get_node_labels(block_workflows, workflow_names)
            # convert workflow relationships to adjacency matrix for each block
            matrix = make_adjacency_matrix(block_workflows, parent_workflows)
            # create figures
            figures = plot_workflow_network(matrix, block_workflow_names)
            # get the samples for each bmpp id
            samples_bmpp = sort_call_ready_samples(project_name, blocks, bmpp_samples, workflow_names)
        
            # get release status of input sequences for each block
            # get the input limskeys for each workflow in project
            limskeys = get_workflow_limskeys(project_name, database)
        
            # get the file swid and release status for each limskey for fastq-generating workflows
            # excluding fastq-import workflows
            status = get_file_release_status(project_name, database)
            release_status = get_block_release_status(block_workflows, limskeys, status)
    
            # check if blocks are complete
            complete = is_block_complete(blocks, expected_workflows, workflow_names)
            # check if blocks have extra workflows
            extra = extra_workflows(block_workflows, expected_workflows)
        
            # get the amount of data for each workflow
            amount_data = get_amount_data(project_name, database)
            # order blocks based on the amount of data
            ordered_blocks = order_blocks(blocks, amount_data)
        
            # name each block according to the selected block order
            names = name_WGS_blocks(ordered_blocks)
               
            for samples in blocks:
                WGS_blocks[samples] = {}
                for block in blocks[samples]:
                    WGS_blocks[samples][block] = {}
                    # record network image
                    WGS_blocks[samples][block]['network'] = figures[samples][block]
                    # record all workflow ids
                    WGS_blocks[samples][block]['workflows'] = block_workflows[samples][block]
                    # record release status
                    WGS_blocks[samples][block]['release'] = release_status[samples][block]
                    # record block date
                    WGS_blocks[samples][block]['date'] = block_date[samples][block]
                    # record complete status
                    WGS_blocks[samples][block]['complete'] = complete[samples][block]
                    # record extra workflow status
                    WGS_blocks[samples][block]['extra'] = extra[samples][block]
                    # reecord block name
                    WGS_blocks[samples][block]['name'] = names[samples][block]
                    # add project and case ids
                    WGS_blocks[samples][block]['project_id'] = project_name
                    WGS_blocks[samples][block]['case_id'] = case
    
    return WGS_blocks


def get_WGS_blocks_info(project_name, case, database):
    '''
    (str, str, str) -> list 
    
    Returns a list of dictionaries containing WGS block information for a given project and case
                
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - case (str): Donor id 
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT samples, bmpp_anchor, workflows, name, date, release_status, \
                        complete, extra, network from WGS_blocks WHERE project_id = '{0}' AND \
                        case_id = '{1}'".format(project_name, case)).fetchall() 
    conn.close()

    L = [dict(i) for i in data]

    D = {}
    # group by samples
    for i in L:
        samples = i['samples']
        if samples not in D:
            D[samples] = []
        # add call ready workflows
        call_ready = list(map(lambda x: x.strip(), i['bmpp_anchor'].split('.')))
        i['call_ready'] = call_ready    
        workflows = list(map(lambda x: x.strip(), i['workflows'].split(';')))
        # add caller workflows
        callers = set(workflows).difference(set(call_ready))
        i['callers'] = callers
        # map each sample to the
        bmpp_samples = map_samples_to_bmpp_runs(project_name, call_ready, 'merged.db')
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
    '''
                
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT wfrun_id, platform FROM {0} WHERE project_id = '{1}'".format(table, project_name)).fetchall() 
    conn.close()

    D = {}
    for i in data:
        D[i['wfrun_id']] = i['platform']
    
    return D
