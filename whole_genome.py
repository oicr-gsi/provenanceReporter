# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:42:36 2023

@author: rjovelin
"""


from utilities import connect_to_db, convert_epoch_time, remove_non_analysis_workflows,\
    get_children_workflows




def get_bmpp_case(project_name, case, platform, library_type):
    '''
    (str, str, str, str) -> list
    
    Returns a list of bmpp workflow Ids corresponding to the specific case from project 
    with input sequences from platform and library_type
    
    Parameters
    ----------
    - project_name (str): Bane of project of interest
    - case (str): Donor id 
    - platform (str): Sequencing platfor.
                      Values are novaseq, miseq, nextseq and hiseq
    - library_type (str): 2-letters code describing the type of the library (eg, WG, WT,..)
    '''
    
    conn = connect_to_db()
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



def get_call_ready_samples(project_name, bmpp_run_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary with normal and tumour samples from project_name processed through bamMergePreprcessing
    workflow with bmpp_run_id 
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - bmpp_run_id (str): BamMergePreprocessing workflow run identifier
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


def map_samples_to_bmpp_runs(project_name, bmpp_ids):
    '''
    (str, list) -> dict
    
    Returns a dictionary with normal, tumor samples for each bmpp run id
      
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - bmpp_ids (list): List of BamMergePreprocessing workflow run identifiers for a single case
    '''

    D = {}
    for i in bmpp_ids:
        # initiate dictionary
        samples = get_call_ready_samples(project_name, i)
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
    

def map_analysis_workflows_to_sample(project_name, sample, platform):
    '''
    (str, str, str) -> list
    
    Returns a list with workflow information for a sample from project_name with input 
    sequences sequenced on a given platform
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - sample (str): Specific sample from project
    - platform (str): Sequencing platform: novaseq, miseq, nextseq or hiseq
    '''
        
    sample = sample.split('_')
    case = '_'.join(sample[0:2])
    tissue_type = sample[2]
    tissue_origin = sample[3]
    library_type = sample[4]
    group_id = '_'.join(sample[5:])
    
    conn = connect_to_db()    
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




def find_common_workflows(project_name, platform, samples):
    '''
    (str, str, list) -> list
    
    Returns a list of workflows processed with the normal and tumor sample pair in samples
    for data from a given sequencing platform and project of interest
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - platform (str): Sequencing platform: novaseq, miseq, hiseq or nextseq
    - samples (list): list with a normal and a tumor sample
    '''
    
    # get the analysis workflows for the normal and tumor sample    
    L1 = map_analysis_workflows_to_sample(project_name, samples[0], platform)
    L2 = map_analysis_workflows_to_sample(project_name, samples[1], platform)
    # get workflows in common, processed with the tumor, normal sample pair
    merged = list(set(L1).intersection(L2))
    
    return merged
    


def map_workflows_to_sample_pairs(project_name, platform, pairs):
    '''
    (str, str, list) -> dict
    
    Returns a dictionary with workflow information for each normal, tumor sample pair
    in pairs for a given project and for sequences done on a given platform
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - platform (str): Sequencing platform: novaseq, miseq, hiseq or nextseq
    - pairs (list): List of normal, tumor sample pairs
    '''
    
    D = {}
    for i in pairs:
        j = ' | '.join(sorted(i))
        L = find_common_workflows(project_name, platform, i)
        if len(L) != 0:
            D[j] = L 
    return D



def get_case_samples(project_name, case, library_type):
    '''
    
    
    
    '''
    conn = connect_to_db()
    data = conn.execute("SELECT Libraries.sample, Libraries.group_id, Libraries.library, \
                        Libraries.tissue_type, Libraries.tissue_origin, Libraries.library_type \
                        FROM Libraries WHERE Libraries.project_id = '{0}' \
                        AND Libraries.sample = '{1}' AND Libraries.library_type = '{2}'".format(project_name, case, library_type)).fetchall()
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



def get_sample_bmpp(project_name, sample, platform):
    '''
    
    
    
    '''
    
    sample = sample.split('_')
    
    case = '_'.join(sample[0:2])
    tissue_type = sample[2]
    tissue_origin = sample[3]
    library_type = sample[4]
    group_id = sample[5]
    
    conn = connect_to_db()
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


def map_sample_pairs_to_bmpp_runs(project_name, platform, pairs):
    
    D = {}
    for samples in pairs:
        bmpp_1 = get_sample_bmpp(project_name, samples[0], platform)
        bmpp_2 = get_sample_bmpp(project_name, samples[1], platform)

        if bmpp_1 and bmpp_2:
            

            print('1', bmpp_1)
            print('2', bmpp_2)

            D['.'.join(samples)] = list(set(bmpp_1  + bmpp_2))
    
    return D





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



def get_parent_workflows(project_name):
    '''
    (str) -> dict
    
    Returns a dictionary with workflow name, list of workflow_ids that are parent
    to all each workflow (i.e immediate upstream workflow) for a given project
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    conn = connect_to_db()
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



def get_workflows_analysis_date(project_name):
    '''
    (str) -> dict
    
    Returns the creation date of any file for each workflow id for the project of interest
       
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    
    # connect to db
    conn = connect_to_db()
    # extract project info
    data = conn.execute("SELECT DISTINCT creation_date, wfrun_id FROM Files WHERE project_id='{0}'".format(project_name)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        D[i['wfrun_id']] = i['creation_date']
        
    return D



def get_block_analysis_date(block_workflows, creation_dates):
    '''
    (dict, dict) -> dict, dict
    
    Returns a dictionary with the most recent analysis date of any workflow downstream
    of bmpp parent workflows for each sample pair, and a dictionary with the analysis date 
    of each workflow downstream  of the bmpp workflows
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
    - creation_dates (dict): Dictionary with creation date of each worklow in project
    '''

    block_date = {}
    workflow_dates = {}

    for block in block_workflows:
        block_date[block] = {}
        workflow_dates[block] = {}
        for bmpp in block_workflows[block]:
            workflow_dates[block][bmpp] = {}
            most_recent = 0
            # get the analysis date for each workflow
            for wf in block_workflows[block][bmpp]:
                if wf in creation_dates:
                    wf_date = creation_dates[wf]
                else:
                    wf_date = 'NA'
                if wf_date != 'NA':
                    workflow_dates[block][bmpp][wf] = convert_epoch_time(wf_date)
                    if wf_date  > most_recent:
                        most_recent = wf_date
                else:
                    workflow_dates[block][bmpp][wf] = 'NA'
            if most_recent:
                block_date[block][bmpp] = convert_epoch_time(most_recent)
            else:
                block_date[block][bmpp] = 'NA'
    return block_date, workflow_dates




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



def get_workflow_file_count(project_name):
    '''
    (str) -> dict
    
    Returns a dictionary with the number of files for each workflow in project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT DISTINCT Files.file, Files.wfrun_id FROM Files WHERE Files.project_id = '{0}'".format(project_name)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        counts[i['wfrun_id']] = counts.get(i['wfrun_id'], 0) + 1
       
    return counts





def get_block_workflow_file_count(block_workflows, file_counts):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with the file count for each workflow of each block
    and bmpp parent_workflow
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
    - file_counts (dict): Dictionary with file count for each workflow in project
    '''
        
    D = {}
    for block in block_workflows:
        for bmpp in block_workflows[block]:
            for workflow in block_workflows[block][bmpp]:
                D[workflow] = file_counts[workflow]
    return D




# def get_workflow_limskeys(workflow_id):
#     '''
#     (str) -> list
    
#     Returns a list of input limskeys to workflow run id 
    
#     Parameters
#     ----------
#     - workflow_id (str): Workflow run identifier
#     '''
    
    
#     conn = connect_to_db()
#     data = conn.execute("SELECT Workflow_Inputs.limskey FROM Workflow_Inputs WHERE Workflow_Inputs.wfrun_id = '{0}'".format(workflow_id)).fetchall()
#     conn.close()

#     data = list(set(data))

#     limskeys = [i['limskey'] for i in data]
#     #limskeys = list(set(limskeys()))                    
        
#     return limskeys



def get_workflow_limskeys(project_name):
    '''
    (str) -> dict
    
    Returns a dictionary with list of limskeys for each workflow id in project
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
        
    conn = connect_to_db()
    data = conn.execute("SELECT Workflow_Inputs.limskey, Workflow_Inputs.wfrun_id FROM Workflow_Inputs WHERE Workflow_Inputs.project_id = '{0}'".format(project_name)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['wfrun_id'] not in D:
            D[i['wfrun_id']] = []
        D[i['wfrun_id']].append(i['limskey'])
            
    return D


def get_file_release_status(project_name):
    '''
    (str) -> dict
    
    Returns a dictionary with file swid and file QC status from Nabu for each limskey
    corresponding to fastq-generating workflows casava and bcl2fastq for the project of interest
       
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    # ignore fastq-import workflows because files from these workflows are not released    
    conn = connect_to_db()
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
            D[block][bmpp] = status
    return D



def get_amount_data(block_workflows):
    '''
    (dict) -> dict
    
    Returns a dictionary with the amount of lane data for each workflow 
    of each sample plair and parent bmpp workflow
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
    '''
    
    D = {}
    for block in block_workflows:
        D[block] = {}
        for bmpp in block_workflows[block]:
            D[block][bmpp] = {}
            for workflow_id in block_workflows[block][bmpp]:
                D[block][bmpp][workflow_id] = len(get_workflow_limskeys(workflow_id))
    return D




def is_block_complete(blocks, expected_workflows):
    '''
    (list, list) -> dict
    
    Returns a dictionary indicating if the downstream workflows of each parent bmpp workflows
    are complete (ie, contains all the expected workflows) for each block (ie,sample pair)
    
    Parameters
    ----------
    - blocks (dict): Dictionary of workflow information organized by sample pair and parent bmpp workflows
    - expected_workflows (list): List of expected generic workflows names
    '''
    
    D = {}
    
    for block in blocks:
        D[block] = {}
        for bmpp in blocks[block]:
            if len(blocks[block][bmpp]) == 0:
                complete = False
            else:
                complete = True
                c = []
        
                workflows = []
                for d in blocks[block][bmpp]:
                    for workflow in d:
                        workflows.append(d[workflow]['parent']['wf'])
                        if d[workflow]['children']:
                            for k in d[workflow]['children']:
                                workflows.append(k['wf'])
                # homogeneize workflow names by removing the matched suffix
                for i in range(len(workflows)):
                    if '_' in workflows[i]:
                        workflows[i] = workflows[i].split('_')[0]
                # check that all workflows are present
                if sorted(list(set(workflows))) != sorted(list(set(expected_workflows))):
                    complete = False
                c.append(complete)
        
                complete = all(c)
            D[block][bmpp] = complete
            
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
        
    # order sub-blocks (anchored by bmpp paret workflows) for each block (ie, sample pair)
    D = {}
        
    for block in blocks:
        L = []
        for bmpp in blocks[block]:
            total = 0
            # sum all lanes for all call-ready workflows within each block
            workflows = list(map(lambda x: x.strip(), bmpp.split('.')))
            for workflow_id in workflows:
                assert block in amount_data
                total += amount_data[block][bmpp][workflow_id]
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



def create_block_json(project_name, blocks, block, bmpp_parent):
    
    # organize the workflows by block and samples
    D = {}
    
    sample_id = '.'.join(list(map(lambda x: x.strip(), block.split('|'))))
    D[sample_id] = []
    D[sample_id].extend(bmpp_parent.split('.'))
    for d in blocks[block][bmpp_parent]:
        for workflow in d:
            D[sample_id].append(d[workflow]['parent']['wfrun_id'])
            if d[workflow]['children']:
                for k in d[workflow]['children']:
                    D[sample_id].append(k['wfrun_id'])
                    
    block_data = {}
            
    conn = connect_to_db()
        
    for sample in D:
        if sample not in block_data:
            block_data[sample] = {}
        for workflow_id in D[sample]:
            data = conn.execute("SELECT Workflows.wfrun_id, Workflows.wf, Workflows.wfv FROM Workflows \
                                WHERE Workflows.project_id = '{0}' AND Workflows.wfrun_id = '{1}';".format(project_name, workflow_id)).fetchall()
            for i in data:
                workflow_name = i['wf']
                wfrun_id = i['wfrun_id']
                workflow_version = i['wfv']
                block_data[sample][workflow_name] = {'workflow_id': wfrun_id, 'workflow_version': workflow_version}
                                            
    conn.close()                
                    
    return block_data                





def get_call_ready_cases(project_name, platform, library_type):
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
    parents = get_children_workflows(project_name)
    
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



