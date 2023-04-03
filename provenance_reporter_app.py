# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:32:40 2022

@author: rjovelin
"""

import sqlite3
import json
from flask import Flask, render_template, request, url_for, flash, redirect, Response, send_file
from werkzeug.exceptions import abort
import requests
import gzip
import os
import time
import pandas as pd




def connect_to_db():
    '''
    (None) -> sqlite3.Connection
    
    Returns a connection to SqLite database prov_report.db.
    This database contains information extracted from FPR
    '''
    
    conn = sqlite3.connect('merged.db')
    conn.row_factory = sqlite3.Row
    return conn



def group_sequences(L):
    '''
    (list) -> list

    Returns a list sequence file information by grouping paired fastqs    
    Pre-condition: all fastqs are paired-fastqs. Non-paired-fastqs are discarded.
    
    Parameters
    ----------
    - L (list): List of sqlite3.Row extracted from the database and containing sequence file information
    '''
    
    # sort list according to files
    L.sort(key = lambda x: x['file'])
    
    F = []
    
    for i in range(len(L) -1):
        # check if adjacent files are paired
        case1, case2 = L[i]['sample'], L[i+1]['sample']
        sample1, sample2 = case1 + '_' + L[i]['ext_id'], case2 + '_' + L[i+1]['ext_id']
        library1, library2 = L[i]['library'] + '_' + L[i]['group_id'], L[i+1]['library'] + '_' + L[i+1]['group_id']
        workflow1, workflow2 = L[i]['workflow'] + '_' + L[i]['version'], L[i+1]['workflow'] + '_' + L[i+1]['version']
        wfrun1, wfrun2 = L[i]['wfrun_id'], L[i+1]['wfrun_id']      
        file1, file2 = L[i]['file'], L[i+1]['file']
        run1, run2 = L[i]['run'] + '_' + str(L[i]['lane']), L[i+1]['run'] + '_' + str(L[i+1]['lane'])
        platform1, platform2 = L[i]['platform'], L[i+1]['platform']
        status1, status2 = L[i]['status'], L[i+1]['status']
        ticket1, ticket2 = L[i]['ticket'], L[i+1]['ticket'] 
        read_count1 = json.loads(L[i]['attributes'])['read_count'] if 'read_count' in json.loads(L[i]['attributes']) else 'NA' 
        read_count2 = json.loads(L[i+1]['attributes'])['read_count'] if 'read_count' in json.loads(L[i+1]['attributes']) else 'NA' 

        if case1 == case2 and run1 == run2 and platform1 == platform2 \
        and library1 == library2 and sample1 == sample2 and wfrun1 == wfrun2:
            assert read_count1 == read_count2 
            assert workflow1 == workflow2
            assert json.loads(L[i]['attributes'])['read_number'] == '1' and json.loads(L[i+1]['attributes'])['read_number'] == '2'            
            if 'GDR' in ticket1:
                ticket1 = os.path.join('https://jira.oicr.on.ca/browse/', os.path.basename(ticket1))    
            else:
                ticket1 = ''
            d = {'case': case1, 'sample': sample1, 'library': library1, 'run': run1,
                 'files': [file1, file2], 'read_count': read_count1, 'workflow': workflow1,
                 'release_status': status1, 'ticket': ticket1}
            F.append(d)
       
    F.sort(key = lambda x: x['case'])
     
    return F


def get_library_design(library_source):
    '''
    (str) -> str
    
    Returns the description of library_source as defined in MISO
    
    Parameters
    ----------
    - library_source (str): Code of the library source as defined in MISO
    '''

    library_design = {'WT': 'Whole Transcriptome', 'WG': 'Whole Genome', 'TS': 'Targeted Sequencing',
                      'TR': 'Total RNA', 'SW': 'Shallow Whole Genome', 'SM': 'smRNA', 'SC': 'Single Cell',
                      'NN': 'Unknown', 'MR': 'mRNA', 'EX': 'Exome', 'CT': 'ctDNA', 'CM': 'cfMEDIP',
                      'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq'}

    if library_source in library_design:
        return library_design[library_source]
    else:
        return None


def get_project_info(project_name):
    '''
    (str) -> list
    
    Returns a list with project information extracted from database for project_name 
    
    Parameters
    ----------
    - project_name 9str): Project of interest
    '''
    
    # connect to db
    conn = connect_to_db()
    # extract project info
    projects = conn.execute('SELECT * FROM Projects').fetchall()
    conn.close()
    
    # keep info for project_name
    project = [i for i in projects if i['project_id'] == project_name][0]
        
    return project
    
    
    
    
def get_pipelines(project_name):
    '''
    (str) -> list
    
    Returns a list of pipeline names based on the library codes extracted from database for project_name
    
    Parameters
    ----------
    - project_name (str) Name of the project of interest
    '''    
    
    # connect to db
    conn = connect_to_db()
    # extract library source
    library_source = conn.execute("SELECT library_type FROM Files WHERE project_id = '{0}';".format(project_name)).fetchall()
    # get the library definitions
    pipelines = [get_library_design(j) for j in sorted(list(set([i['library_type'] for i in library_source]))) if get_library_design(j)]
    conn.close()
    
    return pipelines
    


def get_call_ready_cases(data):
    '''
    (list) -> dict

    Returns a dictionary with samples and libraries and bamMergePreprocessing
    workflow id for each case in a project

    Parameters
    ----------
    - data (list): List of sqlite3.Row extracted from the database and sample information for bamMergePreprocessing workflow iterations
    '''
    
    L = []
    for i in data:
        if dict(i) not in L:
            L.append(dict(i))
    
    # group info for each case
    D = {}
    for i in L:
        case = i['sample']
        creation_date = i['creation_date']
        wfrun_id = i['wfrun_id']
        sample = '_'.join([case, i['tissue_origin'], i['tissue_type'], i['group_id']])
        if i['group_id']:
            sample = '_'.join([case, i['tissue_origin'], i['tissue_type'], i['group_id']])
        else:
            sample = '_'.join([case, i['tissue_origin'], i['tissue_type']])
                
        library = i['library'] 
        # keep only bmpp workflows processing sequences from Novaseq instruments
        platform = i['platform']
        if 'novaseq' in platform.lower():
            if case in D:
                # get the sample name and library name
                # compare creation time and workflow run_id
                if D[case]['creation_date'] < creation_date:
                    # update with most recent information
                    D[case]['samples'] = [sample]
                    D[case]['libraries'] = [library]
                    D[case]['wfrun_id'] = wfrun_id
                    D[case]['creation_sate'] = creation_date
                elif D[case]['creation_date'] == creation_date:
                    if D[case]['wfrun_id'] == wfrun_id:
                        D[case]['samples'].append(sample)
                        D[case]['libraries'].append(library)
            else:
                D[case] = {'libraries': [library],
                           'samples': [sample],
                           'wfrun_id': wfrun_id,
                           'creation_date': creation_date}
    for i in D:
        D[i]['libraries'] = list(set(D[i]['libraries']))
        D[i]['samples'] = list(set(D[i]['samples']))
        
    return D


def get_bmpp_files(data):
    '''
    (list) -> (str, str, list)

    Returns a tuple with the bmpp worflow run id, bmpp files and info about merged samples and libraries
    for a single case in a project

    Parameters
    ----------
    - data (list): List of sqlite3.Row extracted from the database for a single case for bamMergePreprocessing workflow iterations
    '''
       
    L = []
    for i in data:
        if dict(i) not in L:
            L.append(dict(i))
    
    # group info for each case
    D = {}
    for i in L:
        case = i['sample']
        creation_date = i['creation_date']
        wfrun_id = i['wfrun_id']
        if i['group_id']:
            sample = '_'.join([case, i['tissue_origin'], i['tissue_type'], i['group_id']])
        else:
            sample = '_'.join([case, i['tissue_origin'], i['tissue_type']])
        tissue_type = i['tissue_type']
        tissue_origin = i['tissue_origin']
        library = i['library'] 
        file = i['file']
        # keep only bmpp workflows processing sequences from Novaseq instruments
        platform = i['platform']
        if 'novaseq' in platform.lower():
            if case in D:
                # get the sample name and library name
                # compare creation time and workflow run_id
                if D[case]['creation_date'] < creation_date:
                    # update with most recent information
                    D[case]['samples'][sample] = {'libraries': [library], 'tissue_type': tissue_type, 'tissue_origin': tissue_origin}
                    D[case]['wfrun_id'] = wfrun_id
                    D[case]['files'] = [file]
                    D[case]['creation_date'] = creation_date
                elif D[case]['creation_date'] == creation_date:
                    if D[case]['wfrun_id'] == wfrun_id:
                        if sample in D[case]['samples']:
                            D[case]['samples'][sample]['libraries'].append(library)
                            assert D[case]['samples'][sample]['tissue_type'] == tissue_type 
                            assert D[case]['samples'][sample]['tissue_origin'] == tissue_origin
                            D[case]['files'].append(file)
                        else:
                            D[case]['samples'][sample] = {'libraries': [library], 'tissue_type': tissue_type, 'tissue_origin': tissue_origin}
                            D[case]['files'].append(file)
            else:
                D[case] = {'samples': {sample: {'libraries': [library], 'tissue_type': tissue_type, 'tissue_origin': tissue_origin}},
                           'wfrun_id': wfrun_id,
                           'creation_date': creation_date,
                           'files': [file]}
    for case in D:
        D[case]['files'] = [os.path.dirname(D[case]['files'][0])] + sorted(list(map(lambda x: os.path.basename(x), list(set(D[case]['files'])))))
        for sample in D[case]['samples']:
            D[case]['samples'][sample]['libraries'] = list(set(D[case]['samples'][sample]['libraries']))
    
    
    case = list(D.keys())
    assert len(case) == 1
    case = case[0]
        
    bmpp_files = '\n'.join(D[case]['files'])
        
    wfrun_id = D[case]['wfrun_id']
    # organize library and sample info for table
    L = []
    for sample in D[case]['samples']:
        for library in D[case]['samples'][sample]['libraries']:
            line = [sample, D[case]['samples'][sample]['tissue_type'],
                    D[case]['samples'][sample]['tissue_origin'], library]
            L.append(line)
    
    return wfrun_id, bmpp_files, L



def get_parent_workflows(project_name, workflow_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary with workflow name, list of workflow_ids that are all parent of 
    workflow_id (i.e immediate upstream workflow) for a given project_name
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT Workflows.wf, Parents.parents_id FROM Parents JOIN Workflows \
                        WHERE Parents.project_id = '{0}' AND Workflows.project_id = '{0}' \
                        AND Parents.children_id = '{1}' AND Workflows.wfrun_id = Parents.parents_id;".format(project_name, workflow_id)).fetchall()
    data= list(set(data))
    
    D = {}
    for i in data:
        if i['wf'] in D:
            D[i['wf']].append(i['parents_id'])
            D[i['wf']] = list(set(D[i['wf']]))
        else:
            D[i['wf']] = [i['parents_id']]
    conn.close()
    
    return D


def get_children_workflows(project_name, workflow_id):
    '''
    (str, str) -> list
    
    Returns a dictionary with workflow name, list of workflow_ids that are all children of 
    workflow_id (i.e immediate downstream workflow) for a given project_name
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    conn = connect_to_db()
    
    data = conn.execute("SELECT Workflows.wf, Parents.children_id FROM Parents JOIN Workflows \
                        WHERE Parents.project_id = '{0}' AND Workflows.project_id = '{0}' \
                        AND Parents.parents_id = '{1}' AND Workflows.wfrun_id = Parents.children_id;".format(project_name, workflow_id)).fetchall()
    data= list(set(data))
    
    D = {}
    for i in data:
        if i['wf'] in D:
            D[i['wf']].append(i['children_id'])
            D[i['wf']] = list(set(D[i['wf']]))
        else:
            D[i['wf']] = [i['children_id']]
    conn.close()
    
    return D


def get_workflow_files(project_name, workflow_id):
    '''
    (str, str) -> list, str
    
    Returns a list of files generated by workflow_id for project_name and the date when 
    the workflow completed
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - worflow_id (str): Workflow run id
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT Files.file, Files.creation_date FROM Files WHERE Files.project_id = '{0}' AND \
                        Files.wfrun_id = '{1}';".format(project_name, workflow_id)).fetchall()
    data = list(set(data))   
    conn.close()

    L = []
    creation_date = ''
    if data:
        L = [i['file'] for i in data]
        creation_date = data[0]['creation_date']
    
    return L, creation_date


def filter_out_QC_workflows(project_name, workflows):
    '''
    (str, dict) -> dict
    
    Returns a dictionary of workflow name, list of workflow ids removing any QC workflows
        
    Parameters
    ----------
    - project_name (str): name of project of interest
    - workflows (dict): Dictionary of workflow, list of workflow ids that are either parent or children of an other workflow
    '''

    to_remove = [i for i in workflows if i.lower().split('_')[0] in ['wgsmetrics', 'insertsizemetrics', 
                         'bamqc', 'calculatecontamination', 'callability']]
    for i in to_remove:
        del workflows[i]
    return workflows

    
def bmpp_input_raw_seq_status(project_name, bmpp_id):
    '''
    (str, str) -> bool

    Returns True if all the input fastqs, excepting fastqs from import workflows, 
    to the bmpp workflow run bmpp_id have been released and False otherwise
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    # get bwamem input workflow ids
    d = get_parent_workflows(project_name, bmpp_id)
    bwamem_ids = d['bwaMem']
    # get the fastq-generating worflow ids
    fastqs_workflow_ids = []
    for workflow_id in bwamem_ids:
        d = get_parent_workflows(project_name, workflow_id)
        for i in d:
            fastqs_workflow_ids.extend(d[i])
    fastqs_workflow_ids = list(set(fastqs_workflow_ids))
    
    conn = connect_to_db()
    
    # track release status of all fastqs 
    D = {}
    
    # get the file swids of the fastq-generating workflows
    for workflow_id in fastqs_workflow_ids:
        data = conn.execute("SELECT Workflows.wf, Files.file_swid, FilesQC.status  \
                              FROM Workflows JOIN Files JOIN FilesQC WHERE Files.project_id = '{0}' \
                              AND Workflows.project_id = '{0}' AND FilesQC.project_id = '{0}' AND  \
                              Files.wfrun_id = '{1}' AND FilesQC.file_swid = Files.file_swid AND \
                              Workflows.wfrun_id = '{1}';".format(project_name, workflow_id)).fetchall()
            
        # # skip import workflows because fastqs from these workflow may not need to be shared back
        for i in data:
            if 'import' not in i['wf'].lower():
                assert i['file_swid'] not in D
                D[i['file_swid']] = i['status']
    
    
    conn.close()
    
    if D:
        return all(map(lambda x: x.lower() == 'pass', list(D.values())))
    else:
        return False
       
        
def get_workflow_info(project_name, workflow_id):
    '''
    (str, str) -> dict

    Returns a dictionary with information about workflow specified by workflow_id in project_name
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    conn = connect_to_db()
    
    D = {}
    
    # grab information about bmpp downstream workflow_id excluding QC workflows
    data = conn.execute("SELECT Workflows.wfrun_id, Workflows.wfv, Workflows.attributes, Workflow_Inputs.library, Libraries.sample, Libraries.ext_id, \
                        Libraries.group_id, Libraries.tissue_type, Libraries.tissue_origin, \
                        Workflows.wf FROM Libraries JOIN Workflows JOIN Workflow_Inputs \
                        WHERE Workflow_Inputs.project_id = '{0}' AND Libraries.project_id = '{0}' AND \
                        Workflows.wfrun_id = '{1}' AND Workflow_Inputs.wfrun_id = '{1}' AND Workflow_inputs.Library = Libraries.library AND \
                        LOWER(Workflows.wf) NOT IN ('wgsmetrics_call_ready', 'insertsizemetrics_call_ready', \
                        'bamqc_call_ready');".format(project_name, workflow_id)).fetchall()
    
    
    data = list(set(data))
        
    if data:
        assert len(data) == 2
        # consider only tumor/normal pairs
        if data[0]['tissue_type'] == 'R':
            normal, tumour = 0, 1
        elif data[1]['tissue_type'] == 'R':
            normal, tumour = 1, 0
                    
        normal_sample = '_'.join([data[normal]['sample'], data[normal]['tissue_origin'], data[normal]['tissue_type'], data[normal]['group_id']])
        tumour_sample = '_'.join([data[tumour]['sample'], data[tumour]['tissue_origin'], data[tumour]['tissue_type'], data[tumour]['group_id']])
        sample = normal_sample +';' + tumour_sample
                    
        libraries = ';'.join([data[normal]['library'], data[tumour]['library']])

        attributes = data[0]['attributes'].replace("\\\"", "").replace('\\', '')
        attributes = json.loads(attributes)
        if 'reference' in attributes:
            attributes = attributes['reference'].replace('"', '')
                
        D[libraries] = {'sample': sample, 'workflow': data[0]['wf'],
                        'workflow_id': data[0]['wfrun_id'], 'version': data[0]['wfv'],
                        'attributes': attributes}
           
    conn.close()
    
    return D
        
           

def get_bmpp_downstream_workflows(project_name, bmpp_id):
    '''
    (str, str) -> dict

    Returns a dictionary with information about dowmstream bmpp_id workflows in project_name
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    # get the bmpp downstream workflows
    downstream_workflows = get_children_workflows(project_name, bmpp_id)
    # filter out QC workflows
    downstream_workflows = filter_out_QC_workflows(project_name, downstream_workflows)

    D = {}
    
    for workflow in downstream_workflows:
        for workflow_id in downstream_workflows[workflow]:
            # group all downstream workflows per library pair
            d = get_workflow_info(project_name, workflow_id)
            assert d
            libraries = list(d.keys())[0]
            if libraries not in D:
                D[libraries] = {}
            D[libraries][d[libraries]['workflow']] = d[libraries] 
            # get the downstream workflow (ie mavis, VEP, delly)
            child_workflow = get_children_workflows(project_name, workflow_id)
            child_workflow = filter_out_QC_workflows(project_name, child_workflow)
            if child_workflow:
                for i in child_workflow:
                    for j in child_workflow[i]:
                        w = get_workflow_info(project_name, j)
                        if w:
                            D[libraries][w[libraries]['workflow']] = w[libraries]
    
    # add parent workflows
    for libraries in D:
        for workflow in D[libraries]:
            parent_workflow = get_parent_workflows(project_name, D[libraries][workflow]['workflow_id'])
            if 'parent' not in D[libraries][workflow]:
                D[libraries][workflow]['parent'] = []
            if parent_workflow not in D[libraries][workflow]['parent']:
                D[libraries][workflow]['parent'].append(parent_workflow)
        
   # get the files for each workflow
    for libraries in D:
        for workflow in D[libraries]:
            files, creation_date = get_workflow_files(project_name, D[libraries][workflow]['workflow_id'])
            D[libraries][workflow]['files'] = files
            # convert epoch time to standard time
            if creation_date:
                creation_date = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(creation_date)))
            D[libraries][workflow]['creation_date'] = creation_date
    
    return D



def get_samples(project_name):
    '''
    
    
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT library, sample, ext_id, group_id, group_id_description, library_type, tissue_origin, tissue_type FROM Libraries WHERE project_id = '{0}'".format(project_name)).fetchall()
    
    data = list(set(data))
    
    return data
    

def get_cases(project_name):
    '''
    
    
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT case_id, donor_id, species, sex, created_date, modified_date, miso FROM Samples WHERE project_id = '{0}'".format(project_name)).fetchall()
    
    # data = sorted([(i['case_id'], i) for i in data])
    #data = [i[1] for i in data]
    
    data = [dict(i) for i in data]
    
    
    return data




def get_last_sequencing(project_name):
    '''
    
    
    '''
    
    conn = connect_to_db()
    sequencing = conn.execute("SELECT Files.creation_date FROM Files WHERE Files.project_id = '{0}' \
                         AND LOWER(Files.workflow) in ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq');".format(project_name)).fetchall()
    conn.close()
    # get the most recent creation date of fastq generating workflows
    if sequencing:
        seq_date = sorted(list(set([i['creation_date'] for i in sequencing])))[-1]
        # convert to readable time
        seq_date = time.strftime('%Y-%m-%d', time.localtime(seq_date))
    else:
        seq_date = 'NA'
    return seq_date


def get_sample_counts(project_name, case):
    '''
    
    
    '''
    
    conn = connect_to_db()
    
    data = conn.execute("SELECT library, sample, tissue_type, group_id FROM Libraries WHERE project_id = '{0}' and sample = '{1}';".format(project_name, case)).fetchall()
    conn.close()

    libraries = len(set([dict(i)['library'] for i in data]))
    samples = [(dict(i)['sample'] + '_' + dict(i)['group_id'], dict(i)['tissue_type']) for i in data]
    normals = len(set([i[0] + '_' + i[1] for i in samples if i[1] == 'R']))    
    tumors = len(set([i[0] + '_' + i[1] for i in samples if i[1] != 'R']))
        
    return normals, tumors, libraries
    



# map pipelines to views
routes = {'Whole Genome': 'whole_genome_sequencing'}

    


app = Flask(__name__)


@app.template_filter()
def find_workflow_id(generic_name, bmpp_children_workflows, library):
    '''
    (str, str, dict, str) -> str
    
    Returns the workflow id of a workflow that has generic name as substring and
    NA if no workflow has generic name as substring.
            
    Parameters
    ----------
    - generic_name (str): Generic workflow name, may be substring of workflow name in bmpp_children_workflows
    - bmpp_children_workflows (dict): Dictionary with downstream workflow information
    - library (str): Libraries of interest
    '''
    
    # make a list of downstream workflows for that bmpp run
    L = list(bmpp_children_workflows[library].keys())
    # create a same size list of generic name workflows
    workflows = [generic_name] * len(L)
    
    # define function to identify generic workflow as subtring of workflow name
    is_workflow = lambda x,y: x.split('_')[0].lower() == y.lower()
    
    # check if generic workflow is substring of bmpp children workflows 
    found = list(map(is_workflow, L, workflows))
    if any(found):
        return bmpp_children_workflows[library][L[found.index(True)]]['workflow_id']
    else:
        return 'NA'
    


@app.route('/')
def index():
    
    # connect to db and extract project info
    conn = connect_to_db()
    projects = conn.execute('SELECT * FROM Projects').fetchall()
    conn.close()
    
    projects = sorted([(i['project_id'], i) for i in projects])
    projects = [i[1] for i in projects]
    
    return render_template('index.html', projects=projects)

@app.route('/<project_name>')
def project_page(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
    # get case information
    cases = get_cases(project_name)
    
    for i in range(len(cases)):
        normals, tumors, libraries = get_sample_counts(project_name, cases[i]['case_id'])
        cases[i]['normals'], cases[i]['tumors'], cases[i]['libraries'] = normals, tumors, libraries        
        
    # get the date of the last sequencing data
    seq_date = get_last_sequencing(project['project_id'])
    
    return render_template('project.html', routes=routes, project=project, pipelines=pipelines, cases=cases, seq_date=seq_date)

@app.route('/<project_name>/sequencing')
def sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
    
    # get sequences    
    conn = connect_to_db()
    files = conn.execute("SELECT Files.file, Files.workflow, Files.version, Files.wfrun_id, Files.attributes, \
                         FilesQC.status, FilesQC.ticket, Workflow_Inputs.run, Workflow_Inputs.lane, Workflow_Inputs.platform, \
                         Libraries.library, Libraries.sample, Libraries.ext_id, Libraries.group_id \
                         from Files JOIN FilesQC JOIN Workflow_Inputs JOIN Libraries WHERE Files.project_id = '{0}' \
                         AND FilesQC.project_id = '{0}' AND FilesQC.file_swid = Files.file_swid \
                         AND Workflow_Inputs.project_id = '{0}' AND Libraries.project_id = '{0}' \
                         AND Files.wfrun_id = Workflow_Inputs.wfrun_id AND Workflow_Inputs.library = Libraries.library \
                         AND LOWER(Files.workflow) in ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq');".format(project_name)).fetchall()
    conn.close()
    # find and group pairs of fastqs
    sequences = group_sequences(files)

    return render_template('sequencing.html', routes=routes, project=project, sequences=sequences, pipelines=pipelines)



@app.route('/<project_name>/whole_genome_sequencing')
def whole_genome_sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
        
    conn = connect_to_db()
    #extract sample, library and workflow information for call ready workflow bamMergePreprocessing
    data = conn.execute("SELECT Files.creation_date, Libraries.library, Libraries.sample, \
                          Libraries.ext_id, Libraries.group_id, Libraries.tissue_type, \
                          Libraries.tissue_origin, Workflow_Inputs.run, Workflow_Inputs.lane, \
                          Workflow_Inputs.platform, Workflows.wf, Workflows.wfv, Workflows.wfrun_id \
                          from Files JOIN Libraries JOIN Workflow_Inputs JOIN Workflows \
                          WHERE Files.project_id = '{0}' AND Libraries.project_id = '{0}' \
                          AND Workflow_Inputs.project_id = '{0}' AND Workflows.project_id = '{0}' \
                          AND Files.wfrun_id = Workflow_Inputs.wfrun_id  AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id AND Workflow_Inputs.library = Libraries.library \
                          AND LOWER(SUBSTR(Workflows.wf, 1, 21)) = 'bammergepreprocessing';".format(project_name)).fetchall()
    conn.close()

    # get samples and libraries for the most recent bmpp run for each case in project
    cases = get_call_ready_cases(data)
    
    samples = sorted(list(cases.keys()))

   
    return render_template('Whole_Genome_Sequencing.html', routes = routes, project=project, samples=samples, cases=cases, pipelines=pipelines)



@app.route('/<project_name>/whole_genome_sequencing/<case>')
def wgs_case(project_name, case):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
        
    # extract sample, library and workflow information for call ready workflow bamMergePreprocessing
    conn = connect_to_db()
    data = conn.execute("SELECT Files.creation_date, Files.file, Libraries.library, Libraries.sample, \
                         Libraries.ext_id, Libraries.group_id, Libraries.tissue_type, \
                         Libraries.tissue_origin, Workflow_Inputs.run, Workflow_Inputs.lane, \
                         Workflow_Inputs.platform, Workflows.wf, Workflows.wfv, Workflows.wfrun_id \
                         from Files JOIN Libraries JOIN Workflow_Inputs JOIN Workflows \
                         WHERE Files.project_id = '{0}' AND Libraries.project_id = '{0}' \
                         AND Workflow_Inputs.project_id = '{0}' AND Workflows.project_id = '{0}' \
                         AND Files.wfrun_id = Workflow_Inputs.wfrun_id  AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id AND Workflow_Inputs.library = Libraries.library \
                         AND LOWER(SUBSTR(Workflows.wf, 1, 21)) = 'bammergepreprocessing' AND Libraries.sample ='{1}'".format(project_name, case)).fetchall()
    conn.close()

    # get sample, library and file info for for the most recent bmpp run for case in project
    bmpp_id, bmpp_files, bmpp_info = get_bmpp_files(data)
    
    # get QC status of bmpp input fastq files
    fastq_status = bmpp_input_raw_seq_status(project_name, bmpp_id)

    # get the bmpp downstream workflows
    bmpp_children_workflows = get_bmpp_downstream_workflows(project_name, bmpp_id)
    
    for libraries in bmpp_children_workflows:
        for workflow in bmpp_children_workflows[libraries]:
            if bmpp_children_workflows[libraries][workflow]['files']:
                files = [os.path.dirname(bmpp_children_workflows[libraries][workflow]['files'][0])] + sorted(map(lambda x: os.path.basename(x), bmpp_children_workflows[libraries][workflow]['files']))
                bmpp_children_workflows[libraries][workflow]['files'] = files
    
    return render_template('WGS_case.html', routes = routes, fastq_status=fastq_status,
                            bmpp_info=bmpp_info, bmpp_id=bmpp_id, bmpp_files=bmpp_files,
                            sample_case=case, project=project, pipelines=pipelines,
                            bmpp_children_workflows=bmpp_children_workflows, case=case)



@app.route('/download_bmpp_data/<project_name>/<case>/<bmpp_id>')
def get_bmpp_data(project_name, case, bmpp_id):
    '''
    
    
    '''
    
    # get bmpp downstream workflow info
    bmpp_children_workflows = get_bmpp_downstream_workflows(project_name, bmpp_id)
    
    # format bmpp downstream workflow info for DARE
    D = {}
    for i in bmpp_children_workflows:
        for workflow in bmpp_children_workflows[i]:
            sample = bmpp_children_workflows[i][workflow]['sample'].replace(';', '.')
            workflow_id = bmpp_children_workflows[i][workflow]['workflow_id']
            version = bmpp_children_workflows[i][workflow]['version']    
            if sample not in D:
                D[sample] = {}
            D[sample][workflow] = {"workflow_id": workflow_id, "workflow_version": version}
    
    # add bmpp workflow info
    conn = connect_to_db()
    data = conn.execute("SELECT Workflows.wfrun_id, Workflows.wfv, Workflows.wf FROM Workflows \
                        WHERE Workflows.project_id = '{0}' AND Workflows.wfrun_id = '{1}';".format(project_name, bmpp_id)).fetchall()
    data = list(set(data))
    d = dict(data[0])
    conn.close()
    
    for sample in D:
        D[sample][d['wf']] = {'workflow_id': d['wfrun_id'], 'workflow_version': d['wfv']}
       
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(D),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}_WGS_{1}_{2}.json".format(project_name, case, bmpp_id)})





@app.route('/download_cases/<project_name>')
def download_cases_table(project_name):
    '''
    
    
    '''
    
    # get case information
    cases = get_cases(project_name)
    
    D = {}
    for i in cases:
        i = dict(i)
        assert i['case_id'] not in D
        D[i['case_id']] = i
    
    data = pd.DataFrame(D.values())
    data.to_excel('{0}_cases.xlsx'.format(project_name), index=False)
   
    return send_file("{0}_cases.xlsx".format(project_name), as_attachment=True)



@app.route('/download_identifiers/<project_name>')
def download_identifiers_table(project_name):
    '''
    
    
    '''
    
    conn = connect_to_db()
    identifiers = conn.execute("SELECT library, sample, ext_id, group_id, group_id_description, \
                               library_type, tissue_origin, tissue_type FROM Libraries \
                               WHERE Libraries.project_id = '{0}';".format(project_name)).fetchall()
    conn.close()
    
    identifiers = list(set(identifiers))
    
    D ={}
    for i in range(len(identifiers)):
        D[i] = dict(identifiers[i])
    data = pd.DataFrame(D.values())
    
    outputfile = '{0}_libraries.xlsx'.format(project_name)
    data.to_excel(outputfile, index=False)
   
    return send_file(outputfile, as_attachment=True)


if __name__ == "__main__":
    app.run(host='0.0.0.0')
    