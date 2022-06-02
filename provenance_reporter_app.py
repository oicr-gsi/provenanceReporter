# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:32:40 2022

@author: rjovelin
"""

import sqlite3
import json
from flask import Flask, render_template, request, url_for, flash, redirect
from werkzeug.exceptions import abort
import requests
import gzip



def get_provenance_data(provenance):
    '''
    (str) -> list
    
    Returns a list of dictionary with lims information for each library
    
    Parameters
    ----------
    - provenance (str): URL of the pinery provenance API 
    '''
    
    response = requests.get(provenance)
    if response.ok:
        L = response.json()
    else:
        L = []
    
    return L



def collect_info(data, names, keys):
    '''
    (dict, list, list) -> dict

    Returns a dictionary with keys by extracting specific information using names from data
    Note, the value of any name in data may be a string or single-element-list

    Parameters
    ----------
    - data (dict): Dictionary with data to be extracted
    - names (list): List of keys of interest in data
    - keys (list): List of renamed keys in output dictionary
    '''
    
    d = {}
    for j in range(len(names)):
        if names[j] in data:
            if type(data[names[j]]) == list:
                d[keys[j]] = data[names[j]][0]         
            else:
                d[keys[j]] = data[names[j]]
        else:
            d[keys[j]] = ''

    return d
    


def extract_lims(sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
    '''
    (str) -> dict
    
    Returns a dictionary with lims information for each library of each project

    Parameters
    ----------
    - sample_provenance (str): URL of the pinery sample_provenance API
    ''' 

    # get sample info from pinery
    L = get_provenance_data(sample_provenance)
    
    # store lims information for each library of each project
    # {project: {sample: {library: info, sequencing: [sequencing info]}}}    
    D = {}

    for i in L:
        project = i['studyTitle']
        sample = i['rootSampleName']  
        if project not in D:
            D[project] = {}
        if sample not in D[project]:
            D[project][sample] = {}
    
        # collect sample information
        d = collect_info(i['sampleAttributes'], ['geo_tissue_type', 'geo_external_name', 'geo_tissue_origin',
                 'geo_library_source_template_type', 'geo_prep_kit',
                 'geo_tissue_preparation', 'geo_receive_date'], ['tissue_type', 'ext_id', 'tissue_origin', 'library_type', 
                         'prep', 'tissue_prep', 'sample_received_date'])
        # add library name
        library = i['sampleName']
        d['library'] = library
                
        # store sample information
        if library not in D[project][sample]:
            D[project][sample][library] = {}
            D[project][sample][library]['library_info'] = d
            D[project][sample][library]['sequencing'] = []
        
        # update sample information for each library if some fields are missing
        for k in d:
            if D[project][sample][library]['library_info'][k] == '':
                D[project][sample][library]['library_info'][k] = d[k]    
        
        # collect sequencing information for each library        
        
        s = collect_info(i, ['sequencerRunName', 'laneNumber', 'sampleProvenanceId', 'iusTag',
                             'sequencerRunPlatformModel'], ['run', 'lane', 'limskey', 'barcode', 'platform'])
        D[project][sample][library]['sequencing'].append(s)
        
    return D    
        


def extract_qc_status_from_nabu(project, nabu_api = 'http://gsi-dcc.oicr.on.ca:3000'):
    '''
    (str, str) -> dict

    Returns a dictionary with QC status of all files in project    
        
    Parameters
    ----------
    - project (str): Name of project of interest
    - nabu_api (str): URL of the nabu API
    '''
        
    response = requests.get(nabu_api + '/fileqcs?project={0}'.format(project))
    
    if response.ok:
        L = response.json()['fileqcs']
    else:
        L = []
    
    D = {}
    
    if L:
        for i in L:
            filepath = i['filepath']
            d = {}
            qc = collect_info(i, ['skip', 'user', 'date', 'qcstatus', 'ref', 'stalestatus'], ['skip', 'user', 'date', 'status', 'ref', 'fresh']) 
            d['qc'] = qc
            d['fid'] = 'f.' + str(i['fileswid'])
            d['filepath'] = filepath
            assert filepath not in D
            D[filepath] = d 
    
    return D


def is_gzipped(file):
    '''
    (str) -> bool

    Return True if file is gzipped

    Parameters
    ----------
    - file (str): Path to file
    '''
    
    # open file in rb mode
    infile = open(file, 'rb')
    header = infile.readline()
    infile.close()
    if header.startswith(b'\x1f\x8b\x08'):
        return True
    else:
        return False


def open_fpr(fpr):
    '''
    (str) -> _io.TextIOWrapper
    
    Returns a file open for reading
    
    Parameters
    ----------
    - fpr (str): Path to File Provenance Report file
    '''

    # open provenance for reading. allow gzipped file or not
    if is_gzipped(fpr):
        infile = gzip.open(fpr, 'rt', errors='ignore')
    else:
        infile = open(fpr)
    return infile



def get_FPR_records(project, fpr):
    '''
    (str, str) -> list
    
    Returns a list with all the records for project from the File Provenance Report.
        
    Parameters
    ----------
    - project (str): Name of project of interest
    - fpr (str): Path to File Provenance Report file
    '''
        
    # get the records for a single project
    records = []
    # open provenance for reading. allow gzipped file or not
    infile = open_fpr(fpr)
    for line in infile:
        if project in line:
            line = line.rstrip().split('\t')
            if line[1] == project:
                records.append(line)
    infile.close()
    return records


def collect_file_info_from_fpr(project, fpr):
    '''
    (str, str) -> dict

    Returns a dictionary with file information extracted from File Provenance Report
    for all files in project     
    
    Parameters
    ----------
    - project (str): Name of project of interest
    - fpr (str): Path to File Provenance Report file
    '''

    records = get_FPR_records(project, fpr)
    
    D = {}
    
    
    for i in records:
        # get file path
        file = i[46]
        # get md5sums
        md5sum = i[47]
        # get file_swid
        file_swid = i[44]
        # get workflow
        workflow = i[30]
        # get workflow version
        version = i[31]        
        # get workflow run accession
        workflow_run = i[36]
        # collect file info
        D[file] = {'md5sum': md5sum, 'fid': file_swid, 'wf': workflow, 'wfv': version,
                   'wfrunid': workflow_run, 'file': file}

    return D


def add_file_info_to_qc(qc_info, file_info):
    '''
    (dict, dict) -> None
    
    Modifies dictionary qc_info in place with file info (workflow, version and md5sum)
    
    Parameters
    ----------
    - qc_info (dict): Dictionary with qc_info generated by function extract_qc_status_from_nabu
    - file_info (dict): Dictionary with file info generated collect_file_info_from_fpr
    '''

    for file in file_info:
        assert file in qc_info
        for i in ['md5sum', 'wfrunid', 'wfv', 'wf']:
            qc_info[file][i] = file_info[file][i]
        


def extract_file_info(project, fpr, nabu_api = 'http://gsi-dcc.oicr.on.ca:3000'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with QC and file information for all files in project
    
    Parameters
    ----------
    - project (str): name of project of interest
    - fpr (str): Path to the File Provenance Reporter 
    - nabu_api (str): URL of the Nabu api
    '''
    # extract file information from File Provenance Report
    file_info = collect_file_info_from_fpr(project, fpr)    
    # extract QC information from Nabu
    qc_info = extract_qc_status_from_nabu(project, nabu_api)
    # add file info to qc info
    add_file_info_to_qc(qc_info, file_info)
    
    return qc_info


def get_parent_workflows(project, fpr):
    '''
    (str, str) -> (dict, dict, dict)

    Returns a tuple with dictionaries with worklow information, parent file ids for each workflow,
    and workflow id for each file and project of interest
    
    Parameters
    ----------
    - project (str): Project name of interest
    - fpr (str): Path to the File Provenance Report
    '''

    F, P, W = {}, {}, {}

    records = get_FPR_records(project, fpr)
    
    for i in records:
        # get project name
        project = i[1]
        # get workflow, workflow version and workflow run accession
        workflow, workflow_version, workflow_run = i[30], i[31], i[36]
        # get file path, md5sum and file id
        file, md5sum, file_swid = i[46], i[47], i[44]
        input_files = i[38]
        if input_files:
            input_files = sorted(input_files.split(';'))
        else:
            input_files = []
        
        if project not in P:
            P[project] = {}
        
        if workflow_run in P[project]:
            assert P[project][workflow_run] == input_files
        else:
            P[project][workflow_run] = input_files
                
        if project not in W:
            W[project] = {}
        if workflow_run in W[project]:
            assert W[project][workflow_run] == {'wfrun_id': 'wf.' + workflow_run, 'wfv': workflow_version, 'wf': workflow}
        else:
            W[project][workflow_run] = {'wfrun_id': 'wf.' + workflow_run, 'wfv': workflow_version, 'wf': workflow}
        
        if project not in F:
            F[project] = {}
        if file_swid in F[project]:
            assert F[project][file_swid] == workflow_run
        else:
            F[project][file_swid] = workflow_run

    return W, P, F        
        


def identify_parent_children_workflows(W, P, F):
    '''
    (dict, dict, dict) -> (dict, dict)     
    
    Returns a tuple with dictionaries of parents: children workflows, and children: parents
    workflows relationsips for a given project
        
    Parameters
    ----------
    - W (dict): Dictionary with workflow information for each workflow run id
    - P (dict): Input file ids for each workflow run id
    - F (dict): Map of file id and workflow id
    '''
    
    parents, children = {}, {}
    
    for project in P:
        if project not in parents:
            parents[project] = {}
        for workflow in P[project]:
            parent_workflows = sorted(list(set([F[project][i] for i in P[project][workflow]])))
            parents[project][workflow] = parent_workflows
    
    for project in parents:
        if project not in children:
            children[project] = {}
        for workflow in parents[project]:
            for parent_workflow in parents[project][workflow]:
                if parent_workflow not in children[project]:
                    children[project][parent_workflow] = [workflow]
                else:
                    children[project][parent_workflow].append(workflow)
    
    for project in children:
        for workflow in children[project]:
            children[project][workflow] = sorted(list(set(children[project][workflow])))
    
    return parents, children
        



def extract_workflow_info(project, fpr):
    '''
    (str, str) -> dict
    
    Returns a dictionary with information about all workflows in project
    
    Parameters
    ----------
    - project (str): Project name of interest
    - fpr (str): Path to the File Provenance Report
    '''


    # get information about workflows, worklow input files and maps of file and workflows
    W, P, F = get_parent_workflows(project, fpr)
    # identify parents, children workflow relationships
    parents, children = identify_parent_children_workflows(W, P, F)

    D = {}

    records = get_FPR_records(project, fpr)
    
    for i in records:
        
        #if i[30] == 'bcl2fastq':
        if i[30]:
            # get project name
            project = i[1]
            # get sample name
            sample = i[7]
            # get workflow and workflow version
            workflow, workflow_version = i[30], i[31]
            # get workflow run accession
            workflow_run = i[36]
        
            if project not in D:
                D[project] = {}
            if workflow not in D[project]:
                D[project][workflow] = {}
            if sample not in D[project][workflow]:
                D[project][workflow][sample] = {}
            if workflow_run not in D[project][workflow][sample]:
                D[project][workflow][sample][workflow_run] = {}
            
        
            # get file path, md5sum and file id
            file, md5sum, file_swid = i[46], i[47], i[44]
            # get lane and run
            run, lane = i[18], i[24]            
            # get library and limskey
            library, limskey  = i[13], i[56]
            # get file attributes
            file_attributes = i[45]
            if file_attributes:
                file_attributes = file_attributes.split(';')
                file_attributes = {k.split('=')[0]: k.split('=')[1] for k in file_attributes}
                if 'read_count' in file_attributes and 'read_number' in file_attributes:
                    read_count = file_attributes['read_count']
                    read_number = file_attributes['read_number']
                else:
                    read_count, read_number = '', ''
                    
            #### check other workflows for information captured. read count and read number may not be needed in other workflows
            
            D[project][workflow][sample][workflow_run]['libraries'] = [{'lane': lane, 'run': run, 'limskey': limskey, 'id': sample, 'lib': library}]
            if 'files' not in D[project][workflow][sample][workflow_run]:
                D[project][workflow][sample][workflow_run]['files'] = []
            file_info = {'md5sum': md5sum, 'wfrunid': 'wf.' + workflow_run, 'path': file,
             'fid': 'f.' + file_swid, 'wfv': workflow_version, 'wf': workflow}
            if file_info not in D[project][workflow][sample][workflow_run]['files']:
                D[project][workflow][sample][workflow_run]['files'].append(file_info)
                
            D[project][workflow][sample][workflow_run]['wfrunid'] = 'wf.' + workflow_run    
            
            ####### what is wfinput_string
            D[project][workflow][sample][workflow_run]['info'] = {'read_number': read_number,
                                                                  'wfv': workflow_version,
                                                                  'wf': workflow,
                                                                  'read_count': read_count,
                                                                  'wfinput_string': ''}
            # get information about children workflows
            children_workflows_info = []
            if workflow_run in children[project]:
                for k in children[project][workflow_run]:
                    children_workflows_info.append(W[project][k])
            D[project][workflow][sample][workflow_run]['children'] = {'workflows': children_workflows_info}
            
            # get information about parent workflows
            parent_workflows_info = []
            if workflow_run in parents[project]:
                for k in parents[project][workflow_run]:
                    parent_workflows_info.append(W[project][k])
            D[project][workflow][sample][workflow_run]['parents'] = {'workflows': parent_workflows_info}

    return D            


def extract_projects(fpr):
    '''
    (str) -> list

    Returns a list of projects present in File Provenance Report

    Parameters
    ----------
    - fpr (str): Path to the File Provenance Report    
    '''
    
    infile = open_fpr(fpr)
    projects = set()
    for line in infile:
        projects.add(line.rstrip().split('\t')[1])
    infile.close()
    projects = sorted(list(projects))
    return projects



def validate_projects(pinery_projects, fpr_projects):
    '''
    (list, list) -> list
    
    Returns a list of projects present both in File Provenance Report and in Pinery
    
    Parameters
    ----------
    - pinery_projects (list): List of projects defined in Pinery
    - fpr_projects (list): List of projects recorded in File Provenance Report
    '''

    valid_projects = set(pinery_projects).intersection(set(fpr_projects))
    valid_projects = sorted(list(valid_projects))
    return valid_projects

    



####################################

def extract_project_info(project_provenance = 'http://pinery.gsi.oicr.on.ca/sample/projects'):
    '''
    (str) -> list
    
    Returns a list of dictionary with project information pulled down from the
    project_provenance Pinary API

    Parameters
    ----------
    - project_provenance (str): Pinery API, http://pinery.gsi.oicr.on.ca/sample/projects
    '''
    
    response = requests.get(project_provenance)
    if response.ok:
        L = response.json()
    else:
        L = []
    
    D = {}
    
    if L:
        for i in L:
            name = i['name']
            assert name not in D
            D[name] = {'name': name}
            for j in ['pipeline', 'description', 'active', 'contact_name', 'contact_email']:
                if j in i:
                    D[name][j] = i[j] 
                else:
                    D[name][j] = ''
                if j == 'active' and j in i:
                    if i[j]:
                        D[name][j] = 'Active'
                    else:
                        D[name][j] = 'Completed'
    
    S = [D[i] for i in sorted(list(D.keys()))]
    return S                            
    
    
    


def collect_sequence_info(data):
    '''    
    
    '''
    
    # create a dict to store sequence info for each case
    # {case: {library: {run: lane: [[file_path, file_id, read_count]]}}}
    D = {}
    
    
    for i in data['pipeline']['workflows']:
        if i['info']['wf'].lower() in ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport']:
            workflow = i['info']['wf'] + '_' + i['info']['wfv']
            assert len(i['libraries']) == 1
            case = i['libraries'][0]['id']
            lib = i['libraries'][0]['lib']
            run = i['libraries'][0]['run']
            lane = i['libraries'][0]['lane']
            if 'read_count' in i['info']:
                read_count = i['info']['read_count']
            else:
                read_count = ''
            assert len(i['files']) == 2
            if case not in D:
                D[case] = {}
            if lib not in D[case]:
                D[case][lib] = {}
            if run not in D[case][lib]:
                D[case][lib][run] = {}
            if lane not in D[case][lib][run]:
                D[case][lib][run][lane] = []
            for j in i['files']:
                D[case][lib][run][lane].append([j['path'], j['fid'].split('.')[1], read_count, workflow])
            D[case][lib][run][lane].sort()
            assert len(D[case][lib][run][lane]) == 2
    return D



def collect_lims_info(data):
    '''
    
    
    
    '''
    
    # create dict {case: {library: group_id}}
    D = {}
    
    for i in data['lims']['samples']:
        case = i['id']
        for j in i['libraries']:
            library = j['lib']
            group_id = j['group_id']
            ext_id = j['ext_id']
            if case not in D:
                D[case] = {}
            D[case]['ext_id'] = ext_id
            assert library not in D[case]
            D[case][library] = group_id
    return D    


def collect_release_status(data):
    '''
    
    
    '''
    
    D = {}
    
    for i in data['pipeline']['files']:
        path = i['path']
        file_id = i['fid']
        status = i['qc']['status']
        assert path not in D
        D[path] = [file_id.split('.')[1], status]
    
    return D



def add_lims_info_to_sequence_data(data):
    '''
    
    
    '''
    
    S = collect_sequence_info(data)
    L = collect_lims_info(data)
    R = collect_release_status(data)
    
    {}
    D = {}
    
    
    # {case: {library: {run: lane: [[file_path, file_id, read_count, workflow]]}}}
    # create dict {case: {sample: {library: group_id, 'ext_id': exit_id}}}
    
    for case in S:
        sample = case + '_' + L[case]['ext_id']   
        for library in S[case]:
            lib = library + '_' + L[case][library]
            for run_id in S[case][library]:
                 for lane in S[case][library][run_id]:
                     run = run_id + '_' + lane
                     files = []
                     file_id = []
                     workflow = []
                     read_count = []
                     release_status = []
                     for file in S[case][library][run_id][lane]:
                         files.append(file[0])
                         file_id.append(file[1])
                         read_count = file[2]
                         workflow = file[3]
                         release_status = R[file[0]][-1]
                         assert file[1] == R[file[0]][0]
                     assert len(files) == 2
                     assert len(file_id) == 2
                     
                     
                     if case not in D:
                         D[case] = {}
                     if sample not in D[case]:
                         D[case][sample] = {}
                     if lib not in D[case][sample]:
                         D[case][sample][lib] = {}
                     
                     assert run not in D[case][sample][lib]
                     D[case][sample][lib][run] = {'files': sorted(files), 'file_id': file_id, 'release': release_status, 'read_count': read_count, 'workflow': workflow}
                     
    return D            
            

infile = open('TGL01MOH.json')
data = json.load(infile)
infile.close()

S = collect_sequence_info(data)
L = collect_lims_info(data)
R = collect_release_status(data)
D = add_lims_info_to_sequence_data(data)





app = Flask(__name__)

@app.route('/')
def index():
    projects = extract_project_info()
    return render_template('index.html', projects=projects)

@app.route('/<project_name>')
def project(project_name):
    projects = extract_project_info()
    project = [i for i in projects if i['name'] == project_name][0]
    return render_template('project.html', project=project)

@app.route('/<project_name>/sequencing')
def sequencing(project_name):
    projects = extract_project_info()
    project = [i for i in projects if i['name'] == project_name][0]
    
    
    infile = open('TGL01MOH.json')
    data = json.load(infile)
    infile.close()

    sequences = add_lims_info_to_sequence_data(data)
    
    
    
    return render_template('sequencing.html', project=project, sequences=sequences)





