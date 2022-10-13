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
import os
import time


# def get_provenance_data(provenance):
#     '''
#     (str) -> list
    
#     Returns a list of dictionary with lims information for each library
    
#     Parameters
#     ----------
#     - provenance (str): URL of the pinery provenance API 
#     '''
    
#     response = requests.get(provenance)
#     if response.ok:
#         L = response.json()
#     else:
#         L = []
    
#     return L



# def collect_info(data, names, keys):
#     '''
#     (dict, list, list) -> dict

#     Returns a dictionary with keys by extracting specific information using names from data
#     Note, the value of any name in data may be a string or single-element-list

#     Parameters
#     ----------
#     - data (dict): Dictionary with data to be extracted
#     - names (list): List of keys of interest in data
#     - keys (list): List of renamed keys in output dictionary
#     '''
    
#     d = {}
#     for j in range(len(names)):
#         if names[j] in data:
#             if type(data[names[j]]) == list:
#                 d[keys[j]] = data[names[j]][0]         
#             else:
#                 d[keys[j]] = data[names[j]]
#         else:
#             d[keys[j]] = ''

#     return d
    


# def extract_lims(sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with lims information for each library of each project

#     Parameters
#     ----------
#     - sample_provenance (str): URL of the pinery sample_provenance API
#     ''' 

#     # get sample info from pinery
#     L = get_provenance_data(sample_provenance)
    
#     # store lims information for each library of each project
#     # {project: {sample: {library: info, sequencing: [sequencing info]}}}    
#     D = {}

#     for i in L:
#         project = i['studyTitle']
#         sample = i['rootSampleName']  
#         if project not in D:
#             D[project] = {}
#         if sample not in D[project]:
#             D[project][sample] = {}
    
#         # collect sample information
#         d = collect_info(i['sampleAttributes'], ['geo_tissue_type', 'geo_external_name', 'geo_tissue_origin',
#                  'geo_library_source_template_type', 'geo_prep_kit',
#                  'geo_tissue_preparation', 'geo_receive_date', 'geo_group_id', 'geo_group_id_description'], ['tissue_type', 'ext_id', 'tissue_origin', 'library_type', 
#                          'prep', 'tissue_prep', 'sample_received_date', 'group_id', 'group_id_description'])
#         # add library name
#         library = i['sampleName']
#         d['library'] = library
                
#         # store sample information
#         if library not in D[project][sample]:
#             D[project][sample][library] = {}
#             D[project][sample][library]['library_info'] = d
#             D[project][sample][library]['sequencing'] = []
        
#         # update sample information for each library if some fields are missing
#         for k in d:
#             if D[project][sample][library]['library_info'][k] == '':
#                 D[project][sample][library]['library_info'][k] = d[k]    
        
#         # collect sequencing information for each library        
        
#         s = collect_info(i, ['sequencerRunName', 'laneNumber', 'sampleProvenanceId', 'iusTag',
#                              'sequencerRunPlatformModel'], ['run', 'lane', 'limskey', 'barcode', 'platform'])
#         D[project][sample][library]['sequencing'].append(s)
        
#     return D    
        


# def extract_qc_status_from_nabu(project, nabu_api = 'http://gsi-dcc.oicr.on.ca:3000'):
#     '''
#     (str, str) -> dict

#     Returns a dictionary with QC status of all files in project    
        
#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - nabu_api (str): URL of the nabu API
#     '''
        
#     response = requests.get(nabu_api + '/fileqcs?project={0}'.format(project))
    
#     if response.ok:
#         L = response.json()['fileqcs']
#     else:
#         L = []
    
#     D = {project: {}}
    
#     if L:
#         for i in L:
#             filepath = i['filepath']
#             d = {}
#             qc = collect_info(i, ['skip', 'user', 'date', 'qcstatus', 'ref', 'stalestatus'], ['skip', 'user', 'date', 'status', 'ref', 'fresh']) 
#             d['qc'] = qc
#             d['fid'] = 'f.' + str(i['fileswid'])
#             d['filepath'] = filepath
#             #assert filepath not in D[project]
#             D[project][filepath] = d 
    
#     return D


# def is_gzipped(file):
#     '''
#     (str) -> bool

#     Return True if file is gzipped

#     Parameters
#     ----------
#     - file (str): Path to file
#     '''
    
#     # open file in rb mode
#     infile = open(file, 'rb')
#     header = infile.readline()
#     infile.close()
#     if header.startswith(b'\x1f\x8b\x08'):
#         return True
#     else:
#         return False


# def open_fpr(fpr):
#     '''
#     (str) -> _io.TextIOWrapper
    
#     Returns a file open for reading
    
#     Parameters
#     ----------
#     - fpr (str): Path to File Provenance Report file
#     '''

#     # open provenance for reading. allow gzipped file or not
#     if is_gzipped(fpr):
#         infile = gzip.open(fpr, 'rt', errors='ignore')
#     else:
#         infile = open(fpr)
#     return infile



# def get_FPR_records(project, fpr):
#     '''
#     (str, str) -> list
    
#     Returns a list with all the records for project from the File Provenance Report.
        
#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - fpr (str): Path to File Provenance Report file
#     '''
        
#     # get the records for a single project
#     records = []
#     # open provenance for reading. allow gzipped file or not
#     infile = open_fpr(fpr)
#     for line in infile:
#         if project in line:
#             line = line.rstrip().split('\t')
#             if line[1] == project:
#                 records.append(line)
#     infile.close()
#     return records


# def collect_file_info_from_fpr(project, fpr):
#     '''
#     (str, str) -> dict

#     Returns a dictionary with file information extracted from File Provenance Report
#     for all files in project     
    
#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - fpr (str): Path to File Provenance Report file
#     '''

#     records = get_FPR_records(project, fpr)
    
#     D = {project: {}}
    
    
#     for i in records:
#         # get file path
#         file = i[46]
#         # get md5sums
#         md5sum = i[47]
#         # get file_swid
#         file_swid = i[44]
        
#         d = collect_info({k.split('=')[0]:k.split('=')[1] for k in i[17].split(';')},
#                          ['geo_library_source_template_type'], ['library_type'])
#         # get workflow
#         workflow = i[30]
#         # get workflow version
#         version = i[31]        
#         # get workflow run accession
#         workflow_run = i[36]
#         # collect file info
#         D[project][file] = {'md5sum': md5sum, 'fid': file_swid, 'wf': workflow, 'wfv': version,
#                    'wfrunid': workflow_run, 'file': file, 'library_type': d['library_type']}

#     return D


# def add_file_info_to_qc(qc_info, file_info):
#     '''
#     (dict, dict) -> None
    
#     Modifies dictionary qc_info in place with file info (workflow, version and md5sum)
    
#     Parameters
#     ----------
#     - qc_info (dict): Dictionary with qc_info generated by function extract_qc_status_from_nabu
#     - file_info (dict): Dictionary with file info generated collect_file_info_from_fpr
#     '''

#     for project in file_info:
#         for file in file_info[project]:
#             # udpate qc info with file info
#             if file not in qc_info[project]:
#                 qc_info[project][file] = {'qc': {'skip': '', 'user': '', 'date': '',
#                                           'status': '', 'ref': '', 'fresh': ''},
#                                           'fid': '', 'filepath': file}
#             for i in file_info[project][file]:
#                 qc_info[project][file][i] = file_info[project][file][i]
            

# def extract_file_info(project, fpr, nabu_api = 'http://gsi-dcc.oicr.on.ca:3000'):
#     '''
#     (str, str, str) -> dict
    
#     Returns a dictionary with QC and file information for all files in project
    
#     Parameters
#     ----------
#     - project (str): name of project of interest
#     - fpr (str): Path to the File Provenance Reporter 
#     - nabu_api (str): URL of the Nabu api
#     '''
#     # extract file information from File Provenance Report
#     file_info = collect_file_info_from_fpr(project, fpr)    
#     # extract QC information from Nabu
#     qc_info = extract_qc_status_from_nabu(project, nabu_api)
#     # add file info to qc info
#     add_file_info_to_qc(qc_info, file_info)
    
#     return qc_info


# def get_parent_workflows(project, fpr):
#     '''
#     (str, str) -> (dict, dict, dict)

#     Returns a tuple with dictionaries with worklow information, parent file ids for each workflow,
#     and workflow id for each file and project of interest
    
#     Parameters
#     ----------
#     - project (str): Project name of interest
#     - fpr (str): Path to the File Provenance Report
#     '''

#     F, P, W = {}, {}, {}

#     records = get_FPR_records(project, fpr)
    
#     for i in records:
#         # get project name
#         assert project == i[1]
#         # get workflow, workflow version and workflow run accession
#         workflow, workflow_version, workflow_run = i[30], i[31], i[36]
#         # get file path, md5sum and file id
#         file, md5sum, file_swid = i[46], i[47], i[44]
#         input_files = i[38]
#         if input_files:
#             input_files = sorted(input_files.split(';'))
#         else:
#             input_files = []
        
#         if project not in P:
#             P[project] = {}
      
#         if workflow_run in P[project]:
#             assert P[project][workflow_run] == input_files
#         else:
#             P[project][workflow_run] = input_files
                
#         if project not in W:
#             W[project] = {}
#         if workflow_run in W[project]:
#             assert W[project][workflow_run] == {'wfrun_id': 'wf.' + workflow_run, 'wfv': workflow_version, 'wf': workflow}
#         else:
#             W[project][workflow_run] = {'wfrun_id': 'wf.' + workflow_run, 'wfv': workflow_version, 'wf': workflow}
        
#         if project not in F:
#             F[project] = {}
#         if file_swid in F[project]:
#             assert F[project][file_swid] == workflow_run
#         else:
#             F[project][file_swid] = workflow_run

#     return W, P, F        
        


# def identify_parent_children_workflows(W, P, F):
#     '''
#     (dict, dict, dict) -> (dict, dict)     
    
#     Returns a tuple with dictionaries of parents: children workflows, and children: parents
#     workflows relationsips for a given project
        
#     Parameters
#     ----------
#     - W (dict): Dictionary with workflow information for each workflow run id
#     - P (dict): Input file ids for each workflow run id
#     - F (dict): Map of file id and workflow id
#     '''
    
#     parents, children = {}, {}
    
#     for project in P:
#         if project not in parents:
#             parents[project] = {}
#         for workflow in P[project]:
#             parent_workflows = sorted(list(set([F[project][i] for i in P[project][workflow] if i in F[project]])))
#             parents[project][workflow] = parent_workflows
    
#     for project in parents:
#         if project not in children:
#             children[project] = {}
#         for workflow in parents[project]:
#             for parent_workflow in parents[project][workflow]:
#                 if parent_workflow not in children[project]:
#                     children[project][parent_workflow] = [workflow]
#                 else:
#                     children[project][parent_workflow].append(workflow)
    
#     for project in children:
#         for workflow in children[project]:
#             children[project][workflow] = sorted(list(set(children[project][workflow])))
    
#     return parents, children
        



# def extract_workflow_info(project, fpr):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with information about all workflows in project
    
#     Parameters
#     ----------
#     - project (str): Project name of interest
#     - fpr (str): Path to the File Provenance Report
#     '''


#     # get information about workflows, worklow input files and maps of file and workflows
#     W, P, F = get_parent_workflows(project, fpr)
#     # identify parents, children workflow relationships
#     parents, children = identify_parent_children_workflows(W, P, F)

#     D = {}

#     records = get_FPR_records(project, fpr)
    
#     for i in records:
#         if project == i[1]:
#             # get sample name
#             sample = i[7]
#             # get workflow and workflow version
#             workflow, workflow_version = i[30], i[31]
#             # get workflow run accession
#             workflow_run = i[36]
        
#             if project not in D:
#                 D[project] = {}
#             if workflow not in D[project]:
#                 D[project][workflow] = {}
#             if sample not in D[project][workflow]:
#                 D[project][workflow][sample] = {}
#             if workflow_run not in D[project][workflow][sample]:
#                 D[project][workflow][sample][workflow_run] = {}
            
        
#             # get file path, md5sum and file id
#             file, md5sum, file_swid = i[46], i[47], i[44]
#             # get lane and run
#             run, lane = i[18], i[24]            
#             # get library and limskey
#             library, limskey  = i[13], i[56]
#             # get file attributes
#             file_attributes = i[45]
#             if file_attributes:
#                 file_attributes = file_attributes.split(';')
#                 file_attributes = {k.split('=')[0]: k.split('=')[1] for k in file_attributes}
#             else:
#                 file_attributes = {}
            
#             # get workflow info
#             info = i[37]
#             if info:
#                 info = info.split(';')
#                 info = {k.split('=')[0]: k.split('=')[1] for k in info if k.split('=')[0] not in ['cromwell-workflow-id', 'major_olive_version']}
#             else:
#                 info = {}                
        
#             # get library type
#             d = collect_info({k.split('=')[0]:k.split('=')[1] for k in i[17].split(';')},
#                              ['geo_library_source_template_type'], ['library_type'])
                    
#             if 'libraries' not in D[project][workflow][sample][workflow_run]:
#                 D[project][workflow][sample][workflow_run]['libraries'] = []
#             lib_info = {'lane': lane, 'run': run, 'limskey': limskey, 'id': sample, 'lib': library, 'library_type': d['library_type']}
#             if lib_info not in D[project][workflow][sample][workflow_run]['libraries']:
#                 D[project][workflow][sample][workflow_run]['libraries'].append(lib_info)
                        
#             if 'files' not in D[project][workflow][sample][workflow_run]:
#                 D[project][workflow][sample][workflow_run]['files'] = []
#             file_info = {'md5sum': md5sum, 'wfrunid': 'wf.' + workflow_run, 'path': file,
#              'fid': 'f.' + file_swid, 'wfv': workflow_version, 'wf': workflow}
#             if file_info not in D[project][workflow][sample][workflow_run]['files']:
#                 D[project][workflow][sample][workflow_run]['files'].append(file_info)
                
#             D[project][workflow][sample][workflow_run]['wfrunid'] = 'wf.' + workflow_run    

#             D[project][workflow][sample][workflow_run]['info'] = {}
#             for k in file_attributes:
#                 D[project][workflow][sample][workflow_run]['info'][k] = file_attributes[k]
#             for k in info:
#                 D[project][workflow][sample][workflow_run]['info'][k] = info[k]
#             D[project][workflow][sample][workflow_run]['info']['wfv'] = workflow_version
#             D[project][workflow][sample][workflow_run]['info']['wf'] = workflow
                        
#             # get information about children workflows
#             children_workflows_info = []
#             if workflow_run in children[project]:
#                 for k in children[project][workflow_run]:
#                     children_workflows_info.append(W[project][k])
#             D[project][workflow][sample][workflow_run]['children'] = {'workflows': children_workflows_info}
            
#             # get information about parent workflows
#             parent_workflows_info = []
#             if workflow_run in parents[project]:
#                 for k in parents[project][workflow_run]:
#                     parent_workflows_info.append(W[project][k])
#             D[project][workflow][sample][workflow_run]['parents'] = {'workflows': parent_workflows_info}

#     return D            


# def extract_project_data(project, fpr, nabu_api='http://gsi-dcc.oicr.on.ca:3000', sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
#     '''
#     (str, str, str, str) -> dict

#     Returns a dictionary with workflow, file and sample information for a given project    
        
#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - fpr (str): Path to the File Provenance Report
#     - nabu_api (str): URL address of the Nabu API
#     - sample_provenance (str): URL address of the Pinery sample provenance API
#     '''

#     # extract information about workflows
#     workflow_info = extract_workflow_info(project, fpr)
#     # extract qc and file information
#     file_info = extract_file_info(project, fpr, nabu_api)
#     # extract sample info
#     lims = extract_lims(sample_provenance)

#     data ={project: {'workflows': {}, 'files': [], 'lims': []}}

#     for workflow in workflow_info[project]:
#         #data[project]['workflows'].append({workflow: workflow_info[project][workflow]})
#         data[project]['workflows'][workflow] = workflow_info[project][workflow]
        
#     for file in file_info[project]:
#         data[project]['files'].append(file_info[project][file])
    
#     for sample in lims[project]:
#         data[project]['lims'].append({'id': sample, 'libraries': lims[project][sample]})
    
#     return data




# def write_json_to_file(outputfile, project, fpr, nabu_api='http://gsi-dcc.oicr.on.ca:3000', sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
    
    
#     data = extract_project_data(project, fpr, nabu_api, sample_provenance)
    
#     newfile = open(outputfile, 'w')
#     json.dump(data, newfile,indent=4) 
#     newfile.close()
    


# # fpr = '/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz'    
# # outputfolder = '/scratch2/groups/gsi/bis/rjovelin/provenance_reporter'

# # for project in ['TGL01MOH', 'HCCCFD', 'MOCHA', 'KLCS', 'ARCH1']:
# #     print(project)
# #     outputfile = os.path.join(outputfolder, project + '.provreport.json')
# #     write_json_to_file(outputfile, project, fpr, nabu_api='http://gsi-dcc.oicr.on.ca:3000', sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance')
        




# def extract_projects(fpr):
#     '''
#     (str) -> list

#     Returns a list of projects present in File Provenance Report

#     Parameters
#     ----------
#     - fpr (str): Path to the File Provenance Report    
#     '''
    
#     infile = open_fpr(fpr)
#     projects = set()
#     for line in infile:
#         projects.add(line.rstrip().split('\t')[1])
#     infile.close()
#     projects = sorted(list(projects))
#     return projects



# def validate_projects(pinery_projects, fpr_projects):
#     '''
#     (list, list) -> list
    
#     Returns a list of projects present both in File Provenance Report and in Pinery
    
#     Parameters
#     ----------
#     - pinery_projects (list): List of projects defined in Pinery
#     - fpr_projects (list): List of projects recorded in File Provenance Report
#     '''

#     valid_projects = set(pinery_projects).intersection(set(fpr_projects))
#     valid_projects = sorted(list(valid_projects))
#     return valid_projects

    





# ##### functions for new json ############

# def extract_project_info(project_provenance = 'http://pinery.gsi.oicr.on.ca/sample/projects'):
#     '''
#     (str) -> list
    
#     Returns a list of dictionary with project information pulled down from the
#     project_provenance Pinary API

#     Parameters
#     ----------
#     - project_provenance (str): Pinery API, http://pinery.gsi.oicr.on.ca/sample/projects
#     '''
    
#     response = requests.get(project_provenance)
#     if response.ok:
#         L = response.json()
#     else:
#         L = []
    
#     D = {}
    
#     if L:
#         for i in L:
#             name = i['name']
#             assert name not in D
#             D[name] = {'name': name}
#             for j in ['pipeline', 'description', 'active', 'contact_name', 'contact_email']:
#                 if j in i:
#                     D[name][j] = i[j] 
#                 else:
#                     D[name][j] = ''
#                 if j == 'active' and j in i:
#                     if i[j]:
#                         D[name][j] = 'Active'
#                     else:
#                         D[name][j] = 'Completed'
    
#     S = [D[i] for i in sorted(list(D.keys()))]
#     return S                            
    
        


# def collect_sequence_info(project, data):
#     '''    
#     (str, dict) -> dict

#     Returns a dictionary with sequencing information by extracting relevant information from data

#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - data (dict): Dictionary with worflow, file and qc info from parsing FPR, Nabu and Pinery
#     '''
    
#     # create a dict to store sequence info for each case
#     # {case: {library: {run: lane: [[file_path, file_id, read_count]]}}}
#     D = {}
    
#     for workflow in data[project]['workflows']:
#         if workflow.lower() in ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport']:
#            for case in data[project]['workflows'][workflow]:
#                for workflow_run in data[project]['workflows'][workflow][case]:
#                    info = data[project]['workflows'][workflow][case][workflow_run]['info']
#                    library = data[project]['workflows'][workflow][case][workflow_run]['libraries']
#                    files = data[project]['workflows'][workflow][case][workflow_run]['files']
#                    workflow_name = info['wf'] + '_' + info['wfv']
#                    assert len(data[project]['workflows'][workflow][case][workflow_run]['libraries']) == 1
#                    assert case == library[0]['id']
#                    lib = library[0]['lib']
#                    run = library[0]['run']
#                    lane = library[0]['lane']
#                    if 'read_count' in info:
#                        read_count = info['read_count']
#                    else:
#                        read_count = ''
                   
#                    # assert len(files) == 2
                                      
#                    if case not in D:
#                        D[case] = {}
#                    if lib not in D[case]:
#                        D[case][lib] = {}
#                    if run not in D[case][lib]:
#                        D[case][lib][run] = {}
#                    if lane not in D[case][lib][run]:
#                        D[case][lib][run][lane] = []
#                    for j in files:
#                        D[case][lib][run][lane].append([j['path'], j['fid'].split('.')[1], read_count, workflow_name])
#                    D[case][lib][run][lane].sort()
                   
#                    if len(D[case][lib][run][lane]) != 2:
#                        print(D[case][lib][run][lane])
                   
#                    #assert len(D[case][lib][run][lane]) == 2
    
    
#     return D



# def collect_lims_info(project, data):
#     '''
#     (str, dict) -> dict

#     Returns a dictionary with lims information by extracting relevant information from data

#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - data (dict): Dictionary with worflow, file and qc info from parsing FPR, Nabu and Pinery
#     '''
    
#     # create dict {case: {library: group_id}}
#     D = {}
    
#     for i in data[project]['lims']:
#         case = i['id']
#         for library in i['libraries']:
#             info = i['libraries'][library]['library_info']
#             group_id = info['group_id']
#             ext_id = info['ext_id']
#             if case not in D:
#                 D[case] = {}
#             D[case]['ext_id'] = ext_id
#             assert library not in D[case]
#             D[case][library] = group_id
#     return D    



# def collect_release_status(project, data):
#     '''
#     (str, dict) -> dict

#     Returns a dictionary with release status by extracting relevant information from data

#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - data (dict): Dictionary with worflow, file and qc info from parsing FPR, Nabu and Pinery
#     '''
    
#     D = {}
    
#     for i in data[project]['files']:
#         path = i['filepath']
#         file_id = i['fid']
#         status = i['qc']['status']
#         assert path not in D
#         D[path] = [file_id.split('.')[-1], status]
   
#     return D


# def add_lims_info_to_sequence_data(project, data):
#     '''
#     (str, dict) -> dict

#     Returns a dictionary with sequencing information for all cases in project by extracting relevant information from data

#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - data (dict): Dictionary with worflow, file and qc info from parsing FPR, Nabu and Pinery
#     '''
    
#     S = collect_sequence_info(project, data)
#     L = collect_lims_info(project, data)
#     R = collect_release_status(project, data)
    
#     {}
#     D = {}
    
    
#     # {case: {library: {run: lane: [[file_path, file_id, read_count, workflow]]}}}
#     # create dict {case: {sample: {library: group_id, 'ext_id': exit_id}}}
    
#     for case in S:
#         sample = case + '_' + L[case]['ext_id']   
#         for library in S[case]:
#             lib = library + '_' + L[case][library]
#             for run_id in S[case][library]:
#                  for lane in S[case][library][run_id]:
#                      run = run_id + '_' + lane
#                      files = []
#                      file_id = []
#                      workflow = ''
#                      read_count = ''
#                      release_status = ''
#                      for file in S[case][library][run_id][lane]:
#                          files.append(file[0])
#                          file_id.append(file[1])
#                          read_count = file[2]
#                          workflow = file[3]
#                          release_status = R[file[0]][-1]
#                          # assert file[1] == R[file[0]][0]
                     
#                      #assert len(files) == 2
#                      #assert len(file_id) == 2
                     
                     
#                      if case not in D:
#                          D[case] = {}
#                      if sample not in D[case]:
#                          D[case][sample] = {}
#                      if lib not in D[case][sample]:
#                          D[case][sample][lib] = {}
                     
#                      assert run not in D[case][sample][lib]
#                      D[case][sample][lib][run] = {'files': sorted(files), 'file_id': file_id, 'release': release_status, 'read_count': read_count, 'workflow': workflow}
                     
#     return D            





# app = Flask(__name__)

# @app.route('/')
# def index():
#     projects = extract_project_info()

#     # filter projects
#     valid = ['ARCH1', 'ARCHCVD', 'DLBCLR', 'HCCCFD', 'KLCS', 'MOCHA', 'SIMONE', 'TGL01MOH', 'TGL56']
#     projects = [i for i in projects if i['name'] in valid]

#     return render_template('index.html', projects=projects)

# @app.route('/<project_name>')
# def project(project_name):
#     projects = extract_project_info()
#     project = [i for i in projects if i['name'] == project_name][0]
#     return render_template('project.html', project=project)

# @app.route('/<project_name>/sequencing')
# def sequencing(project_name):
#     projects = extract_project_info()
#     project = [i for i in projects if i['name'] == project_name][0]

#     prov_data = {'ARCH1': 'ARCH1.prov.report.json',
#                  'ARCHCVD': 'ARCHCVD.prov.report.json',
#                  'DLBCLR': 'DLBCLR.prov.report.json',
#                  'HCCCFD': 'HCCCFD.prov.report.json',
#                  'KLCS': 'KLCS.prov.report.json',
#                  'MOCHA': 'MOCHA.prov.report.json',
#                  'SIMONE': 'SIMONE.prov.report.json',
#                  'TGL01MOH': 'TGL01MOH.prov.report.json',
#                  'TGL56': 'TGL56.prov.report.json'}

#     infile = open(prov_data[project_name])
#     data = json.load(infile)
#     infile.close()

#     sequences = add_lims_info_to_sequence_data(project_name, data)
        
#     return render_template('sequencing.html', project=project, sequences=sequences)



######################### db dvelopment


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
                 'geo_tissue_preparation', 'geo_receive_date', 'geo_group_id', 'geo_group_id_description'], ['tissue_type', 'ext_id', 'tissue_origin', 'library_type', 
                         'prep', 'tissue_prep', 'sample_received_date', 'group_id', 'group_id_description'])
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
    
    D = {project: {}}
    
    if L:
        for i in L:
            filepath = i['filepath']
            d = {}
            qc = collect_info(i, ['skip', 'user', 'date', 'qcstatus', 'ref', 'stalestatus'], ['skip', 'user', 'date', 'status', 'ref', 'fresh']) 
            d['qc'] = qc
            d['fid'] = 'f.' + str(i['fileswid'])
            d['filepath'] = filepath
            #assert filepath not in D[project]
            D[project][filepath] = d 
    
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
    
    D = {project: {}}
    
    
    for i in records:
        # get file path
        file = i[46]
        # get md5sums
        md5sum = i[47]
        # get file_swid
        file_swid = i[44]
        
        d = collect_info({k.split('=')[0]:k.split('=')[1] for k in i[17].split(';')},
                         ['geo_library_source_template_type'], ['library_type'])
        # get workflow
        workflow = i[30]
        # get workflow version
        version = i[31]        
        # get workflow run accession
        workflow_run = i[36]
        # collect file info
        D[project][file] = {'md5sum': md5sum, 'fid': file_swid, 'wf': workflow, 'wfv': version,
                   'wfrunid': workflow_run, 'file': file, 'library_type': d['library_type']}

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

    for project in file_info:
        for file in file_info[project]:
            # udpate qc info with file info
            if file not in qc_info[project]:
                qc_info[project][file] = {'qc': {'skip': '', 'user': '', 'date': '',
                                          'status': '', 'ref': '', 'fresh': ''},
                                          'fid': '', 'filepath': file}
            for i in file_info[project][file]:
                qc_info[project][file][i] = file_info[project][file][i]
            

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


def find_parent_workflows(project, fpr):
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
        assert project == i[1]
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
            parent_workflows = sorted(list(set([F[project][i] for i in P[project][workflow] if i in F[project]])))
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
    W, P, F = find_parent_workflows(project, fpr)
    # identify parents, children workflow relationships
    parents, children = identify_parent_children_workflows(W, P, F)

    D = {}

    records = get_FPR_records(project, fpr)
    
    for i in records:
        if project == i[1]:
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
            else:
                file_attributes = {}
            
            # get workflow info
            info = i[37]
            if info:
                info = info.split(';')
                info = {k.split('=')[0]: k.split('=')[1] for k in info if k.split('=')[0] not in ['cromwell-workflow-id', 'major_olive_version']}
            else:
                info = {}                
        
            # get library type
            d = collect_info({k.split('=')[0]:k.split('=')[1] for k in i[17].split(';')},
                             ['geo_library_source_template_type'], ['library_type'])
                    
            if 'libraries' not in D[project][workflow][sample][workflow_run]:
                D[project][workflow][sample][workflow_run]['libraries'] = []
            lib_info = {'lane': lane, 'run': run, 'limskey': limskey, 'id': sample, 'lib': library, 'library_type': d['library_type']}
            if lib_info not in D[project][workflow][sample][workflow_run]['libraries']:
                D[project][workflow][sample][workflow_run]['libraries'].append(lib_info)
                        
            if 'files' not in D[project][workflow][sample][workflow_run]:
                D[project][workflow][sample][workflow_run]['files'] = []
            file_info = {'md5sum': md5sum, 'wfrunid': 'wf.' + workflow_run, 'path': file,
             'fid': 'f.' + file_swid, 'wfv': workflow_version, 'wf': workflow}
            if file_info not in D[project][workflow][sample][workflow_run]['files']:
                D[project][workflow][sample][workflow_run]['files'].append(file_info)
                
            D[project][workflow][sample][workflow_run]['wfrunid'] = 'wf.' + workflow_run    

            D[project][workflow][sample][workflow_run]['info'] = {}
            for k in file_attributes:
                D[project][workflow][sample][workflow_run]['info'][k] = file_attributes[k]
            for k in info:
                D[project][workflow][sample][workflow_run]['info'][k] = info[k]
            D[project][workflow][sample][workflow_run]['info']['wfv'] = workflow_version
            D[project][workflow][sample][workflow_run]['info']['wf'] = workflow
                        
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


def extract_project_data(project, fpr, nabu_api='http://gsi-dcc.oicr.on.ca:3000', sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
    '''
    (str, str, str, str) -> dict

    Returns a dictionary with workflow, file and sample information for a given project    
        
    Parameters
    ----------
    - project (str): Name of project of interest
    - fpr (str): Path to the File Provenance Report
    - nabu_api (str): URL address of the Nabu API
    - sample_provenance (str): URL address of the Pinery sample provenance API
    '''

    # extract information about workflows
    workflow_info = extract_workflow_info(project, fpr)
    # extract qc and file information
    file_info = extract_file_info(project, fpr, nabu_api)
    # extract sample info
    lims = extract_lims(sample_provenance)

    data ={project: {'workflows': {}, 'files': [], 'lims': []}}

    for workflow in workflow_info[project]:
        #data[project]['workflows'].append({workflow: workflow_info[project][workflow]})
        data[project]['workflows'][workflow] = workflow_info[project][workflow]
        
    for file in file_info[project]:
        data[project]['files'].append(file_info[project][file])
    
    for sample in lims[project]:
        data[project]['lims'].append({'id': sample, 'libraries': lims[project][sample]})
    
    return data




def write_json_to_file(outputfile, project, fpr, nabu_api='http://gsi-dcc.oicr.on.ca:3000', sample_provenance='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance'):
    
    
    data = extract_project_data(project, fpr, nabu_api, sample_provenance)
    
    newfile = open(outputfile, 'w')
    json.dump(data, newfile,indent=4) 
    newfile.close()
    


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

    




##### functions for new json ############


def connect_to_db():
    '''
    (None) -> sqlite3.Connection
    
    Returns a connection to SqLite database prov_report.db.
    This database contains information extracted from FPR
    '''
    
    #conn = sqlite3.connect('prov_report_test.db')
    conn = sqlite3.connect('prov.test.db')
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
        read_count1 = json.loads(L[i]['attributes'])['read_count'] if 'read_count' in json.loads(L[i]['attributes']) else 'NA' 
        read_count2 = json.loads(L[i+1]['attributes'])['read_count'] if 'read_count' in json.loads(L[i+1]['attributes']) else 'NA' 

        if case1 == case2 and run1 == run2 and platform1 == platform2 \
        and library1 == library2 and sample1 == sample2 and wfrun1 == wfrun2:
            assert read_count1 == read_count2 
            assert workflow1 == workflow2
            assert json.loads(L[i]['attributes'])['read_number'] == '1' and json.loads(L[i+1]['attributes'])['read_number'] == '2'            
            
            d = {'case': case1, 'sample': sample1, 'library': library1, 'run': run1,
                 'files': [file1, file2], 'read_count': read_count1, 'workflow': workflow1,
                 'release_status': 'NA'}
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
        sample = '_'.join([case, i['tissue_origin'], i['tissue_type'], i['group_id']])
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
        D[case]['files'] = sorted(list(map(lambda x: os.path.basename(x), list(set(D[case]['files'])))))
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
    (str, str) -> list
    
    Returns a list of parent workflow ids (i.e immediate upstream workflow) upstream of the 
    workflow defined by workflow_id for a given project_name
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    conn = connect_to_db()
    
    data = conn.execute("SELECT Parents.parents_id FROM Parents WHERE Parents.project_id = '{0}' \
                        AND Parents.children_id = '{1}'".format(project_name, workflow_id)).fetchall()
    
    L = list(set([i['parents_id'] for i in data]))
    
    conn.close()
    
    return L


def get_children_workflows(project_name, workflow_id):
    '''
    (str, str) -> list
    
    Returns a list of children workflows (i.e workflow immediately downstream) 
    of workflow identified by workflow_id. 
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - bmpp_id (str): bamMergePreprocessing workflow id 
    '''
    
    conn = connect_to_db()
    
    data = conn.execute("SELECT Parents.children_id FROM Parents WHERE Parents.project_id = '{0}' \
                        AND Parents.parents_id = '{1}';".format(project_name, workflow_id)).fetchall()
    
    L = list(set([i['children_id'] for i in data]))
        
    conn.close()
    
    return L


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
    (str, list) -> list
    
    Returns a list of workflow ids removed from QC workflows
    
    Parameters
    ----------
    - project_name (str): name of project of interest
    - workflows (list): List of workflow ids
    '''

    conn = connect_to_db()

    L = []
    for workflow_id in workflows:
        data = conn.execute("SELECT Workflows.wfrun_id FROM Workflows WHERE Workflows.project_id = '{0}'\
                            AND Workflows.wfrun_id = '{1}' AND LOWER(Workflows.wf) NOT IN ('wgsmetrics_call_ready', 'insertsizemetrics_call_ready', \
                            'bamqc_call_ready');".format(project_name, workflow_id)).fetchall()
        if data:
            assert len(data) == 1
            L.append(data[0]['wfrun_id'])
        
    conn.close()
    
    return L



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
    bwamem_ids = get_parent_workflows(project_name, bmpp_id)
    # get the fastq-generating worflow ids
    fastqs_workflow_ids = []
    for workflow_id in bwamem_ids:
        fastqs_workflow_ids.extend(get_parent_workflows(project_name, workflow_id))
    
    fastqs_workflow_ids = list(set(fastqs_workflow_ids))
    
    conn = connect_to_db()
    
    # track release status of all fastqs 
    D = {}
    
    
    
    ##### COMMENTED CODE TO WORK AROUND NABU ####
    
    # get the file swids of the fastq-generating workflows
    for workflow_id in fastqs_workflow_ids:
        # data = conn.execute("SELECT Workflows.wf, Files.file_swid, FilesQC.status  \
        #                       FROM Workflows JOIN Files JOIN FilesQC WHERE Files.project_id = '{0}' \
        #                       AND Workflows.project_id = '{0}' AND FilesQC.project_id = '{0}' AND  \
        #                       Files.wfrun_id = '{1}' AND FilesQC.file_swid = Files.file_swid AND \
        #                       Workflows.wfrun_id = '{1}';".format(project_name, workflow_id)).fetchall()
            
        # # skip import workflows because fastqs from these workflow may not need to be shared back
        # for i in data:
        #     if 'import' not in i['wf'].lower():
        #         assert i['file_swid'] not in D
        #         D[i['file_swid']] = i['status']
    
        
        data = conn.execute("SELECT Workflows.wf, Files.file_swid FROM Workflows JOIN Files \
                            WHERE Files.project_id = '{0}' AND Workflows.project_id = '{0}' AND \
                            Files.wfrun_id = '{1}' AND Workflows.wfrun_id = '{1}';".format(project_name, workflow_id)).fetchall()
        for i in data:
            if 'import' not in i['wf'].lower():
                assert i['file_swid'] not in D
                D[i['file_swid']] = 'NA'
             
    if D:
        if all(map(lambda x: x.lower() == 'pass', D.values())):
            return True
        else:
            return False
    else:
        return False
       
        
    conn.close()
    










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
                    
        D[libraries] = {'sample': sample, 'workflow': data[0]['wf'].split('_')[0],
                        'workflow_id': data[0]['wfrun_id'], 'version': data[0]['wfv'], 'attributes': data[0]['attributes']}
        
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

    # make a list of expected workflows
    expected_workflows = {'mutect2': 'variantEffectPredictor',
                          'varscan': 'sequenza',
                          'delly': 'mavis'}
    
    D = {}
    
    for workflow_id in downstream_workflows:
        # group all downstream workflows per library pair
        d = get_workflow_info(project_name, workflow_id)
        assert d
        libraries = list(d.keys())[0]
        # get the downstream workflow (ie mavis, VEP, delly)
        child_workflow = get_children_workflows(project_name, workflow_id)
        child_workflow = filter_out_QC_workflows(project_name, child_workflow)
        if child_workflow:
            child_workflow = child_workflow[0]    
            w = get_workflow_info(project_name, child_workflow)
        else:
            w = {}
            w[libraries] = {'sample': d[libraries]['sample'], 
                            'workflow': expected_workflows[d[libraries]['workflow']],
                            'workflow_id': 'NA', 'version': 'NA', 'attributes': 'NA'}       
        if libraries not in D:
            D[libraries] = []
        D[libraries].append(d[libraries])
        D[libraries].append(w[libraries])
    
    # sort workflows   
    expected_workflows = ['mutect2', 'variantEffectPredictor', 'varscan', 'sequenza', 'delly', 'mavis']
    for libraries in D:
        a = [0 for i in range(len(D[libraries]))]
        for i in D[libraries]:
            j = expected_workflows.index(i['workflow'])
            a[j] = i
        D[libraries] = a
    
    # get the files for each workflow
    for libraries in D:
        for i in range(len(D[libraries])):
            files, creation_date = get_workflow_files(project_name, D[libraries][i]['workflow_id'])
            D[libraries][i]['files'] = files
            # convert epoch time to standard time
            if creation_date:
                creation_date = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(creation_date)))
            D[libraries][i]['creation_date'] = creation_date
    
    return D


def format_bmpp_dowmstream_workflows(D):
    '''
    (dict) -> list
    
    Returns a list with information about workflows downstream a specific bmpp workflow run 
    to display in the html table
    
    Parameters
    ----------
    - D (dict): Dictionary with information about bmpp downstream workflows
    '''

    downstream_workflows = []
    for libraries in D:
        normal, tumour = D[libraries][0]['sample'].split(';')
        d = [normal, tumour]
        for i in D[libraries]:
            d.append(i['workflow_id'])
        downstream_workflows.append(d)      
    return downstream_workflows



def format_dowmstream_workflow_files(D, bmpp_id):
    '''
    (dict, str) -> list    
    
    Returns a list with information and files generated by workflows downstream pecific bmpp workflow run bmpp_id
    to display in the html table
    
    Parameters
    ----------
    - D (dict): Dictionary with information about bmpp downstream workflows
    - bmpp_id (str): Specific bmpp workflow run id
    '''

    # get the workflow files
    downstream_workflow_files = []
    for libraries in D:
        for i in D[libraries]:
            if i['workflow_id'] != 'NA':
                if i['attributes']:
                    attributes = i['attributes'].replace("\\\"", "").replace('\\', '')
                    attributes = json.loads(attributes)
                    if 'reference' in attributes:
                        attributes = attributes['reference'].replace('"', '')
                else:
                    attributes = 'NA'
                L = [i['workflow_id'], i['workflow'], i['creation_date'],
                     attributes, ['bamMergePreprocessing', bmpp_id], [os.path.dirname(i['files'][0])] + sorted(map(lambda x: os.path.basename(x), i['files']))]
                downstream_workflow_files.append(L)
    return downstream_workflow_files






    

# map pipelines to views
routes = {'Whole Genome': 'whole_genome_sequencing'}

    


app = Flask(__name__)

@app.route('/')
def index():
    
    # connect to db and extract project info
    conn = connect_to_db()
    projects = conn.execute('SELECT * FROM Projects').fetchall()
    conn.close()
    
    return render_template('index.html', projects=projects)

@app.route('/<project_name>')
def project(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
    
    return render_template('project.html', routes=routes, project=project, pipelines=pipelines)

@app.route('/<project_name>/sequencing')
def sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
    
    
    conn = connect_to_db()

    # files = conn.execute("SELECT Files.file, Files.workflow, Files.version, Files.wfrun_id, \
    #                      Files.attributes, FilesQC.status from Files JOIN FilesQC WHERE Files.project_id = '{0}' \
    #                      AND FilesQC.project_id = '{0}' AND Files.file_swid = FilesQC.file_swid ;".format(project_name)).fetchall()
     
    # extract file information for fastq-generating workflows
    files = conn.execute("SELECT Files.file, Files.workflow, Files.version, Files.wfrun_id, Files.attributes, \
                         Workflow_Inputs.run, Workflow_Inputs.lane, Workflow_Inputs.platform, \
                         Libraries.library, Libraries.sample, Libraries.ext_id, Libraries.group_id \
                         from Files JOIN Workflow_Inputs JOIN Libraries WHERE Files.project_id = '{0}' \
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

    # extract sample, library and workflow information for call ready workflow bamMergePreprocessing
    data = conn.execute("SELECT Files.creation_date, Libraries.library, Libraries.sample, \
                         Libraries.ext_id, Libraries.group_id, Libraries.tissue_type, \
                         Libraries.tissue_origin, Workflow_Inputs.run, Workflow_Inputs.lane, \
                         Workflow_Inputs.platform, Workflows.wf, Workflows.wfv, Workflows.wfrun_id \
                         from Files JOIN Libraries JOIN Workflow_Inputs JOIN Workflows \
                         WHERE Files.project_id = '{0}' AND Libraries.project_id = '{0}' \
                         AND Workflow_Inputs.project_id = '{0}' AND Workflows.project_id = '{0}' \
                         AND Files.wfrun_id = Workflow_Inputs.wfrun_id  AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id AND Workflow_Inputs.library = Libraries.library \
                         AND LOWER(Workflows.wf) = 'bammergepreprocessing'".format(project_name)).fetchall()
    
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
                         AND LOWER(Workflows.wf) = 'bammergepreprocessing'AND Libraries.sample ='{1}'".format(project_name, case)).fetchall()
    conn.close()

    # get sample, library and file info for for the most recent bmpp run for case in project
    bmpp_id, bmpp_files, bmpp_info = get_bmpp_files(data)
    
    # get QC status of bmpp input fastq files
    fastq_status = bmpp_input_raw_seq_status(project_name, bmpp_id)

    # get the bmpp downstream workflows
    D = get_bmpp_downstream_workflows(project_name, bmpp_id)
    downstream_workflows = format_bmpp_dowmstream_workflows(D)
    # get the workflow files
    downstream_workflow_files = format_dowmstream_workflow_files(D, bmpp_id)
    
       
    return render_template('WGS_case.html', routes = routes, fastq_status=fastq_status,
                           bmpp_info=bmpp_info, bmpp_id=bmpp_id, bmpp_files=bmpp_files,
                           sample_case=case, project=project, pipelines=pipelines,
                           downstream_workflows=downstream_workflows,
                           downstream_workflow_files=downstream_workflow_files)


