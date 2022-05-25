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
            assert len(i['libraries']) == 1
            case = i['libraries'][0]['id']
            lib = i['libraries'][0]['lib']
            run = i['libraries'][0]['run']
            lane = i['libraries'][0]['lane']
            read_count = i['info']['read_count']
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
                D[case][lib][run][lane].append([j['path'], j['fid'].split('.')[1], read_count])
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
    
    
    # {case: {library: {run: lane: [[file_path, file_id, read_count]]}}}
    # create dict {case: {library: group_id, 'ext_id': exit_id}}
    
    for case in S:
        sample = case + '_' + L[case]['ext_id']   
        for library in S[case]:
            lib = library + '_' + L[case][library]
            for run_id in S[case][library]:
                 for lane in S[case][library][run_id]:
                     run = run_id + '_' + lane
                     files = []
                     file_id = []
                     for file in S[case][library][run_id][lane]:
                         files.append(file[0])
                         file_id.append(file[1])
                         read_count = file[2]
                         release_status = R[file[0]][-1]
                         assert file[1] == R[file[0]][0]
                     assert len(files) == 2
                     assert len(file_id) == 2
                     
                     if sample not in D:
                         D[sample] = {}
                     if lib not in D[sample]:
                         D[sample][lib] = {}
                     
                     assert run not in D[sample][lib]
                     D[sample][lib][run] = {'files': files, 'file_id': file_id, 'release': release_status, 'read_count': read_count}
                     
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