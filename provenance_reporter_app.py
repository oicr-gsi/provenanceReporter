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
import itertools

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import io
import base64

from utilities import connect_to_db, get_children_workflows, get_miso_sample_link,\
    get_pipelines, get_workflow_names
from whole_genome import get_call_ready_cases, get_bmpp_case, get_case_call_ready_samples, group_normal_tumor_pairs, \
    find_analysis_blocks, map_workflows_to_sample_pairs, map_workflows_to_parent, list_block_workflows, \
    get_block_analysis_date, sort_call_ready_samples, get_block_workflow_file_count, get_block_release_status, \
    get_amount_data, is_block_complete, order_blocks, name_WGS_blocks, create_block_json, map_samples_to_bmpp_runs, \
    get_parent_workflows, get_workflows_analysis_date, get_workflow_file_count, \
    get_workflow_limskeys, get_file_release_status    
from networks import get_node_labels, make_adjacency_matrix, plot_workflow_network
from whole_transcriptome import get_WT_call_ready_cases, get_star_case, get_WT_case_call_ready_samples, \
    map_workflows_to_samples, find_WT_analysis_blocks, name_WT_blocks
from project import get_project_info, get_cases, get_sample_counts, add_missing_donors, get_last_sequencing
from sequencing import get_sequences, collect_sequence_info


   

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
        
           

def get_samples(project_name):
    '''
    
    
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT library, sample, ext_id, group_id, group_id_description, library_type, tissue_origin, tissue_type FROM Libraries WHERE project_id = '{0}'".format(project_name)).fetchall()
    
    data = list(set(data))
    
    return data
    




    
    






# map pipelines to views
routes = {'Whole Genome': 'whole_genome_sequencing', 'Whole Transcriptome': 'whole_transcriptome'}

    


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
    


@app.template_filter()
def shorten_workflow_id(workflow_run_id):
    '''
    (str) -> str
    
    Shorten the workflow run id to 8 characters + trailing dots
             
    Parameters
    ----------
    - workflow_run_id (str): Workflow unique run identifier
    '''
    
    return workflow_run_id[:8] + '...'
    
    
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
    # get library and sample counts
    counts = get_sample_counts(project_name)
    # add missing donors to counts (ie, when counts are 0)
    counts = add_missing_donors(cases, counts)
    # get the date of the last sequencing data
    seq_date = get_last_sequencing(project['project_id'])
    
    return render_template('project.html', routes=routes, project=project, pipelines=pipelines, cases=cases, counts=counts, seq_date=seq_date)


@app.route('/<project_name>/sequencing')
def sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
    # get sequence file information
    files = collect_sequence_info(project_name)
    # re-organize sequence information
    sequences = get_sequences(files)

    return render_template('sequencing.html', routes=routes, project=project, sequences=sequences, pipelines=pipelines)



@app.route('/<project_name>/whole_genome_sequencing')
def whole_genome_sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
        
    # get samples and libraries and workflow ids for each case
    cases = get_call_ready_cases(project_name, 'novaseq', 'WG')
    samples = sorted(list(cases.keys()))

   
    return render_template('Whole_Genome_Sequencing.html', routes = routes, project=project, samples=samples, cases=cases, pipelines=pipelines)


@app.route('/<project_name>/whole_genome_sequencing/<case>')
def wgs_case(project_name, case):
    
    tasks = []
    
    start = time.time()
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    end1 = time.time()
    print('project', end1 - start)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
    end2 = time.time()
    print('pipeline', end2 - end1)
    
    # get all the bmpp runs for WG library type and Novaseq platform
    bmpp = get_bmpp_case(project_name, case, 'novaseq', 'WG')
    end3 = time.time()
    print('bmpp', end3 - end2)
       
    
    # get the normal and tumor samples for each bmpp id
    bmpp_samples = map_samples_to_bmpp_runs(project_name, bmpp)
    end4 = time.time()
    print('bmpp_samples', end4 - end3)
    
    # identify all the samples processed
    samples = get_case_call_ready_samples(project_name, bmpp_samples)
    end5 = time.time()
    print('samples', end5 - end4)
    
    # get all pairs N/T samples
    pairs = group_normal_tumor_pairs(samples)
    end6 = time.time()
    print('pairs', end6 - end5)
    
    
    # find analysis workflows for each N/T pairs
    # remove sample pairs without analysis workflows
    D = map_workflows_to_sample_pairs(project_name, 'novaseq', pairs)
    end7 = time.time()
    print('analysis workflows', end7 - end6)
    
       
    # find the parents of each workflow
    parents = get_parent_workflows(project_name)
    end8 = time.time()
    print('parents', end8 - end7)
    
    # get the parent workflows for each block
    parent_workflows = map_workflows_to_parent(D, parents)
    end9 = time.time()
    print('parent workflows', end9 - end8)
    
    # find the blocks by mapping the analysis workflows to their parent workflows    
    blocks = find_analysis_blocks(D, parents, parent_workflows, bmpp)
    end10 = time.time()
    print('blocks', end10 - end9)
    
    # list all workflows for each block
    block_workflows = list_block_workflows(blocks)
    end11 = time.time()
    print('block workflows', end11 - end10)
    
    # get the workflow creation date for all the workflows in project
    creation_dates = get_workflows_analysis_date(project_name)
    # assign date to each block. most recent file creation date from all workflows within block 
    # get the date of each workflow within block
    block_date, workflow_date = get_block_analysis_date(block_workflows, creation_dates)
    end12 = time.time()
    print('workflow date', end12 - end11)
    
    # map each workflow run id to its workflow name
    workflow_names = get_workflow_names(project_name)
    # get the workflow names
    block_workflow_names = get_node_labels(block_workflows, workflow_names)
    end13 = time.time()
    print('workflow names', end13 - end12)
    
    # convert workflow relationships to adjacency matrix for each block
    matrix = make_adjacency_matrix(block_workflows, parent_workflows)
    end14 = time.time()
    print('matrix', end14 - end13)
    
    # create figures
    figures = plot_workflow_network(matrix, block_workflow_names)
    end15 = time.time()
    print('figures', end15 - end14)
    
     
    # get the samples for each bmpp id
    samples_bmpp = sort_call_ready_samples(project_name, blocks, bmpp_samples, workflow_names)
    end16 = time.time()
    print('samples bmpp', end16 - end15)
    
    
    # get the file count of each workflow in project
    file_counts = get_workflow_file_count(project_name)
    # get the workflow file counts
    file_counts = get_block_workflow_file_count(block_workflows, file_counts)
    end17 = time.time()
    print('file counts', end17 - end16)
    
    # get release status of input sequences for each block
    # get the input limskeys for each workflow in project
    limskeys = get_workflow_limskeys(project_name)
    end18 = time.time()
    print('limskeys', end18 - end17)
    
      
    # get the file swid and release status for each limskey for fastq-generating workflows
    # excluding fastq-import workflows
    status = get_file_release_status(project_name)
    release_status = get_block_release_status(block_workflows, limskeys, status)
    end19 = time.time()
    print('release status', end19 - end18)
    
    # get the amount of data for each workflow
    amount_data = get_amount_data(block_workflows, limskeys)
    end20 = time.time()
    print('amount data', end20 - end19)
    
    
    
    # check if blocks are complete
    expected_workflows = sorted(['mutect2', 'variantEffectPredictor', 'delly', 'varscan', 'sequenza', 'mavis'])           
    complete = is_block_complete(blocks, expected_workflows)
    end21 = time.time()
    print('complete', end21 - end20)
    
    
    
    # order blocks based on the amount of data
    ordered_blocks = order_blocks(blocks, amount_data)
    end22 = time.time()
    print('order blocks', end22 - end21)


    # name each block according to the selected block order
    names = name_WGS_blocks(ordered_blocks)
    end23 = time.time()
    print('name blocks', end23 - end22)    
    
    # get miso link
    miso_link = get_miso_sample_link(project_name, case)
    end24 = time.time()
    print('miso link', end24 - end23)
    
    
    # sort sample pairs names
    sample_pairs_names = sorted(list(blocks.keys()))
    end25 = time.time()
    print('sort sample pairs', end25 - end24)
    
    
    print('all tasks', end25 -  start)
    
    tasks.extend([('project', end1 - start), ('pipeline', end2 - end1),
                 ('bmpp', end3 - end2), ('bmpp_samples', end4 - end3),
                 ('samples', end5 - end4), ('pairs', end6 - end5),
                 ('analysis workflows', end7 - end6), ('parents', end8 - end7),
                 ('parent workflows', end9 - end8), ('blocks', end10 - end9),
                 ('block workflows', end11 - end10), ('workflow date', end12 - end11),
                 ('workflow names', end13 - end12), ('matrix', end14 - end13),
                 ('figures', end15 - end14), ('samples bmpp', end16 - end15),
                 ('file counts', end17 - end16), ('limskeys', end18 - end17),
                 ('release status', end19 - end18), ('amount data', end20 - end19),
                 ('complete', end21 - end20), ('order blocks', end22 - end21),
                 ('name blocks', end23 - end22), ('miso link', end24 - end23),
                 ('sort sample pairs', end25 - end24)])
    tasks.sort(key=lambda x: x[1])
    print(tasks)
    
    
    
    
    
    
    
    return render_template('WGS_case.html', project=project, routes = routes,
                           case=case, pipelines=pipelines, sample_pairs_names=sample_pairs_names,
                           blocks=blocks, names=names, ordered_blocks=ordered_blocks,
                           miso_link=miso_link, complete=complete, release_status=release_status,
                           block_date=block_date, workflow_date=workflow_date,
                           figures=figures, samples_bmpp=samples_bmpp, 
                           file_counts=file_counts, amount_data=amount_data, )



@app.route('/<project_name>/whole_transcriptome')
def whole_transcriptome(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
        
    # get samples and libraries and workflow ids for each case
    cases = get_WT_call_ready_cases(project_name, 'novaseq', 'WT')
    samples = sorted(list(cases.keys()))

    return render_template('Whole_transcriptome.html', routes = routes, project=project,
                           samples=samples, cases=cases, pipelines=pipelines)



@app.route('/<project_name>/whole_transcriptome/<case>')
def wt_case(project_name, case):
    
    # get the project info for project_name from db
    project = get_project_info(project_name)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name)
        
    # build the somatic calling block

    # identify all call ready star runs for novaseq
    star = get_star_case(project_name, case, 'novaseq', 'WT')
    
    # identify the samples processed
    samples = get_WT_case_call_ready_samples(project_name, star)
    
    # remove samples without analysis workflows
    D = map_workflows_to_samples(project_name, 'novaseq', samples)

    # get the parent workflows for each block
    parent_workflows = map_workflows_to_parent(project_name, D)
        
    blocks = find_WT_analysis_blocks(project_name, D, parent_workflows, star)
    
    # list all workflows for each block
    block_workflows = list_block_workflows(blocks)
    
    # assign date to each block. most recent file creation date from all workflows within block 
    # get the date of each workflow within block
    block_date, workflow_date = get_block_analysis_date(block_workflows)
    
    # get the workflow names
    workflow_names = get_node_labels(block_workflows)

    # convert workflow relationships to adjacency matrix for each block
    matrix = make_adjacency_matrix(block_workflows, parent_workflows)
                                
    # create figures
    figures = plot_workflow_network(matrix, workflow_names)
    
    # get the samples for each star id
    samples_star = sort_call_ready_samples(project_name, blocks)
    
    # get the workflow file counts
    file_counts = get_block_workflow_file_count(block_workflows)
    
    # get release status of input sequences for each block
    release_status = get_block_release_status(block_workflows)
    
    # get the amount of data for each workflow
    amount_data = get_amount_data(block_workflows)
    
    # check if blocks are complete
    expected_workflows = sorted(['arriba', 'rsem', 'star', 'starfusion', 'mavis'])           
    complete = is_block_complete(blocks, expected_workflows)
    
    # order blocks based on the amount of data
    ordered_blocks = order_blocks(blocks, amount_data)

    # name each block according to the selected block order
    names = name_WT_blocks(ordered_blocks)
    
    # get miso link
    miso_link = get_miso_sample_link(project_name, case)
    
    # sort sample pairs names
    sample_pairs_names = sorted(list(blocks.keys()))
    
    
    return render_template('WT_case.html', project=project, routes = routes,
                           case=case, pipelines=pipelines, sample_pairs_names=sample_pairs_names,
                           blocks=blocks, names=names, ordered_blocks=ordered_blocks,
                           miso_link=miso_link, complete=complete, release_status=release_status,
                           block_date=block_date, workflow_date=workflow_date,
                           figures=figures, samples_star=samples_star, 
                           file_counts=file_counts, amount_data=amount_data, )


@app.route('/download_wgs_block/<project_name>/<case>/<block>/<bmpp_parent>')
def download_block_data(project_name, case, block, bmpp_parent):
    '''
    
    
    '''
    
    # get all the bmpp runs for WG library type and Novaseq platform
    bmpp = get_bmpp_case(project_name, case, 'novaseq', 'WG')
    # get the normal and tumor samples for each bmpp id
    bmpp_samples = map_samples_to_bmpp_runs(project_name, bmpp)
    # identify all the samples processed
    samples = get_case_call_ready_samples(project_name, bmpp_samples)
    # get all pairs N/T samples
    pairs = group_normal_tumor_pairs(samples)
    # find analysis workflows for each N/T pairs
    # remove sample pairs without analysis workflows
    D = map_workflows_to_sample_pairs(project_name, 'novaseq', pairs)
    # get the parents of each workflow
    parents = get_parent_workflows(project_name)
    # get the parent workflows for each block
    parent_workflows = map_workflows_to_parent(D, parents)
    # find the blocks by mapping the analysis workflows to their parent workflows    
    blocks = find_analysis_blocks(D, parents, parent_workflows, bmpp)
    # map each workflow run id to its workflow name
    workflow_names = get_workflow_names(project_name)
    # create json with workflow information for block for DARE
    block_data = create_block_json(project_name, blocks, block, bmpp_parent, workflow_names)
    
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(block_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}_WGS_{1}_{2}.json".format(project_name, case, block)})


@app.route('/download_wt_block/<project_name>/<case>/<block>/<star_parent>')
def download_WT_block_data(project_name, case, block, star_parent):

    

    # build the somatic calling block

    # identify all call ready star runs for novaseq
    star = get_star_case(project_name, case, 'novaseq', 'WT')
    
    # identify the samples processed
    samples = get_WT_case_call_ready_samples(project_name, star)
    
    # remove samples without analysis workflows
    D = map_workflows_to_samples(project_name, 'novaseq', samples)

    # get the parent workflows for each block
    parent_workflows = map_workflows_to_parent(project_name, D)
    
    # find the blocks
    blocks = find_WT_analysis_blocks(project_name, D, parent_workflows, star)
    
    # create json with workflow information for block for DARE
    block_data = create_block_json(project_name, blocks, block, star_parent)

    # send the json to outoutfile                    
    return Response(
        response=json.dumps(block_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}_WT_{1}_{2}.json".format(project_name, case, block)})

    

@app.route('/download_cases/<project_name>')
def download_cases_table(project_name):
    '''
    (str) -> None
    
    Download a table with project information in Excel format
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    # get case information
    cases = get_cases(project_name)
    
    # get library and sample counts
    counts = get_sample_counts(project_name)
    # add missing donors to counts (ie, when counts are 0)
    counts = add_missing_donors(cases, counts)
    
    D = {}
    for i in cases:
        donor = i['case_id']
        D[donor] = i
        D[donor]['library'] = counts[donor]['library']
        D[donor]['normal'] = counts[donor]['normal']
        D[donor]['tumor'] = counts[donor]['tumor']
        
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


# if __name__ == "__main__":
#     app.run(host='0.0.0.0')
    