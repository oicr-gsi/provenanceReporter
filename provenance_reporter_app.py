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

from utilities import connect_to_db, get_miso_sample_link,\
    get_pipelines, get_workflow_names, get_library_design
from whole_genome import get_call_ready_cases, get_bmpp_case, get_case_call_ready_samples, group_normal_tumor_pairs, \
    find_analysis_blocks, map_workflows_to_sample_pairs, map_workflows_to_parent, list_block_workflows, \
    get_block_analysis_date, sort_call_ready_samples, get_block_release_status, \
    get_amount_data, is_block_complete, order_blocks, create_block_json, map_samples_to_bmpp_runs, \
    get_parent_workflows, get_workflows_analysis_date, get_workflow_file_count, \
    get_workflow_limskeys, get_file_release_status, get_WGS_blocks_info, get_sequencing_platform    
from networks import get_node_labels, make_adjacency_matrix, plot_workflow_network
from whole_transcriptome import get_WT_call_ready_cases, get_star_case, get_WT_case_call_ready_samples, \
    map_workflows_to_samples, find_WT_analysis_blocks, name_WT_blocks, map_samples_to_star_runs
from project import get_project_info, get_cases, get_sample_counts, count_libraries, \
     get_library_types, add_missing_donors, get_last_sequencing
from sequencing import get_sequences, collect_sequence_info


   
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
def readable_time(date):
    '''
    (str) -> str
    
    Returns epoch time in readable format
    
    Parameters
    ----------
    - date (str): Epoch time
    '''
    
    #return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(date)))
    return time.strftime('%Y-%m-%d', time.localtime(int(date)))


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


@app.template_filter()
def format_created_time(created_time):
    '''
    (str) -> str
    
    Remove time in created time and keep only the date
                 
    Parameters
    ----------
    - created_time (str): Time a sample is created in the format Year-Month-DayTHour:Mn:SecZ
    '''
    
    return created_time[:created_time.index('T')]
    
    
@app.route('/')
def index():
    
    # connect to db and extract project info
    conn = connect_to_db('merged.db')
    projects = conn.execute('SELECT * FROM Projects').fetchall()
    conn.close()
    
    projects = sorted([(i['project_id'], i) for i in projects])
    projects = [i[1] for i in projects]
    
    return render_template('index.html', projects=projects)

@app.route('/<project_name>')
def project_page(project_name):
    
    
    # get the project info for project_name from db
    project = get_project_info(project_name, 'merged.db')
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, 'merged.db')
    # get case information
    cases = get_cases(project_name, 'merged.db')
    # get the species
    species = ', '.join(sorted(list(set([i['species'] for i in cases]))))
    # get library and sample counts
    counts = get_sample_counts(project_name, 'merged.db')
    # add missing donors to counts (ie, when counts are 0)
    counts = add_missing_donors(cases, counts)
    # count libraries for each library type
    # get the library types
    library_types =  get_library_types(project_name, 'merged.db')
    # get the library names
    library_names = {i: get_library_design(i) for i in library_types}
    libraries = count_libraries(project_name, library_types, cases, 'merged.db')
    # get the date of the last sequencing data
    seq_date = get_last_sequencing(project['project_id'], 'merged.db')
    
    return render_template('project.html', routes=routes, project=project,
                           pipelines=pipelines, cases=cases, counts=counts,
                           seq_date=seq_date, species=species, libraries=libraries,
                           library_types = library_types, library_names=library_names)


@app.route('/<project_name>/sequencing')
def sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name, 'merged.db')
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, 'merged.db')
    # get sequence file information
    files = collect_sequence_info(project_name, 'merged.db')
    # re-organize sequence information
    sequences = get_sequences(files)

    return render_template('sequencing.html', routes=routes, project=project, sequences=sequences, pipelines=pipelines)



@app.route('/<project_name>/whole_genome_sequencing')
def whole_genome_sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name, 'merged.db')
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, 'merged.db')
        
    # get samples and libraries and workflow ids for each case
    cases = get_call_ready_cases(project_name, 'novaseq', 'WG', 'merged.db')
    samples = sorted(list(cases.keys()))

   
    return render_template('Whole_Genome_Sequencing.html', routes = routes, project=project, samples=samples, cases=cases, pipelines=pipelines)



@app.route('/<project_name>/whole_genome_sequencing/<case>')
def wgs_case(project_name, case):
    
    
    # get the project info for project_name from db
    project = get_project_info(project_name, 'merged.db')
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, 'merged.db')
    # get miso link
    miso_link = get_miso_sample_link(project_name, case, 'merged.db')
    # get the WGS blocks
    blocks = get_WGS_blocks_info(project_name, case, 'merged.db')
    # sort sample pairs names
    sample_pairs_names = sorted(list(blocks.keys()))
    # get the workflow names
    workflow_names = get_workflow_names(project_name, 'merged.db')
    # get the file count of each workflow in project
    file_counts = get_workflow_file_count(project_name, 'merged.db')
    # get the amount of data for each workflow
    amount_data = get_amount_data(project_name, 'merged.db')
    # get the creation date of all workflows
    creation_dates = get_workflows_analysis_date(project_name, 'merged.db')
    # get the sequencing platform of each workflow
    platforms = get_sequencing_platform(project_name, 'merged.db')
    # find the parents of each workflow
    parents = get_parent_workflows(project_name, 'merged.db')
    
    return render_template('WGS_case.html',
                           project=project,
                           routes = routes,
                           case=case,
                           pipelines=pipelines,
                           blocks=blocks,
                           sample_pairs_names=sample_pairs_names,
                           workflow_names=workflow_names,
                           file_counts=file_counts,
                           amount_data=amount_data,
                           creation_dates=creation_dates,
                           platforms=platforms,
                           parents=parents
                           )


@app.route('/<project_name>/whole_transcriptome')
def whole_transcriptome(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(project_name, 'merged.db')
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, 'merged.db')
    # get samples and libraries and workflow ids for each case
    cases = get_WT_call_ready_cases(project_name, 'novaseq', 'merged.db', 'WT')
    samples = sorted(list(cases.keys()))

    return render_template('Whole_transcriptome.html', routes = routes, project=project,
                           samples=samples, cases=cases, pipelines=pipelines)


@app.route('/<project_name>/whole_transcriptome/<case>')
def wt_case(project_name, case):
    
    
    tasks = []
    start = time.time()
       
    
    # get the project info for project_name from db
    project = get_project_info(project_name, 'merged.db')
    end1 = time.time()
    print('project', end1 - start)
    
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, 'merged.db')
    end2 = time.time()
    print('pipeline', end2 - end1)
           
    # build the somatic calling block

    # identify all call ready star runs for novaseq
    star = get_star_case(project_name, case, 'novaseq', 'WT', 'merged.db')
    end3 = time.time()
    print('star', end3 - end2)
    
    # get the tumor samples for each star id
    star_samples = map_samples_to_star_runs(project_name, star, 'merged.db')
    end4 = time.time()
    print('star_samples', end4 - end3)
    
    # identify all the samples processed
    samples = get_WT_case_call_ready_samples(project_name, star_samples)
    end5 = time.time()
    print('samples', end5 - end4)
    
    # remove samples without analysis workflows
    D = map_workflows_to_samples(project_name, 'novaseq', samples, 'merged.db')
    end6 = time.time()
    print('analysis workflows', end6 - end5)

    # find the parents of each workflow
    parents = get_parent_workflows(project_name, 'merged.db')
    end7 = time.time()
    print('parents', end7 - end6)

    # get the parent workflows for each block
    parent_workflows = map_workflows_to_parent(D, parents)
    end8 = time.time()
    print('parent workflows', end8 - end7)     

    
    # find the blocks by mapping the analysis workflows to their parent workflows    
    blocks = find_WT_analysis_blocks(D, parents, parent_workflows, star)
    end9 = time.time()
    print('blocks', end9 - end8)
    
        
    # list all workflows for each block
    block_workflows = list_block_workflows(blocks)
    end10 = time.time()
    print('block workflows', end10 - end9)
    
    # get the workflow creation date for all the workflows in project
    creation_dates = get_workflows_analysis_date(project_name)
    # assign date to each block. most recent file creation date from all workflows within block 
    # get the date of each workflow within block
    block_date, workflow_date = get_block_analysis_date(block_workflows, creation_dates)
    end11 = time.time()
    print('workflow date', end11 - end10)  
    
    # map each workflow run id to its workflow name
    workflow_names = get_workflow_names(project_name, 'merged.db')
    # get the workflow names
    block_workflow_names = get_node_labels(block_workflows, workflow_names)
    end12 = time.time()
    print('workflow names', end12 - end11)
    
    # convert workflow relationships to adjacency matrix for each block
    matrix = make_adjacency_matrix(block_workflows, parent_workflows)
    end13 = time.time()
    print('matrix', end13 - end12)
 
    # create figures
    figures = plot_workflow_network(matrix, block_workflow_names)
    end14 = time.time()
    print('figures', end14 - end13)
 
    # get the samples for each star id
    samples_star = sort_call_ready_samples(project_name, blocks, star_samples, workflow_names)
    end15 = time.time()
    print('samples bmpp', end15 - end14)
   
    # # get the file count of each workflow in project
    file_counts = get_workflow_file_count(project_name, 'merged.db')
    # get the workflow file counts
    #file_counts = get_block_workflow_file_count(block_workflows, file_counts)
    end16 = time.time()
    print('file counts', end16 - end15)

    # get release status of input sequences for each block
    # get the input limskeys for each workflow in project
    limskeys = get_workflow_limskeys(project_name, 'merged.db')
    end17 = time.time()
    print('limskeys', end17 - end16)

    # get the file swid and release status for each limskey for fastq-generating workflows
    # excluding fastq-import workflows
    status = get_file_release_status(project_name, 'merged.db')
    release_status = get_block_release_status(block_workflows, limskeys, status)
    end18 = time.time()
    print('release status', end18 - end17)

    # get the amount of data for each workflow
    amount_data = get_amount_data(project_name, 'merged.db')
    end19 = time.time()
    print('amount data', end19 - end18)
    
    # check if blocks are complete
    expected_workflows = sorted(['arriba', 'rsem', 'star', 'starfusion', 'mavis'])           
    complete = is_block_complete(blocks, expected_workflows)
    end20 = time.time()
    print('complete', end20 - end19)
    
    # order blocks based on the amount of data
    ordered_blocks = order_blocks(blocks, amount_data)
    end21 = time.time()
    print('order blocks', end21 - end20)

    # name each block according to the selected block order
    names = name_WT_blocks(ordered_blocks)
    end22 = time.time()
    print('name blocks', end22 - end21)    
    
    # get miso link
    miso_link = get_miso_sample_link(project_name, case, 'merged.db')
    end23 = time.time()
    print('miso link', end23 - end22)
    
    # sort sample pairs names
    sample_pairs_names = sorted(list(blocks.keys()))
    end24 = time.time()
    print('sort sample pairs', end24 - end23)
    
    
    print('all tasks', end24 -  start)

    
    
    return render_template('WT_case.html', project=project, routes = routes,
                           case=case, pipelines=pipelines, sample_pairs_names=sample_pairs_names,
                           blocks=blocks, names=names, ordered_blocks=ordered_blocks,
                           miso_link=miso_link, complete=complete, release_status=release_status,
                           block_date=block_date, workflow_date=workflow_date,
                           figures=figures, samples_star=samples_star, 
                           file_counts=file_counts, amount_data=amount_data)



@app.route('/download_wgs_block/<project_name>/<case>/<block>/<bmpp_parent>')
def download_block_data(project_name, case, block, bmpp_parent):
    '''
    
    
    '''
    
    # get the WGS blocks
    blocks = get_WGS_blocks_info(project_name, case, 'merged.db')
    # get the workflow names
    workflow_names = get_workflow_names(project_name, 'merged.db')
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
    star = get_star_case(project_name, case, 'novaseq', 'WT', 'merged.db')
    # get the tumor samples for each star id
    star_samples = map_samples_to_star_runs(project_name, star, 'merged.db')
    samples = get_WT_case_call_ready_samples(project_name, star_samples)
    # remove samples without analysis workflows
    D = map_workflows_to_samples(project_name, 'novaseq', samples, 'merged.db')
    # find the parents of each workflow
    parents = get_parent_workflows(project_name, 'merged.db')
    # get the parent workflows for each block
    parent_workflows = map_workflows_to_parent(D, parents)
    # find the blocks by mapping the analysis workflows to their parent workflows    
    blocks = find_WT_analysis_blocks(D, parents, parent_workflows, star)
    # map each workflow run id to its workflow name
    workflow_names = get_workflow_names(project_name, 'merged.db')
    # create json with workflow information for block for DARE
    block_data = create_block_json(project_name, blocks, block, star_parent, workflow_names)


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
    cases = get_cases(project_name, 'merged.db')
    # get library and sample counts
    counts = get_sample_counts(project_name, 'merged.db')
    # add missing donors to counts (ie, when counts are 0)
    counts = add_missing_donors(cases, counts)
    # count libraries for each library type
    # get the library types
    library_types =  get_library_types(project_name)
    libraries = count_libraries(project_name, library_types, cases, 'merged.db')
    
    D = {}
    for i in cases:
        donor = i['case_id']
        D[donor] = i
        D[donor]['normal'] = counts[donor]['normal']
        D[donor]['tumor'] = counts[donor]['tumor']
        for library_type in library_types:
            D[donor][library_type] = len(libraries[donor][library_type])
        
    data = pd.DataFrame(D.values())
    data.to_excel('{0}_cases.xlsx'.format(project_name), index=False)
   
    return send_file("{0}_cases.xlsx".format(project_name), as_attachment=True)



@app.route('/download_identifiers/<project_name>')
def download_identifiers_table(project_name):
    '''
    
    
    '''
    
    # get sequence file information
    files = collect_sequence_info(project_name, 'merged.db')
    # re-organize sequence information
    sequences = get_sequences(files)
    
    D = {}
    for i in sequences:
        d = {'Case': i['case'],
             'Donor': i['sample'],
             'SampleID': i['group_id'],
             'Sample': i['sample_id'],
             'Description': i['group_description'],
             'Library': i['library'],
             'Library Type': i['library_type'],
             'Tissue Origin': i['tissue_origin'],
             'Tissue Type': i['tissue_type'],
             'File Prefix': i['prefix']}
        D[i['case']] = d     
             
    data = pd.DataFrame(D.values())
     
    outputfile = '{0}_libraries.xlsx'.format(project_name)
    data.to_excel(outputfile, index=False)
    
    return send_file(outputfile, as_attachment=True)


# if __name__ == "__main__":
#     app.run(host='0.0.0.0')
    