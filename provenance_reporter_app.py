# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:32:40 2022

@author: rjovelin
"""

#import sqlite3
import json
from flask import Flask, render_template, request, url_for, flash, redirect, Response, send_file
#from werkzeug.exceptions import abort
import time
import pandas as pd
import matplotlib
matplotlib.use('agg')
from utilities import connect_to_db, get_miso_sample_link,\
    get_pipelines, get_workflow_names, get_library_design, secret_key_generator, \
    get_children_workflows
from whole_genome import get_call_ready_cases, get_amount_data, create_WG_block_json, \
    get_parent_workflows, get_workflows_analysis_date, get_workflow_file_count, \
    get_WGTS_blocks_info, get_sequencing_platform, get_selected_workflows, \
    review_wgs_blocks, get_case_workflows, update_wf_selection, get_block_counts, \
    get_wgs_blocks, create_WGS_project_block_json, get_workflow_output, get_release_status, \
    get_workflow_limskeys, get_file_release_status, map_fileswid_to_filename, \
    map_limskey_to_library, map_library_to_sample, get_WGS_standard_deliverables, \
    get_block_level_contamination, map_library_to_case, get_sample_sequencing_amount    
from whole_transcriptome import get_WT_call_ready_cases, get_WT_standard_deliverables, \
    create_WT_project_block_json, create_WT_block_json
from project import get_project_info, get_cases, get_sample_counts, count_libraries, \
     get_library_types, add_missing_donors, get_last_sequencing
from sequencing import get_sequences, collect_sequence_info, platform_name
from swg_ts import get_swg_ts, review_data, get_input_release_status, \
    create_swg_ts_sample_json, create_swg_ts_project_json, get_swg_ts_standard_deliverables, \
    order_workflows    
   


# map pipelines to views
routes = {'Whole Genome': 'whole_genome_sequencing',
          'Whole Transcriptome': 'whole_transcriptome',
          'Shallow Whole Genome': 'shallow_whole_genome',
          'Targeted Sequencing': 'targeted_sequencing'}


app = Flask(__name__)
#app.config['SECRET_KEY'] = secret_key_generator(10)
app.secret_key = secret_key_generator(10)



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
    
    database = 'merged.db'
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get case information
    cases = get_cases(project_name, database)
    # sort by donor id
    cases = sorted(cases, key=lambda d: d['case_id']) 
    # get the species
    species = ', '.join(sorted(list(set([i['species'] for i in cases]))))
    # get library and sample counts
    counts = get_sample_counts(project_name, database)
    # add missing donors to counts (ie, when counts are 0)
    counts = add_missing_donors(cases, counts)
    # count libraries for each library type
    # get the library types
    library_types =  get_library_types(project_name, database)
    # get the library names
    library_names = {i: get_library_design(i) for i in library_types}
    libraries = count_libraries(project_name, library_types, cases, database)
    # get the date of the last sequencing data
    seq_date = get_last_sequencing(project['project_id'], database)
    
    return render_template('project.html', routes=routes, project=project,
                           pipelines=pipelines, cases=cases, counts=counts,
                           seq_date=seq_date, species=species, libraries=libraries,
                           library_types = library_types, library_names=library_names)


@app.route('/<project_name>/sequencing', methods = ['GET', 'POST'])
def sequencing(project_name):
    
    database = 'merged.db'
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get sequence file information
    files = collect_sequence_info(project_name, database)
    # re-organize sequence information
    sequences = get_sequences(files)
    # sort data
    sequences.sort(key=lambda x: (x['case'], x['sample_id'], x['library'], x['platform']))
    # map the instrument short name to sequencing platform
    platform_names = platform_name(project_name, database)
 
    if request.method == 'POST':
        
        platforms = request.form.getlist('platform')
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
    
            # download all information if platforms are not selected
            if platforms:
                # check that platform is selected
                if platform_names[i['platform']] in platforms:
                    D[i['case']] = d     
            else:
                D[i['case']] = d
                    
        data = pd.DataFrame(D.values())
        outputfile = '{0}_libraries.xlsx'.format(project_name)
        data.to_excel(outputfile, index=False)
            
        return send_file(outputfile, as_attachment=True)

    else:
        return render_template('sequencing.html', routes=routes,
                           project=project, sequences=sequences,
                           pipelines=pipelines, platform_names=platform_names)



@app.route('/<project_name>/whole_genome_sequencing/', methods=['POST', 'GET'])
def whole_genome_sequencing(project_name):
    
    database = 'merged.db'
    workflow_db = 'workflows.db'
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get samples and libraries and workflow ids for each case
    cases = get_call_ready_cases(project_name, 'novaseq', 'WG', database)
    samples = sorted(list(cases.keys()))
    # get the block counts
    blocks = get_wgs_blocks(project_name, database, 'WGS_blocks')
    block_counts = get_block_counts(blocks)
       
    # get analysis block status
    # extract selected status of each workflow
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    block_status = review_wgs_blocks(blocks, selected)
    # make a list of donor ids with block status
    
    if request.method == 'POST':
        deliverable = request.form.get('deliverable')
        # get the workflow names
        workflow_names = get_workflow_names(project_name, database)
        if deliverable == 'selected':
            block_data = create_WGS_project_block_json(project_name, database, blocks, block_status, selected, workflow_names)
        elif deliverable == 'standard':
            # get the pipeline deliverables       
            deliverables = get_WGS_standard_deliverables()
            block_data = create_WGS_project_block_json(project_name, database, blocks, block_status, selected, workflow_names, deliverables)
        else:
            block_data = {}
                
        return Response(
            response=json.dumps(block_data),
            mimetype="application/json",
            status=200,
            headers={"Content-disposition": "attachment; filename={0}.WGS.json".format(project_name)})

    else:
        return render_template('Whole_Genome_Sequencing.html',
                           routes = routes,
                           project=project,
                           samples=samples,
                           cases=cases,
                           pipelines=pipelines,
                           block_status = block_status,
                           block_counts = block_counts
                           )
 

@app.route('/<project_name>/whole_genome_sequencing/<case>/<sample_pair>', methods = ['POST', 'GET'])
def wgs_case(project_name, case, sample_pair):
    
    database = 'merged.db'
    workflow_db = 'workflows.db'
    
    # get the number of lane sequence per sequence and platform and the corresponding release status
    lanes = get_sample_sequencing_amount(project_name, case, sample_pair, database,
                                         'Workflows', 'Workflow_Inputs', 'Libraries')
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get miso link
    miso_link = get_miso_sample_link(project_name, case, database)
    # get the WGS blocks
    blocks = get_WGTS_blocks_info(project_name, case, database, 'WGS_blocks')
    # sort sample pairs names
    sample_pairs_names = sorted(list(blocks.keys()))
    # get the workflow names
    workflow_names = get_workflow_names(project_name, database)
    # get the file count of each workflow in project
    file_counts = get_workflow_file_count(project_name, database)
    # get the amount of data for each workflow
    amount_data = get_amount_data(project_name, database)
    # get the creation date of all workflows
    creation_dates = get_workflows_analysis_date(project_name, database)
    # get the sequencing platform of each workflow
    platforms = get_sequencing_platform(project_name, database)
    # find the parents of each workflow
    parents = get_parent_workflows(project_name, database)
    # extract selected status of each workflow
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the contamination for each anchor workflow 
    contamination = get_block_level_contamination(project_name, database, blocks, sample_pair)
        
    if request.method == 'POST':
        # get the list of checked workflows        
        selected_workflows = request.form.getlist('workflow')
        # get the workflows of each block for sample pair and case
        case_workflows = get_case_workflows(case, database, 'WGS_blocks')
        # list all the workflows for a given sample pair
        # may include workflows from different blocks for a sample sample pair
        # this ensures blocks are mutually exclusive within a sample pair but not within a case
        workflows = []
        for i in case_workflows[sample_pair]:
            workflows.extend(case_workflows[sample_pair][i])
        update_wf_selection(workflows, selected_workflows, selected, workflow_db, 'Workflows')
        return redirect(url_for('wgs_case', case=case, project_name=project_name, sample_pair=sample_pair))
    else:
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
                           parents=parents,
                           selected = selected,
                           sample_pair=sample_pair,
                           contamination=contamination, 
                           lanes=lanes
                           )



@app.route('/<project_name>/<pipeline>/<case>/<sample_pair>/<workflow_id>')
def workflow(project_name, pipeline, case, sample_pair, workflow_id):
    
    database = 'merged.db'
    
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the workflow names
    workflow_names = get_workflow_names(project_name, database)
       
    # find the parents of each workflow
    parents = get_parent_workflows(project_name, database)
    if workflow_id in parents:
        parents = parents[workflow_id]
    else:
        parents = {}
    
    # find the children of each workflow
    D = get_children_workflows(project_name, database)
    children = {}
    if workflow_id in D:
        D = D[workflow_id]
        for i in D:
            if i['wf'] in children:
                children[i['wf']].append(i['children_id'])
            else:
                children[i['wf']] = [i['children_id']]
    
    # get the number of rows in table
    rows, parent_rows, children_rows = 0, 0, 0
    for i in parents:
        parent_rows += len(parents[i])
    for i in children:
        children_rows += len(children[i])
    rows = max([parent_rows, children_rows])
    if parent_rows > children_rows:
        parent_rows = rows
    elif parent_rows < children_rows:
        children_rows = rows
    
    # get input worflow sequences
    limskeys = get_workflow_limskeys(project_name, database, 'Workflow_Inputs')
    limskeys = limskeys[workflow_id]
    # get release status of input sequences, excluding fastq-import workflows
    D = get_file_release_status(project_name, database)
    sequence_status = {i:D[i] for i in limskeys}
    # map file swids to file names
    fastqs = map_fileswid_to_filename(project_name, database, 'Files')
    
    # map library to limskey
    all_libraries = map_limskey_to_library(project_name, database, table='Workflow_Inputs')
    libraries = all_libraries[workflow_id]
    
    
    # map libraries to samples
    all_samples = map_library_to_sample(project_name, database, table = 'Libraries')
    samples = all_samples[case]    
    
    
    sequences = []
    for i in limskeys:
        library = libraries[i]
        sample = samples[library]
        status = sequence_status[i]
        swid1, status1 = status[0][0], status[0][1]
        swid2, status2 = status[1][0], status[1][1]
        file1, file2 = fastqs[swid1], fastqs[swid2]
        seq = [[swid1, file1, status1], [swid2, file2, status2]]
        seq.sort(key=lambda x: x[1])
        for j in seq:
            sequences.append([sample, library, i, j[0], j[1], j[2]])
    sequences.sort(key=lambda x: x[0])
    
    # get workflow output files
    donors = map_library_to_case(project_name, database, table = 'Libraries')
    workflow_outputfiles = get_workflow_output(project_name, database, all_libraries, all_samples, donors, 'Files')
    files = workflow_outputfiles[workflow_id]
     
    
    # get the file release status
    release_status = get_release_status(project_name, database, 'FilesQC')
    
    return render_template('workflow.html',
                           project=project,
                           pipeline=pipeline,
                           case=case,
                           sample_pair=sample_pair,
                           parents=parents,
                           children=children,
                           parent_rows=parent_rows,
                           children_rows=children_rows,
                           rows=rows,
                           workflow_id=workflow_id,
                           workflow_names=workflow_names,
                           files=files,
                           release_status=release_status,
                           sequences=sequences
                           )


@app.route('/<project_name>/whole_transcriptome', methods = ['POST', 'GET'])
def whole_transcriptome(project_name):
    
    database = 'merged.db'
    workflow_db = 'workflows.db'
        
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get samples and libraries and workflow ids for each case
    cases = get_WT_call_ready_cases(project_name, 'novaseq', database, 'WT')
    samples = sorted(list(cases.keys()))
    # get the block counts
    blocks = get_wgs_blocks(project_name, database, 'WT_blocks')
    block_counts = get_block_counts(blocks)
    
    # get analysis block status
    # extract selected status of each workflow
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    block_status = review_wgs_blocks(blocks, selected)
    
    if request.method == 'POST':
        deliverable = request.form.get('deliverable')
        # get the workflow names
        workflow_names = get_workflow_names(project_name, database)
        
        if deliverable == 'selected':
            block_data = create_WT_project_block_json(project_name, database, blocks, block_status, selected, workflow_names)
        elif deliverable == 'standard':
            # get the pipeline deliverables       
            deliverables = get_WT_standard_deliverables()
            block_data = create_WT_project_block_json(project_name, database, blocks, block_status, selected, workflow_names, deliverables)
        else:
            block_data = {}
            
        return Response(
               response=json.dumps(block_data),
               mimetype="application/json",
               status=200,
               headers={"Content-disposition": "attachment; filename={0}.WT.json".format(project_name)})
    
    else:
        return render_template('Whole_transcriptome.html',
                         routes = routes, project=project,
                         samples=samples,
                         cases=cases,
                         pipelines=pipelines,
                         blocks=blocks,
                         block_counts=block_counts,
                         block_status=block_status)


@app.route('/<project_name>/whole_transcriptome/<case>/<tumor_sample>', methods=['POST', 'GET'])
def wt_case(project_name, case, tumor_sample):
        
    
    database = 'merged.db'
    workflow_db = 'workflows.db'
    
    # get the number of lane sequence per sequence and platform and the corresponding release status
    lanes = get_sample_sequencing_amount(project_name, case, tumor_sample, database,
                                         'Workflows', 'Workflow_Inputs', 'Libraries')

    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get miso link
    miso_link = get_miso_sample_link(project_name, case, database)
    # get the WT blocks
    blocks = get_WGTS_blocks_info(project_name, case, database, 'WT_blocks')
    # sort sample pairs names
    sample_names = sorted(list(blocks.keys()))
    # get the workflow names
    workflow_names = get_workflow_names(project_name, database)
    # get the file count of each workflow in project
    file_counts = get_workflow_file_count(project_name, database)
    # get the amount of data for each workflow
    amount_data = get_amount_data(project_name, database)
    # get the creation date of all workflows
    creation_dates = get_workflows_analysis_date(project_name, database)
    # find the parents of each workflow
    parents = get_parent_workflows(project_name, database)
    # extract selected status of each workflow
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    
    if request.method == 'POST':
        # get the list of checked workflows        
        selected_workflows = request.form.getlist('workflow')
        # get the workflows of each block and sample pair for case
        case_workflows = get_case_workflows(case, database, 'WT_blocks')
        # list all the workflows for a given sample that need an update
        # may include workflows from different blocks for a sample
        # this ensures blocks are mutually exclusive within a sample but not within a case
        workflows = []
        for i in case_workflows[tumor_sample]:
            workflows.extend(case_workflows[tumor_sample][i])
        update_wf_selection(workflows, selected_workflows, selected, workflow_db, 'Workflows')
        return redirect(url_for('wt_case', case=case, project_name=project_name, tumor_sample=tumor_sample))
    else:
        return render_template('WT_case.html',
                               project=project,
                               routes = routes,
                               case=case,
                               pipelines=pipelines,
                               blocks=blocks,
                               sample_names=sample_names,
                               workflow_names=workflow_names,
                               file_counts=file_counts,
                               amount_data=amount_data,
                               creation_dates=creation_dates,
                               parents=parents,
                               selected=selected,
                               tumor_sample=tumor_sample,
                               lanes=lanes
                               )


@app.route('/<project_name>/shallow_whole_genome/', methods=['POST', 'GET'])
def shallow_whole_genome(project_name):
        
    
    database = 'merged.db'
    workflow_db = 'workflows.db'
    
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get the shallow whole genome data
    swg = get_swg_ts(project_name, database, 'ichorcna', workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries')
    # get the selection status of workflows
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the input fastqs release status
    release_status = get_file_release_status(project_name, database)
    workflow_release_status = get_input_release_status(swg, release_status)
    status = review_data(swg, selected, workflow_release_status)
    
    row_counts = {}
    for i in swg:
        for j in swg[i]:
            row_counts[i] = len(swg[i][j])
    
    if request.method == 'POST':
        deliverable = request.form.get('deliverable')
        # get the workflow names
        workflow_names = get_workflow_names(project_name, database)
        if deliverable == 'selected':
            data = create_swg_ts_project_json(database, project_name, swg, workflow_names, selected)
        elif deliverable == 'standard':
            deliverables = get_swg_ts_standard_deliverables('swg')
            data = create_swg_ts_project_json(database, project_name, swg, workflow_names, selected, deliverables)
        else:
            data = {}
                         
        return Response(
                response=json.dumps(data),
                mimetype="application/json",
                status=200,
                headers={"Content-disposition": "attachment; filename={0}.SWG.json".format(project_name)})
    
    else:
        return render_template('shallow_whole_genome.html',
                               project=project,
                               routes = routes,
                               pipelines=pipelines,
                               swg=swg,
                               status=status,
                               row_counts=row_counts
                               )


@app.route('/<project_name>/shallow_whole_genome/<case>/<sample>', methods=['POST', 'GET'])
def swg_sample(project_name, case, sample):
    
    database = 'merged.db'
    workflow_db = 'workflows.db'

    # get the number of lane sequence per sequence and platform and the corresponding release status
    lanes = get_sample_sequencing_amount(project_name, case, sample, database,
                                         'Workflows', 'Workflow_Inputs', 'Libraries')

    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get miso link
    miso_link = get_miso_sample_link(project_name, case, database)
    # get the shallow whole genome data
    swg = get_swg_ts(project_name, database, 'ichorcna', workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries')
    # get the selection status of workflows
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the input fastqs release status
    release_status = get_file_release_status(project_name, database)
    status = get_input_release_status(swg, release_status)
    # get the workflow names
    workflow_names = get_workflow_names(project_name, database)
    file_counts = get_workflow_file_count(project_name, database)
    # get the amount of data for each workflow
    amount_data = get_amount_data(project_name, database)
    # get the creation date of all workflows
    creation_dates = get_workflows_analysis_date(project_name, database)
    # sort workflows according to amount of data, release status and creation date
    ordered_workflows = order_workflows(swg, amount_data, status, creation_dates)
      
    if request.method == 'POST':
        # get the selected workflow        
        selected_workflow = request.form.getlist('workflow')
        # get the workflows for the given sample
        workflows = list(swg[case][sample].keys())
        update_wf_selection(workflows, selected_workflow, selected, workflow_db, 'Workflows')
        return redirect(url_for('swg_sample', case=case, project_name=project_name, sample=sample))
    else:
        return render_template('SWG_sample.html',
                           project=project,
                           routes = routes,
                           pipelines=pipelines,
                           swg=swg,
                           status=status,
                           case=case,
                           sample=sample,
                           workflow_names=workflow_names,
                           file_counts=file_counts,
                           amount_data=amount_data,
                           creation_dates=creation_dates,
                           selected = selected,
                           ordered_workflows = ordered_workflows,
                           lanes=lanes
                           )




@app.route('/<project_name>/targeted_sequencing/', methods=['POST', 'GET'])
def targeted_sequencing(project_name):
    
    database = 'merged.db'
    workflow_db = 'workflows.db'
    
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get the consensus cruncher genome data
    CC = get_swg_ts(project_name, database, 'consensuscruncher', workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries')
    # get the selection status of workflows
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the input fastqs release status
    release_status = get_file_release_status(project_name, database)
    workflow_release_status = get_input_release_status(CC, release_status)
    status = review_data(CC, selected, workflow_release_status)
    
    row_counts = {}
    for i in CC:
        for j in CC[i]:
            row_counts[i] = len(CC[i][j])
    
    if request.method == 'POST':
        deliverable = request.form.get('deliverable')
        # get the workflow names
        workflow_names = get_workflow_names(project_name, database)
        if deliverable == 'selected':
            data = create_swg_ts_project_json(database, project_name, CC, workflow_names, selected)
        elif deliverable == 'standard':
            deliverables = get_swg_ts_standard_deliverables('ts')
            data = create_swg_ts_project_json(database, project_name, CC, workflow_names, selected, deliverables)
        else:
            data = {}
    
        return Response(
                response=json.dumps(data),
                mimetype="application/json",
                status=200,
                headers={"Content-disposition": "attachment; filename={0}.TS.json".format(project_name)})
    
    else:
        return render_template('targeted_sequencing.html',
                               project=project,
                               routes = routes,
                               pipelines=pipelines,
                               CC=CC,
                               status=status,
                               row_counts=row_counts
                               )


@app.route('/<project_name>/targeted_sequencing/<case>/<sample>', methods=['POST', 'GET'])
def TS_sample(project_name, case, sample):
    
    database = 'merged.db'
    workflow_db = 'workflows.db'
    
    # get the number of lane sequence per sequence and platform and the corresponding release status
    lanes = get_sample_sequencing_amount(project_name, case, sample, database,
                                         'Workflows', 'Workflow_Inputs', 'Libraries')
  
    # get the project info for project_name from db
    project = get_project_info(project_name, database)
    # get the pipelines from the library definitions in db
    pipelines = get_pipelines(project_name, database)
    # get miso link
    miso_link = get_miso_sample_link(project_name, case, database)
    # get the consensus cruncher data
    CC = get_swg_ts(project_name, database, 'consensuscruncher', workflow_table = 'Workflows', wf_input_table = 'Workflow_Inputs', library_table='Libraries')
    # get the selection status of workflows
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the input fastqs release status
    release_status = get_file_release_status(project_name, database)
    status = get_input_release_status(CC, release_status)
    # get the workflow names
    workflow_names = get_workflow_names(project_name, database)
    file_counts = get_workflow_file_count(project_name, database)
    # get the amount of data for each workflow
    amount_data = get_amount_data(project_name, database)
    # get the creation date of all workflows
    creation_dates = get_workflows_analysis_date(project_name, database)
    # sort workflows according to amount of data, release status and creation date
    ordered_workflows = order_workflows(CC, amount_data, status, creation_dates)
      
    if request.method == 'POST':
        # get the selected workflow        
        selected_workflow = request.form.getlist('workflow')
        # get the workflows for the given sample
        workflows = list(CC[case][sample].keys())
        update_wf_selection(workflows, selected_workflow, selected, workflow_db, 'Workflows')
        return redirect(url_for('TS_sample', case=case, project_name=project_name, sample=sample))
    else:
        return render_template('TS_sample.html',
                           project=project,
                           routes = routes,
                           pipelines=pipelines,
                           CC=CC,
                           status=status,
                           case=case,
                           sample=sample,
                           workflow_names=workflow_names,
                           file_counts=file_counts,
                           amount_data=amount_data,
                           creation_dates=creation_dates,
                           selected = selected,
                           ordered_workflows = ordered_workflows,
                           lanes=lanes
                           )



@app.route('/download_block/<project_name>/<pipeline>/<case>/<pair>/<anchor_wf>/<table>/<selection>')
def download_block_data(project_name, pipeline, case, pair, anchor_wf, table, selection):
        
    database = 'merged.db'
    workflow_db = 'workflows.db'
    
    
    # get the WGS blocks
    blocks = get_WGTS_blocks_info(project_name, case, database, table)
    # get the workflow names
    workflow_names = get_workflow_names(project_name, database)
    # get selected workflows
    selected_workflows = get_selected_workflows(project_name, workflow_db)
    # create json with workflow information for block for DARE
    #block_data = create_block_json(project_name, blocks, pair, anchor_wf, workflow_names, selected_workflows, selection)
    
    if pipeline == 'WG':
        block_data = create_WG_block_json(database, project_name, case, blocks, pair, anchor_wf, workflow_names, selected_workflows, selection)
    elif pipeline == 'WT':
        block_data = create_WT_block_json(database, project_name, case, blocks, pair, anchor_wf, workflow_names, selected_workflows, selection)
    
    pair_name = '.'.join(map(lambda x: x.strip(), pair.split('|')))
    
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(block_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}.{1}.{2}.{3}.{4}.{5}.json".format(project_name, pipeline, case, pair_name, anchor_wf, selection)})


@app.route('/download_swg/<project_name>/<datatype>/<case>/<sample>/<workflow_id>/<selection>')
def download_swg_ts_data(project_name, datatype, case, sample, workflow_id, selection):
        
    database = 'merged.db'
    workflow_db = 'workflows.db'
    
    # get the selection status of workflows
    selected_workflows = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the workflow names
    workflow_names = get_workflow_names(project_name, database)
    
    data = create_swg_ts_sample_json(datatype, database, project_name, case, sample, workflow_id, workflow_names, selected_workflows, selection)
            
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}.{1}.{2}.{3}.{4}.{5}.json".format(project_name, datatype.upper(), case, sample, workflow_id, selection)})



@app.route('/download_cases/<project_name>')
def download_cases_table(project_name):
    '''
    (str) -> None
    
    Download a table with project information in Excel format
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    database = 'merged.db'
    
    # get case information
    cases = get_cases(project_name, database)
    # get library and sample counts
    counts = get_sample_counts(project_name, database)
    # add missing donors to counts (ie, when counts are 0)
    counts = add_missing_donors(cases, counts)
    # count libraries for each library type
    # get the library types
    library_types =  get_library_types(project_name, database)
    libraries = count_libraries(project_name, library_types, cases, database)
    
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



# @app.route('/download_identifiers/<project_name>')
# def download_identifiers_table(project_name):
#     '''
    
    
#     '''
    
#     # get sequence file information
#     files = collect_sequence_info(project_name, 'merged.db')
#     # re-organize sequence information
#     sequences = get_sequences(files)
    
#     D = {}
#     for i in sequences:
#         d = {'Case': i['case'],
#              'Donor': i['sample'],
#              'SampleID': i['group_id'],
#              'Sample': i['sample_id'],
#              'Description': i['group_description'],
#              'Library': i['library'],
#              'Library Type': i['library_type'],
#              'Tissue Origin': i['tissue_origin'],
#              'Tissue Type': i['tissue_type'],
#              'File Prefix': i['prefix']}
#         D[i['case']] = d     
             
#     data = pd.DataFrame(D.values())
     
#     outputfile = '{0}_libraries.xlsx'.format(project_name)
#     data.to_excel(outputfile, index=False)
    
#     return send_file(outputfile, as_attachment=True)


# if __name__ == "__main__":
#     app.run(host='0.0.0.0')
    