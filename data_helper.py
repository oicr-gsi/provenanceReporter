# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 14:47:15 2025

@author: rjovelin
"""

import json
from commons import list_case_workflows
import itertools



def load_data(provenance_data_file):
    '''
    (str) -> list
    
    Returns the list of data contained in the provenance_data_file
    
    Parameters
    ----------
    - provenance_data_file (str): Path to the file with production data extracted from Shesmu
    '''

    infile = open(provenance_data_file, encoding='utf-8')
    provenance_data = json.load(infile)
    infile.close()
    
    return provenance_data


def clean_up_data(provenance_data):
    '''
    (list) -> list
    
    Returns the list of case information removing any case for which information is missing
            
    Parameters
    ----------
    - provenance_data (list): List of dictionaries, each representing the data of a single case
    '''
    
    to_remove = [i for i in provenance_data if len(i['project_info']) == 0]
    for i in to_remove:
        provenance_data.remove(i)
    
    return provenance_data


def is_case_info_complete(case_data):
    '''
    (dict) -> bool
    
    Returns True if the case information is complete
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''
    
    complete = True
    for i in case_data:
        if len(case_data[i]) == 0:
            complete = False
            break
    
    return complete
    
    
def is_signoff_complete(case_data):
    '''
    (dict) -> bool

    Returns True if signoff is complete for all the lims Ids that pass QC
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''

    
    # evaluate all lims ids --> signoff is indicated by complete sequencing status   
    seq = json.loads(case_data['case_info']['sequencing'])
    # evaluate only lims ids for sequencing
    complete = []
    for i in seq:
        if i['type'] == 'FULL_DEPTH_SEQUENCING':
            complete.append(i['complete'])
    
    return all(complete)        
        

def map_lims_to_tests(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all sequencing lims ids passing QC for
    each tests in the case. 
    Assumption: sign off is complete (ie. sequencing status complete for all limds)
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''
    
    D = {}
        
    seq = json.loads(case_data['case_info']['sequencing'])
    # evaluate only lims ids for sequencing
    for i in seq:
        if i['type'] == 'FULL_DEPTH_SEQUENCING':
            # assumes signoff is complete
            assert i['complete']
            test = i['test']
            for j in i['limsIds']:
                if j['qcFailed'] == False:
                    limsid = j['id']
                    if test not in D:
                        D[test] = [limsid]
                    else:
                        D[test].append(limsid)
                        D[test].sort()
    return D
   



def list_library_qualif_lims(case_data):
    '''
    
    
    
    '''

    L = []
    
    
    
    seq = json.loads(case_data['case_info']['sequencing'])
    
    # evaluate only lims ids for sequencing
    for i in seq:
        if i['type'] == 'LIBRARY_QUALIFICATION':
            for j in i['limsIds']:
                L.append(j['id'])
            
    return L



def map_samples_to_lims(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all lims id for each sample in a case    
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''


    # list lims used in library qualification
    qualif = list_library_qualif_lims(case_data)
    
    D = {}
        
    for i in case_data['sample_info']:
        sampleid = i['sampleId']
        limsid = i['limsId']
        # do not include lims for library qualification
        if limsid not in qualif:
            if sampleid not in D:
                D[sampleid] = [limsid]
            else:
                D[sampleid].append(limsid)
                D[sampleid].sort()
            
    return D
    

   
def map_workflows_to_lims(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with lims ids for each workflow in a case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''
        
    D = {}
        
    for d in case_data['workflow_runs']:
        wfrunid = d['wfrunid']
        limsids = d['limsIds'].split(',')
        if wfrunid not in D:
            D[wfrunid] = limsids
        else:
            D[wfrunid].extend(limsids)
        D[wfrunid] = list(set(D[wfrunid]))

    return D




def map_expected_workflows_to_runid(workflow_info, case_workflows):
    '''
    (dict, list) -> bool
    
    Returns a dictionary with the workflow run ids, if they exist, for each expected
    workflow from the assay of a case
        
    Parameters
    ----------
    - workflow_info (dict): Dictionary mapping worfklow run ids to workflow names
    - case_workflows (list): List of expected workflows from the assay
    '''
    
    # map the workflow run ids of production workflows to the expected workflows
    D = {}
    for workflow in case_workflows:
        D[workflow] = []
        for wfrunid in workflow_info:
            if workflow_info[wfrunid] == workflow:
                D[workflow].append(wfrunid)
    
    return D    
    

def get_assay_expected_workflows(pipeline_workflows, tests_samples, samples_lims):
    '''
    (dict, dict, dict) -> dict
    
    Returns a dictionary with expected workflows and corresponding lims according to the defined assays, pipeline
    and case information
    
    Parameters
    ----------
    - pipeline_workflows (dict): Dictionary with pipeline workflows
    - tests_samples (dict): Dictionary with all the samples mapping each test
    - samples_limns (dict): Dictionary with the lims mapping each assay 
    '''
       
    workflows = {}

    # loop over workflows in pipeline
    for workflow in pipeline_workflows:
        # get the expected tests for each workflow from the assay definition
        level = pipeline_workflows[workflow]['level']
        # collect the expected lims for each workflow depending on the case information
        if level == 'lane':
            # each lims id of each test should have a separate workflow run id
            tests = pipeline_workflows[workflow]['tests'].split(',')
            # get the parent workflows
            parents = pipeline_workflows[workflow]['parent_workflows']
            if parents:
                parents = list(map(lambda x: x.strip(), parents.split(',')))
            # collect the limsids for each test using case data
            for test in tests:
                # get the sample id - each test can have multiple samples
                for sampleid in tests_samples[test]:
                    # get the corresponding lims 
                    limsids = samples_lims[sampleid]
                    for lims in limsids:
                        if workflow in workflows:
                            workflows[workflow].append({'workflow': workflow,'test': [test], 'sampleid': sampleid, 'limsids': lims, 'parents': parents, 'parent_workflows': []})
                        else:
                            workflows[workflow] = [{'workflow': workflow, 'test': [test], 'sampleid': sampleid, 'limsids': lims, 'parents': parents, 'parent_workflows': []}]
        
        elif level == 'merge':
            # each workflow has a set of lims ids combined from each corresponding test
            
            
            ####### NEED TO IMPLEMENT THIS - WHICH WORKFLOW HAS 
            
            
            ## testIds can be both , and | separated
            ## if the tests are , separated, then each is a test_set, which will have their own workflow run
            ## within a test_set, testIds are | separated, then each workflow will have limids from all of the test, FOR each test with that

            
            # test_sets=assay[pipeline][wf]['testIds'].split(",")
            # #print("test_sets:",test_sets)
            # for test_set in test_sets:
            #     ### each test in the test_Set might be represented multiple times in the case, and each needs to be accounted fror in combing the data
            #     TestIds=test_set.split("|")


            
            

            if ',' in pipeline_workflows[workflow]['tests']:
                tests = pipeline_workflows[workflow]['tests'].split(',')  
                # get the parent workflows
                parents = pipeline_workflows[workflow]['parent_workflows']
                if parents:
                    parents = list(map(lambda x: x.strip(), parents.split(',')))
                # collect the limsids for each test using case data
                for test in tests:
                    # get the sample id - each test can have multiple samples
                    for sampleid in tests_samples[test]:
                        # get the corresponding lims 
                        limsids = samples_lims[sampleid]
                        # the workflow has all the lims
                        limsids = ','.join(sorted(list(limsids)))
                        if workflow in workflows:
                            workflows[workflow].append({'workflow': workflow, 'test': [test], 'sampleid': sampleid, 'limsids': limsids, 'parents': parents, 'parent_workflows': []})
                        else:
                            workflows[workflow] = [{'workflow': workflow, 'test': [test], 'sampleid': sampleid, 'limsids': limsids, 'parents': parents, 'parent_workflows': []}]
                

            elif '|' in pipeline_workflows[workflow]['tests']:
                tests = pipeline_workflows[workflow]['tests'].split('|')
                # get the parent workflows
                parents = pipeline_workflows[workflow]['parent_workflows']
                if parents:
                    parents = list(map(lambda x: x.strip(), parents.split(',')))
                # collect the limsids for each test using case data
                # get the expected combinations of lims for each combination of test samples
                L = []
                S = []
                for test in tests:
                    l = []
                    s = []
                    # get the sample id - each test can have multiple samples
                    for sampleid in tests_samples[test]:
                        # get the corresponding lims 
                        limsids = samples_lims[sampleid]
                        l.append(limsids)
                        s.append(sampleid)
                    L.append(l)
                    S.append(s)
                
                combined_lims = list(itertools.product(*L))
                combined_samples = list(itertools.product(*S))
                
                
                # merge and sort each set of lims for each set of combined tests 
                for i in range(len(combined_lims)):
                    merged_lims = []
                    merged_samples = []
                    for j in combined_lims[i]:
                        merged_lims.extend(j)
                    for k in combined_samples[i]:
                        merged_samples.append(k)
                    
                    merged_lims = ','.join(sorted(merged_lims))
                    merged_samples = ','.join(sorted(merged_samples))
                               
                    if workflow in workflows:
                        workflows[workflow].append({'workflow': workflow, 'test': tests, 'sampleid': merged_samples, 'limsids': merged_lims, 'parents': parents, 'parent_workflows': []})
                    else:
                        workflows[workflow] = [{'workflow': workflow, 'test': tests, 'sampleid': merged_samples, 'limsids': merged_lims, 'parents': parents, 'parent_workflows': []}]
           
    # add parent workflow information 
    for workflow in workflows:
        for d in workflows[workflow]:
            if d['parents']:
                # find the corresponding parents
                for parent in d['parents']:
                    for k in workflows[parent]:
                        # check that all parent samples are in the children samples
                        # check that all parent tests are in the children tests
                        # check that all parent lims are in the children lims
                        if set(k['sampleid'].split(',')).issubset(set(d['sampleid'].split(','))) and \
                           set(k['test']).issubset(set(d['test'])) and \
                           set(k['limsids'].split(',')).issubset(set(d['limsids'].split(','))):
                           d['parent_workflows'].append(k)       
                        
    return workflows        
    

def get_production_workflows(samples_workflows, workflow_lims):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary mapping each workflow run id, their lims and sample to each workflow
    
    Parameters
    ----------
    - samples_workflows (dict): Dictionary mapping all lims for each sample
    - workflow_lims (dict): Dictionary mapping the lims to each workflow
    '''        
        
    # reorganize data: {workflow: {{'wfrunid': ,'limsids':, 'samples':}}}    
        
    D = {}
    
    for wfrunid in workflow_lims:
        lims = ','.join(sorted(workflow_lims[wfrunid]))
        samples = []
        names = []
        for sample in samples_workflows:
            for d in samples_workflows[sample]['workflows']:
                if d['wfrun_id'] == wfrunid:
                    samples.append(sample)
                    names.append(d['workflow'])
                    break
        names = list(set(names)) 
        assert len(names) == 1
        name = names[0]
        samples = ','.join(sorted(samples))
        if name not in D:
            D[name] = {}
        D[name][wfrunid] = {'limsids': lims, 'samples': samples}
            
    return D        
     
    
    
def find_production_workflow(production_workflows, d):
    '''
    (dict, dict) -> dict    
    
    Returns a dictionary with prodcution data mapping the expected data for a specific workflow
                
    Parameters
    ----------
    - production_workflows (dict): Dictionary with case data extracted from the provenance reporter
    - d (dict): Dictionary with expected workflow information based on assay and case info
    '''
    
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']

    data = {'workflow': None, 'limsids': None, 'wfrunid': None, 'tests': None, 'samples': None, 'parents': []}
    workflow = d['workflow']
    expected_lims = d['limsids']
    expected_samples = d['sampleid']
    test = d['test']
    #  find the workflow in production with the expected limsids and samples
    if workflow not in sequencing_workflows:
        if workflow in production_workflows:
            for wfrunid in production_workflows[workflow]:
                limsids = production_workflows[workflow][wfrunid]['limsids']
                samples = production_workflows[workflow][wfrunid]['samples']
                if expected_samples == samples and expected_lims == limsids:
                    ### check that only 1 wfrunids match the requirement
                    assert data['wfrunid'] is None 
                    # update data collector
                    data['limsids'] = limsids
                    data['samples'] = samples
                    data['wfrunid'] = wfrunid
                    data['tests'] = test
                    data['workflow'] = workflow
     
    else:
        # find the actual sequencing workflow as it may differ from assay
        for key in sequencing_workflows:
            if key in production_workflows:
                for wfrunid in production_workflows[key]:
                    limsids = production_workflows[key][wfrunid]['limsids']
                    samples = production_workflows[key][wfrunid]['samples']
                    if expected_samples == samples and expected_lims == limsids:
                        ### check that only 1 wfrunids match the requirement
                        assert data['wfrunid'] is None 
                        # update data collector
                        data['limsids'] = limsids
                        data['samples'] = samples
                        data['wfrunid'] = wfrunid
                        data['tests'] = test
                        data['workflow'] = key
    return data         
    





    
def map_expected_production_workflows(expected_workflow_lims, production_workflows):
    '''
    (dict, dict) -> dict    
    
    Returns a dictionary with prodcution data mapping the expected data from the assay ans case info
    with the production data available for a case in the provenance reporter
            
    Parameters
    ----------
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    - production_workflows (dict): Dictionary with case data extracted from the provenance reporter
    '''

    D = {}
        
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
        
    # for sequencing workflows, the assay may indicate bcl2fastq but the 
    # sequencing workflows may be diferent if data is injected
        
    for workflow in expected_workflow_lims:
        D[workflow] = []
        for d in expected_workflow_lims[workflow]:
            data = find_production_workflow(production_workflows, d)
            # find the parent workflows
            if d['parents']:
                for k in d['parent_workflows']:
                    # map expected parents to workflows in production
                    parent_workflows = find_production_workflow(production_workflows, k)
                    # update parent if parent workflow found                                    
                    if parent_workflows['wfrunid'] is not None:
                        parent_workflows = {i:j for i,j in parent_workflows.items() if i != 'parents'}
                        data['parents'].append(parent_workflows)
            D[workflow].append(data)
            
            
    return D            
     





       
            
def is_incomplete_workflow_run(d):
    '''
    (dict) -> bool
    
    Returns True is all keys in d have been populated with data
    
    Parameters
    ----------
    - d (dict): Dictionary with workflow run id information in case_analysis
    '''
        
    # analysis is incomplete if any workflow in assay has missing information
    return any(map(lambda x: x is None or len(x) == 0, list(d.values())))
    
    
def is_incomplete_sequencing_workflow_run(d):
    '''
    (dict) -> bool
    
    Returns True is all keys in d have been populated with data
    
    Parameters
    ----------
    - d (dict): Dictionary with workflow run id information in case_analysis
    '''
    
    # analysis is incomplete if any workflow in assay has missing information expect parents
    vals = [d[i] for i in d.keys() if i != 'parents']
    return any(map(lambda x: x is None or len(x) == 0, vals))
 
    
def is_data_complete(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> bool    
    
    Returns True if each workflow in case_analysis has complete information
        
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''
        
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
            
    complete = True
        
    if cases_analysis.keys() != expected_workflow_lims.keys():
        complete = False
    
    # check if there are extra workflows
    for workflow in cases_analysis:
        if len(cases_analysis[workflow]) < len(expected_workflow_lims[workflow]):
            complete = False
    
    # check that all workflows have been identified
    for workflow in cases_analysis:
        # check if workflow if sequencing workflow (not expecting parents)
        if workflow not in sequencing_workflows:
            for d in cases_analysis[workflow]:
                # analysis is incomplete if any workflow in assay has missing information
                if is_incomplete_workflow_run(d):
                    complete = False
                # check if parents are defined
                for parent in d['parents']:
                    if is_incomplete_workflow_run(parent):
                        complete = False
        else:
            for d in cases_analysis[workflow]:
                if is_incomplete_sequencing_workflow_run(d):
                    complete = False
                
    return complete
                



def identify_workflows_with_missing_data(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> list

    Returns a list of workflows with missing data
            
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''

    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
            
    missing = [workflow for workflow in cases_analysis if workflow not in expected_workflow_lims]
        
    # check if there are missing iterations
    for workflow in cases_analysis:
        if workflow in expected_workflow_lims:
            if len(cases_analysis[workflow]) < len(expected_workflow_lims[workflow]):
                missing.append(workflow)
                
    # check that all workflows have been identified
    for workflow in cases_analysis:
        # check if workflow if sequencing workflow (not expecting parents)
        if workflow not in sequencing_workflows:
            for d in cases_analysis[workflow]:
                # analysis is incomplete if any workflow in assay has missing information
                if is_incomplete_workflow_run(d):
                    missing.append(workflow)
                # check if parents are defined
                for parent in d['parents']:
                    if is_incomplete_workflow_run(parent):
                        missing.append(workflow)
        else:
            for d in cases_analysis[workflow]:
                if is_incomplete_sequencing_workflow_run(d):
                    missing.append(workflow)
                
                    
    missing = list(set(missing))

    return missing                    


    
    
def no_extra_data(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> bool    
    
    Returns True is each workflow in case_analysis have a single workflow run id
    matching the lims requirements
        
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''
    
    no_extra = True
        
    # check if there are extra workflows
    for workflow in cases_analysis:
        if len(cases_analysis[workflow]) > len(expected_workflow_lims[workflow]):
            no_extra = False
    
    return no_extra
    
       
    
def identify_extra_workflows(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> list    
    
    Returns a list of workflow with multiple run ids matching the lims requirements
    
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''
    
    no_extra = []
       
    # check if there are extra workflows
    for workflow in cases_analysis:
        if len(cases_analysis[workflow]) > len(expected_workflow_lims[workflow]):
            extra.append(workflow)
    
    extra = list(set(extra))
    
    return extra
    
    

def complete_expected_workflows(workflow_info, case_workflows):
    '''
    (dict, list) -> bool
    
    Returns True if all the expected workflows in case workflows have run in 
    in production and have assigned workflow run ids
    
    
    Parameters
    ----------
    - workflow_info (dict): Dictionary mapping worfklow run ids to workflow names
    - case_workflows (list): List of expected workflows from the assay
    '''

    # map the workflow run ids of production workflows to the expected workflows
    expected_workflows = map_expected_workflows_to_runid(workflow_info, case_workflows)

    complete = True
    for workflow in expected_workflows:
        if len(expected_workflows[workflow]) == 0:
            complete = False
    
    return complete
   


def identify_missing_workflows(workflow_info, case_workflows):
    '''
    (dict, list) -> bool
    
    Returns True if all the expected workflows in case workflows have run in 
    in production and have assigned workflow run ids
        
    Parameters
    ----------
    - workflow_info (dict): Dictionary mapping worfklow run ids to workflow names
    - case_workflows (list): List of expected workflows from the assay
    '''

    # map the workflow run ids of production workflows to the expected workflows
    expected_workflows = map_expected_workflows_to_runid(workflow_info, case_workflows)

    missing = []
    for workflow in expected_workflows:
        if len(expected_workflows[workflow]) == 0:
            missing.append(workflow)
    
    missing = list(set(missing))    
    
    return missing
    
  

def map_tests_to_samples(case_data, tests):
    '''
    (dict, dict) -> dict

    Returns a dictionary matching all samples for each test
    Assumption: sign off is complete (ie. sequencing status complete for all limds)
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    - tests (dict): Dictionary with lims ids for each tests
    '''

    D = {}

    for test in tests:
        # find the sample corresponding to each lims id
        for i in case_data['sample_info']:
            limsid = i['limsId']
            sampleid = i['sampleId']
            if limsid in tests[test]:
                if test not in D:
                    D[test] = [sampleid]
                else:
                    D[test].append(sampleid)
                D[test] = sorted(list(set(D[test]))) 

    return D



def sort_lims_by_samples(tests_samples, samples_lims):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with lists of lims for each sample, if multiple samples
    exist, for each test in a case
        
    Parameters
    ----------
    - test_samples (dict): Dictionart mapping tests with their samples
    - samples_lims (dict): Dictionary mapping samples with their lims ids
    '''
    
    D = {}
    
    for test in tests_samples:
        for sample in tests_samples[test]:
            limsids = samples_lims[sample]
            if test not in D:
                D[test] = [limsids]
            else:
                D[test].append(limsids)
     
    return D



def clean_up_workflows(case_data):
    '''
    (dict) -> dict    
    
    Remove children and parent workflows of case workflows that are not in a case
    
    Parameters
    ----------
    case_data (dict): Dictionary with production data of a given case
    '''

    # make a list of case workflows
    workflows = list_case_workflows(case_data)
    
    # evaluate all children and parents workflows
    for k in ['children', 'parents']:
        for i in range(len(case_data['workflow_runs'])):
            L = json.loads(case_data['workflow_runs'][i][k])
            to_remove = [j for j in L if j[1] not in workflows]
            if len(to_remove) >= 1:
                for j in to_remove:
                    assert j in L
                    L.remove(j)
            case_data['workflow_runs'][i][k] = json.dumps(L)            

    return case_data



def check_workflow_relationships(cases_analysis, parent_workflows, workflow_info):
    '''
    (dict, dict, dict) -> bool
    
    Returns True if all the identified workflows and their parents in case_analysis
    have indded a parent-child workflow relationship
    
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - parent_workflows (dict): Dictionary with parent-children workflow relationships
    - workflow_info (dict): Dictionary with workflow names mapped to workflow run ids
    '''
        
    correct = True
       
    for workflow in cases_analysis:
        for d in cases_analysis[workflow]:
            # check parent if the exist
            wfrunid = d['wfrunid']
            if d['parents']:
                for k in d['parents']:
                    if k['workflow'] != workflow_info[k['wfrunid']]:
                         correct = False
                    if wfrunid not in parent_workflows[k['wfrunid']]:
                        correct = False
                        
    return correct


def reformat_pipeline_workflows(L):
    '''
    
    
    '''
    
    D = {}
    
    
    for d in L:
        workflow = d['workflows']
        assert workflow not in D
        D[workflow] = d
    
    return D
    
