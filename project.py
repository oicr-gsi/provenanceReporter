# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 21:11:22 2023

@author: rjovelin
"""


from utilities import connect_to_db


def get_project_info(project_name):
    '''
    (str) -> list
    
    Returns a list with project information extracted from database for project_name 
    
    Parameters
    ----------
    - project_name (str): Project of interest
    '''
    # connect to db
    conn = connect_to_db()
    # extract project info
    project = conn.execute('SELECT * FROM Projects WHERE project_id=\"{0}\"'.format(project_name)).fetchall()[0]
    conn.close()
    
    return project


def get_cases(project_name):
    '''
    (str) -> list
    
    Returns a list of dictionaries with case information
    
    Paramaters
    -----------
    - project_name (str): Project of interest
    '''
    
    conn = connect_to_db()
    data = conn.execute("SELECT DISTINCT case_id, donor_id, species, sex, created_date, modified_date, miso, parent_project FROM Samples WHERE project_id = '{0}'".format(project_name)).fetchall()
    
    data = [dict(i) for i in data]
       
    return data


def get_sample_counts(project_name):
    '''
    (str) - > dict
    
    Returns a dictionary with sample counts for each donor of a project of interest
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    conn = connect_to_db()
    
    data = conn.execute("SELECT DISTINCT sample, tissue_type, group_id FROM Libraries WHERE project_id = '{0}';".format(project_name)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        donor = i['sample']
        if i['tissue_type'] == 'R':
            normal = i['sample'] + '_' + i['group_id']      
            tumor = ''
        else:
            normal = ''
            tumor = i['sample'] + '_' + i['group_id']
        if donor not in counts:
            counts[donor] = {}
        if 'normal' not in counts[donor]:
            counts[donor]['normal'] = set()
        if 'tumor' not in counts[donor]:
            counts[donor]['tumor'] = set()
        if normal:
            counts[donor]['normal'].add(normal)
        elif tumor:
            counts[donor]['tumor'].add(tumor)


    for i in counts:
        counts[i]['normal'] = len(counts[i]['normal'])
        counts[i]['tumor'] = len(counts[i]['tumor'])
               
                    
    return counts            


def add_missing_donors(cases, counts):
    '''
    (list, dict) -> dict
    
    Update the sample counts with 0 values when donor_id found in cases is not already in counts
    
    Parameters
    ----------
    - cases (list): List of dictionary with case information
    - counts (dict): Dictionary with library and sample counts for each donor
    '''
        
    for i in cases:
        if i['case_id'] not in counts:
            counts[i['case_id']] = {'normal': 0, 'tumor': 0}
    return counts    


def get_library_types(project_name):
    '''
    (str) -> list
    
    Returns a list of different library types for a given project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    # connect to db
    conn = connect_to_db()
    # extract library types
    data = conn.execute("SELECT DISTINCT library_types FROM Projects WHERE project_id = '{0}';".format(project_name)).fetchall()
    conn.close()
    
    library_types = sorted(list(map(lambda x: x.strip(), data[0]['library_types'].split(','))))
   
    return library_types


def count_libraries(project_name, library_types, cases):
    '''
    (str, list, list) -> dict
    
    Returns a dictionary with libraries for each library type and sample for a given project
       
    Parameters
    ----------
    - project_name (str) Name of the project of interest
    - library_types (list): List of library types recorded for project
    - cases (list): List of dictionary with case information  
    '''
    
    # connect to db
    conn = connect_to_db()
    # extract library source
    data = conn.execute("SELECT DISTINCT sample, library_type, library FROM Libraries WHERE project_id = '{0}';".format(project_name)).fetchall()
    conn.close()
    
    libraries= {}
    
    # initiate the libraries dict
    for i in cases:
        libraries[i['case_id']] = {}
        for j in library_types:
            libraries[i['case_id']][j] = set()
    
    # record libraries for each library type
    for i in data:
        libraries[i['sample']][i['library_type']].add(i['library'])
    
    return libraries



def get_last_sequencing(project_name):
    '''
    (str) -> str
    
    Returns the date of the last sequencing for the project of interest
    
    Paramaters
    ----------
    - project_name (str): Project of interest
    '''
    
    conn = connect_to_db()
    sequencing = conn.execute("SELECT DISTINCT Workflow_Inputs.run FROM Workflow_Inputs JOIN Files \
                              WHERE Workflow_Inputs.project_id = '{0}' AND Files.project_id = '{0}' \
                              AND Files.wfrun_id = Workflow_Inputs.wfrun_id \
                              AND LOWER(Files.workflow) in ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq');".format(project_name)).fetchall()
    conn.close()
    
    # get the most recent creation date of fastq generating workflows
    
    if sequencing:
        sequencing = list(set(sequencing))
        sequencing = [i['run'] for i in sequencing]
        sequencing = map(lambda x: x.split('_'), sequencing)
        seq_dates = [i for i in sequencing if any(list(map(lambda x: x.isdigit(), i)))]
        
        F = lambda y: list(map(lambda x: x.isdigit(), y)).index(True)
        date_index = list(map(lambda x: F(x), seq_dates))
        seq_dates = [seq_dates[i][date_index[i]] for i in range(len(date_index))]
        seq_date = sorted(list(set(seq_dates)))[-1]
        seq_date = '20' + str(seq_date)[:2] + '-' + str(seq_date)[2:4] + '-' + str(seq_date)[4:]
        
    else:
        seq_date = 'NA'
    return seq_date
