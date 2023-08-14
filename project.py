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
    
    Returns a dictionary with library and sample counts for each donor of a project of interest
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    conn = connect_to_db()
    
    data = conn.execute("SELECT DISTINCT library, sample, tissue_type, group_id FROM Libraries WHERE project_id = '{0}';".format(project_name)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        donor = i['sample']
        library = i['library']
        if i['tissue_type'] == 'R':
            normal = i['sample'] + '_' + i['group_id']      
            tumor = ''
        else:
            normal = ''
            tumor = i['sample'] + '_' + i['group_id']
        if donor not in counts:
            counts[donor] = {}
        if 'library' not in counts[donor]:
            counts[donor]['library'] = set()
        if 'normal' not in counts[donor]:
            counts[donor]['normal'] = set()
        if 'tumor' not in counts[donor]:
            counts[donor]['tumor'] = set()
        counts[donor]['library'].add(library)
        if normal:
            counts[donor]['normal'].add(normal)
        elif tumor:
            counts[donor]['tumor'].add(tumor)


    for i in counts:
        counts[i]['library'] = len(counts[i]['library'])
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
            counts[i['case_id']] = {'library': 0, 'normal': 0, 'tumor': 0}
    return counts    

