# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:10:24 2023

@author: rjovelin
"""


import json
import os


def get_sequences(L):
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
    
    for i in range(len(L)):
        # keep only read1
        if json.loads(L[i]['attributes'])['read_number'] == '1':
            case = L[i]['sample']
            sample = L[i]['ext_id']
            library =  L[i]['library']
            library_type =  L[i]['library_type']
            tissue_origin =  L[i]['tissue_origin']
            tissue_type =  L[i]['tissue_type']
            group_id = L[i]['group_id']
            workflow = L[i]['workflow'] + '_' + L[i]['version']
            wfrun = L[i]['wfrun_id']
            file = L[i]['file']
            run = L[i]['run'] + '_' + str(L[i]['lane'])
            platform = L[i]['platform']
            status = L[i]['status']
            ticket = L[i]['ticket']
            read_count = json.loads(L[i]['attributes'])['read_count'] if 'read_count' in json.loads(L[i]['attributes']) else 'NA' 
            
            if 'GDR' in ticket:
                ticket = os.path.join('https://jira.oicr.on.ca/browse/', os.path.basename(ticket))    
            else:
                ticket = ''
            readcount = '{:,}'.format(int(read_count)) if read_count != 'NA' else 'NA'
            fileprefix = os.path.basename(file)
            fileprefix = '_'.join(fileprefix.split('_')[:-1])
            d = {'case': case, 'sample': sample, 'library': library, 'run': run,
                 'read_count': readcount, 'workflow': workflow, 'release_status': status,
                 'ticket': ticket, 'prefix':fileprefix, 'platform': platform,
                 'group_id': group_id, 'tissue_type': tissue_type, 'library_type': library_type,
                 'tissue_origin': tissue_origin}
            F.append(d)
       
    F.sort(key = lambda x: x['case'])
     
    return F
