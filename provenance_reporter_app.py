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





# def extract_info():
#     '''
    
    
#     '''
    






# my %opts = ();
# GetOptions(
#         "project|p=s"    => \$opts{"project"},
#         "fpr|f:s"        => \$opts{"fpr"},
#         "workflows|w:s"  => \$opts{"workfkows"},
#         "pipelines|l:s"  => \$opts{"pipelines"},
#         "out|o:s"        => \$opts{"out"},
#         "qc|q"           => \$opts{"qc"},     ## boolean options
#         "release|r"      => \$opts{"release"} ## boolean options
# );
# validate_options(\%opts);
# use JSON::PP;

# my $nabu_api="http://gsi-dcc.oicr.on.ca:3000";




# if($project eq "all"){
#         %{$provenance{pipeline}}=load_fpr($project,$rpt,\%nabu);
# }else{

#         print STDERR "loading pinery info from $pinery_api for project:$project\n";
#         %{$provenance{lims}}=load_pinery($project,$pinery_api);
#         ### bypassing nabu load due to problems with nabu
#         ### empty hash should still work, but no nabu info
#         print STDERR "loading nabu info from $nabu_api for project:$project\n";
#         %nabu=load_nabu($project,$nabu_api);

#         print STDERR "loading file provenance from $rpt and collating information\n";

#         %{$provenance{pipeline}}=load_fpr($project,$rpt,\%nabu);
# }

# #$provenance{pipeline}{files}=load_nabu($project,$nabu_api,$provenance{pipeline}{files});
# #print Dumper($provenance{lims});exit;


# print STDERR "generating provenance json for project $project\n";
# my $json=encode_json(\%provenance);
# (open my $JSON,">","$project.json");
# print $JSON $json;
# close $JSON;



# sub load_fpr{
#         my($project,$fpr,$nabu)=@_;

#         my @recs=`zcat $fpr`;chomp @recs;
#         ### FIRST PASS of FPR
#         my %files;   ### informaitn on files, by fileid
#         my %workflows;  ### inofmraiton on workflows, by workflow run id
#         my %libraries;    ### information on libraries, by library
#         my %data;

#         my $rec_count=0;
#         for my $rec(@recs){

#                 my @f=split /\t/,$rec;
#                 unshift(@f,1);   ### shift index to column number


#                 next if( ($project ne "all") && ($f[2] ne $project));   ### only use record from this project
#                 $rec_count++;

#                 my $wfrunid="wf.".$f[37];
#                 my $fid="f.".$f[45];

#                 $workflows{$wfrunid}{info}{wf}=$f[31];
#                 $workflows{$wfrunid}{info}{wfv}=$f[32];
#                 $workflows{$wfrunid}{info}{wfinput_string}=$f[39];

#                 $workflows{$wfrunid}{files}{$fid}=1;

#                 ### add info about the objects sequenced
#                 my $limskey=$f[57];
#                 $workflows{$wfrunid}{libraries}{$limskey}=1;

#                 $libraries{$limskey}{id}=$f[8];
#                 $libraries{$limskey}{lib}=$f[14];
#                 $libraries{$limskey}{run}=$f[19];
#                 $libraries{$limskey}{lane}=$f[25];

#                 $files{$fid}={path=>$f[47],wf=>$f[31],wfv=>$f[32],wfrunid=>$wfrunid,md5sum=>$f[48]};

#                 ### additional info about the run
#                 my %info;
#                 map{
#                         my($key,$value)=split /=/,$_;
#                         $info{$key}=$value;
#                 }split /;/,$f[38];
#                 map{
#                         my($key,$value)=split /=/,$_;
#                         $info{$key}=$value;
#                 }split /;/,$f[46];

#                 ### collecting the following fields
#                 map{
#                         if($info{$_}){
#                                 $workflows{$wfrunid}{info}{$_}=$info{$_};
#                         }
#                 }qw/mode reference read_count read_number/;
#                 #print $f[31],"\n",Dumper(%info);<STDIN>;


#         }
#         print STDERR "FPR : $rec_count records\n";

#         #print Dumper(%workflows);<STDIN>;

#         ### add children and parentsto the data
#         for my $wfrunid(sort keys %workflows){
#                 my $wf=$workflows{$wfrunid}{info}{wf};
#                 my @input_fids=split /;/,$workflows{$wfrunid}{info}{wfinput_string};
#                 for my $input_fid(@input_fids){
#                         $input_fid="f." . $input_fid;
#                         my $input_wfrun_id=$files{$input_fid}{wfrunid};
#                         #print "$input_wfrun_id";<STDIN>;
#                         my $input_wf=$files{$input_fid}{wf};
#                         ### children
#                         $workflows{$input_wfrun_id}{children}{workflows}{$wfrunid}=1;
#                         ### parents
#                         $workflows{$wfrunid}{parents}{files}{$input_fid}=1;
#                         my $input_wfrunid=$files{$input_fid}{wfrunid};
#                         $workflows{$wfrunid}{parents}{workflows}{$input_wfrunid}=1;
#                 }
#         }
#         #print Dumper(%workflows);<STDIN>;

#         ### organize provenance hash for json import
#         ### workflows first
#         for my $wfrunid(sort keys %workflows){

#                 #print "wfrunid=$wfrunid";<STDIN>;
#                 #print Dumper($workflows{$wfrunid});<STDIN>;

#                 my %workflow=(wfrunid=>$wfrunid);
#                 map{ $workflow{info}{$_}=$workflows{$wfrunid}{info}{$_} }keys %{$workflows{$wfrunid}{info}};

#                 #print Dumper(%workflow);<STDIN>;
#                 my @fids=sort keys %{$workflows{$wfrunid}{files}};
#                 for my $fid(@fids){
#                         my %file=(fid=>$fid);
#                         map{ $file{$_}=$files{$fid}{$_} } keys %{$files{$fid}};
#                         push(@{$workflow{files}},\%file);
#                 }
#                 #print Dumper(%workflow);<STDIN>;

#                 my @limskeys=sort keys %{$workflows{$wfrunid}{libraries}};
#                 for my $limskey(@limskeys){
#                         my %library=(limskey=>$limskey);
#                         map{ $library{$_}=$libraries{$limskey}{$_}} keys %{$libraries{$limskey}};
#                         push(@{$workflow{libraries}},\%library);
#                 }


#                 if($workflows{$wfrunid}{children}){
#                         my @ids=sort keys %{$workflows{$wfrunid}{children}{workflows}};
#                         for my $id(@ids){
#                                 my %child_workflow=(wfrun_id=>$id,wf=>$workflows{$id}{info}{wf},wfv=>$workflows{$id}{info}{wfv});
#                                 push(@{$workflow{children}{workflows}},\%child_workflow);
#                         }
#                 }

#                 if($workflows{$wfrunid}{parents}){
#                         my @ids=sort keys %{$workflows{$wfrunid}{parents}{workflows}};
#                         for my $id(@ids){
#                                 my %parent_workflow=(wfrun_id=>$id,wf=>$workflows{$id}{info}{wf},wfv=>$workflows{$id}{info}{wfv});
#                                 push(@{$workflow{parents}{workflows}},\%parent_workflow);
#                         }
#                         @ids=sort keys %{$workflows{$wfrunid}{parents}{files}};
#                         for my $id(@ids){
#                                 my %parent_file=(fid=>$id,path=>$files{$id}{path});
#                                 push(@{$workflow{parents}{files}},\%parent_file);
#                         }
#                 }

#                 #### pusht the workflow into the data structure
#                 push(@{$data{workflows}},\%workflow);
#         }


#         for my $fid(sort keys %files){
#                 my %file=(fid=>$fid);
#                 map{ $file{$_}=$files{$fid}{$_} } sort keys %{$files{$fid}};

#                 if($$nabu{$fid}){
#                         %{$file{qc}}=%{$$nabu{$fid}};
#                 }

#                 push(@{$data{files}},\%file);
#         }
#         return %data;

# }





# sub validate_options{
#         my ($opts)=@_;

#         usage("Help requested.") if($opts{help});

#         if(! $opts{project}){
#                 usage("ERROR : project must be indicated");
#         }
#         if(! $opts{fpr}){
#                 $opts{fpr}="/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz";
#         }
#         if(! -e $opts{fpr}){
#                 usage("ERROR : file provenance not found at $opts{fpr}");
#         }

#         #$opts{pipelines}="WG";
#         ### set for now to WG
#     #if($opts{pipelines} eq "WG"){
#         #       $opts{workflows}="bcl2fastq,bwaMem,ichorCNA,bamMergePreprocessing,mutect2,variantEffectPredictor,varscan,sequenza,delly,mavis";
#         #}

#         if($opts{out}){
#                 usage("ERROR: output directory $opts{out} does not exist") if(! -d $opts{out});
#         }else{
#                 $opts{out}=".";
#         }
# }

# sub usage{

#         print "\nFPReporter : generation of HTML reports from File Provenance to monitor pipeline output\n";
#         print "\nperl fpreporter.pl [options]\n";
#         print "Options are as follows:\n";
#         print "\t--project.  The project code to summarize [required]\n";
#         print "\t--fpr.  The file provenence report to use, defaults to the system current fpr\n";
#         print "\t--workflows. An ordered, comma separated list of workflows.  Workflows will be shown in the order indicated\n";
#         print "\t--pipelines. A list of standard pipelines to show.  Valid options are WG WT TS ALL\n";
#         print "\t--qc.  Boolean argument to include QC workflows\n";
#         print "\t--release. Boolean argument to show data release information from nabu\n";
#         print "\t--out. Output directory. where to put the html reports\n";
#         print "\t--help displays this usage message.\n";


#         die "\n@_\n\n";
# }



















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
            

# infile = open('TGL01MOH.json')
# data = json.load(infile)
# infile.close()

# S = collect_sequence_info(data)
# L = collect_lims_info(data)
# R = collect_release_status(data)
# D = add_lims_info_to_sequence_data(data)





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





