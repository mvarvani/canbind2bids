configfile: 'config.yml'

from snakebids.utils.snakemake_io import glob_wildcards
from glob import glob
from os.path import join
import os
import json



wildcard_constraints:
    site='[A-Z]{3}',
    subjnum='[0-9]{4}'

#tmp_dir=os.environ['SLURM_TMPDIR']
tmp_dir='/tmp/akhanf'

rule all_unzip:
    input: expand(join(tmp_dir,'{modality}/{group}'),modality=config['modalities'],group=config['groups'])



rule unzip_to_tmp:
    input: config['raw_zip']
    output: directory(join(tmp_dir,'{modality}/{group}'))
    shell: 'unzip -d {output} {input}'


def get_all_t1(wildcards):
    group,study,site,subjnum,sesid1,sesid2 = glob_wildcards(join(tmp_dir,config['in_t1']))
    return expand('bids_{group}/sub-{site}{subjnum}/anat/sub-{site}{subjnum}_T1w.nii.gz',
                zip,group=group,site=site,subjnum=subjnum)

rule all_t1:
    input: get_all_t1
 
      
def get_input_t1(wildcards):
    return glob(join(tmp_dir,config['glob_t1']).format(**wildcards))[-1]
    
rule cp_to_bids_t1:
    input: 
        get_input_t1
    output: 
        'bids_{group}/sub-{site}{subjnum}/anat/sub-{site}{subjnum}_T1w.nii.gz'
    shell:
        'cp {input} {output}'

    
def get_all_dwi(wildcards):
    group,site,subjnum,sesid = glob_wildcards(join(tmp_dir,config['in_dwi']))
    return expand('bids_{group}/sub-{site}{subjnum}/dwi/sub-{site}{subjnum}_dwi.nii.gz',
                zip,group=group,site=site,subjnum=subjnum)

def get_all_dwi_json(wildcards):
    group,site,subjnum,sesid = glob_wildcards(join(tmp_dir,config['in_dwi']))
    return expand('bids_{group}/sub-{site}{subjnum}/dwi/sub-{site}{subjnum}_dwi.json',
                zip,group=group,site=site,subjnum=subjnum)


rule all_dwi:
    input: get_all_dwi
 
rule all_dwi_json:
    input: get_all_dwi_json
      
def get_input_dwi(wildcards):
    dwi_nii=glob(join(tmp_dir,config['glob_dwi']).format(**wildcards))[0]
    dwi_bvec=dwi_nii.replace('.nii.gz','.bvec')
    dwi_bval=dwi_nii.replace('.nii.gz','.bval')
    
    return [dwi_nii,dwi_bvec,dwi_bval]
    
rule cp_to_bids_dwi:
    input: 
        get_input_dwi
    params:
        cp_cmd = lambda wildcards,input,output:  ' && '.join([f'cp {infile} {outfile}' for infile,outfile in zip(input,output)])
    output: 
        nii='bids_{group}/sub-{site}{subjnum}/dwi/sub-{site}{subjnum}_dwi.nii.gz',
        bvec='bids_{group}/sub-{site}{subjnum}/dwi/sub-{site}{subjnum}_dwi.bvec',
        bval='bids_{group}/sub-{site}{subjnum}/dwi/sub-{site}{subjnum}_dwi.bval'
    shell:
        '{params.cp_cmd}'
        
rule cp_to_bids_dwi_json:
    input: 
        'resources/template_dwi.json'
    output: 
        json='bids_{group}/sub-{site}{subjnum}/dwi/sub-{site}{subjnum}_dwi.json',
    shell:
        'cp {input} {output}'
 
    

def get_all_bold(wildcards):
    group,study,site,subjnum,sesid1,sesid2,runid = glob_wildcards(join(tmp_dir,config['in_bold']))
    return expand('bids_{group}/sub-{site}{subjnum}/func/sub-{site}{subjnum}_task-rest_bold.nii.gz',
                zip,group=group,site=site,subjnum=subjnum)

rule all_bold:
    input: get_all_bold
 
      
def get_input_bold(wildcards):
    return glob(join(tmp_dir,config['glob_bold']).format(**wildcards))
    
rule cp_to_bids_bold:
    input: 
        get_input_bold
    output: 
        'bids_{group}/sub-{site}{subjnum}/func/sub-{site}{subjnum}_task-rest_bold.nii.gz'
    shell:
        'cp {input} {output}'

rule all_toplevel:
    input:
        expand('bids_{group}/dataset_description.json',group=config['groups']),
        expand('bids_{group}/task-rest_bold.json',group=config['groups'])

rule all_bids:
    input: 
        rules.all_t1.input,
        rules.all_dwi.input,
        rules.all_dwi_json.input,
        rules.all_bold.input,
        rules.all_toplevel.input

rule dataset_description:
    output:
        json='bids_{group}/dataset_description.json'
    run:
        with open(output.json,'w') as jsonfile:
            dd_dict=config['dataset_description']
            for key in dd_dict.keys():
                dd_dict[key] = dd_dict[key].format(**wildcards)
            jsonfile.write(json.dumps(dd_dict, indent=4))

rule bold_json:
    output:
        json='bids_{group}/task-rest_bold.json'
    run:
        with open(output.json,'w') as jsonfile:
            dd_dict=config['bold_json']
            jsonfile.write(json.dumps(dd_dict, indent=4))


#clinical:

rule all_unzip_clinicalvars:
    input: expand(join(tmp_dir,'clinical/{clinicalvar}_{group}.csv'),clinicalvar=config['clinicalvars'].keys(),group=config['groups'])


rule unzip_clinicalvars_to_tmp:
    input: config['clinicalvar_zip']
    output: directory(join(tmp_dir,'clinical/{clinicalvar}/{group}'))
    shell: 'unzip -d {output} {input}'

rule rename_csv:
    input: 
        join(tmp_dir,'clinical/{clinicalvar}/{group}')
    output:
        join(tmp_dir,'clinical/{clinicalvar}_{group}.csv')
    shell:
        'cp {input}/*/*.csv {output}'

rule split_csv:
    input:
        csv = join(tmp_dir,'clinical/{clinicalvar}_{group}.csv')
    output:
        tsv = 'test/sub-{site}{subjnum}/beh/sub-{site}{subjnum}_task-{task}_beh.tsv'
    script:
        'scripts/split_csv.py'

