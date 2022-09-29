import os
import yaml
import pandas as pd

def _retrieve_path(paths, isabl_sample_id, result_type):
    data_paths = paths.loc[[(isabl_sample_id, result_type)]]
    assert data_paths.shape[0] == 1, f'data_paths.shape = {data_paths.shape}, data_paths=\n{data_paths.values}'
    assert 'result_filepath' in data_paths.columns, data_paths.columns
    path = data_paths['result_filepath'][0]
    assert os.path.exists(path), path
    return path

def _get_scdna_genotyping_data(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.SCDNA-SNVGENOTYPING.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[(paths["sample_category"] == "TUMOR")]
    paths.set_index(['isabl_sample_id', 'result_type'], inplace=True)
    counts = _retrieve_path(paths, wildcards.sample, 'counts')
    inputs = os.path.join(os.path.split(counts)[0], 'input.yaml')
    clusters = os.path.join(config['results_dir'], '{sample}.signals.clone_id.tsv')
    variants = os.path.join(config['results_dir'], '{sample}.pseudobulk_snv.tsv')
    return {
        'counts': counts,
        'inputs': inputs,
        'clusters': clusters,
        'variants': variants,
    }

def _get_scdna_signals(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.SCDNA-SIGNALS.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[(paths["sample_category"] == "TUMOR")]
    paths.set_index(['isabl_sample_id', 'result_type'], inplace=True)
    path = _retrieve_path(paths, wildcards.sample, 'rdatafile')
    return path

rule get_clusters: # cluster: clone_id
    input: 
        _get_scdna_signals # rdatafile
    output: 
        os.path.join(config['results_dir'], '{sample}.signals.clone_id.tsv')
    log:
        os.path.join(config['log_dir'], '{sample}.signals.clone_id.tsv.log')
    singularity:
        "docker://rocker/tidyverse"
    shell:
        'Rscript scripts/get_clone_id.R '
        '{input} {output} &> {log}'
        
rule label_clusters:
    input:
        unpack(_get_scdna_genotyping_data)
    output:
        os.path.join(config['results_dir'], '{sample}.labelled_snv.tsv')
    log:
        os.path.join(config['log_dir'], '{sample}.labelled_snv.tsv.log')
    singularity:
        "docker://soymintc/dask"
    shell:
        'python scripts/label_clusters.py '
        '--counts {input.counts} '
        '--inputs {input.inputs} '
        '--clusters {input.clusters} '
        '--variants {input.variants} '
        '--output {output} '
        '&> {log}'
