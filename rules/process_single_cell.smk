import os
import pandas as pd

def _retrieve_path(paths, isabl_sample_id, result_type):
    data_paths = paths.loc[[(isabl_sample_id, result_type)]]
    assert data_paths.shape[0] == 1, data_paths
    assert 'result_filepath' in data_paths.columns, data_paths.columns
    path = data_paths['result_filepath'][0]
    assert os.path.exists(path), path
    return path

def _get_scdna_signals(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.SCDNA-SIGNALS.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[(paths["sample_category"] == "TUMOR")]
    paths.set_index(['isabl_sample_id', 'result_type'], inplace=True)

    result_types = ['rdatafile']
    dlp_paths = {result_type: _retrieve_path(paths, wildcards.sample, result_type)
        for result_type in result_types}
    return dlp_paths

rule process_single_cell_copy_number:
    input: 
        unpack(_get_scdna_signals)
    output: 
        cna_bins_consensus=os.path.join(config['results_dir'], '{sample}.cna_bins_consensus.tsv'),
        cna_bins_consensus_mutationtimer=os.path.join(config['results_dir'], '{sample}.cna_bins_consensus_mutationtimer.tsv'),
        purity_ploidy=os.path.join(config['results_dir'], '{sample}.purity_ploidy.csv')
    log:
        os.path.join(config['log_dir'], '{sample}.single_cell_copy_number.tsv.log')
    singularity: 
        "/juno/work/shah/vazquezi/images/singularity/spectrum_latest.sif"
    shell:
        'Rscript scripts/proc_copy_number.R '
        '--hscn {input.rdatafile} ' # Rdata file
        '--cna_bins_consensus {output.cna_bins_consensus} '
        '--cna_bins_consensus_mutationtimer {output.cna_bins_consensus_mutationtimer} '
        '--purity_ploidy {output.purity_ploidy} '
        '&> {log}'

rule proc_signals_cna_for_mt:
    input:
        os.path.join(config['results_dir'], '{sample}.cna_bins_consensus_mutationtimer.tsv'),
    output:
        os.path.join(config['results_dir'], '{sample}.cn.tsv'),
    log:
        os.path.join(config['log_dir'], '{sample}.cn.tsv.log'),
    singularity: 
        "docker://soymintc/clickpdvcf:0.0.2"
    shell:
        'python scripts/proc_signals_cna_for_mt.py '
        '--in_cna {input} --out_cn {output} '
        '&> {log}'

