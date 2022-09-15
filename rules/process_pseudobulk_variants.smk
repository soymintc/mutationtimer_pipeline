import os
import pandas as pd

def _retrieve_path(paths, isabl_sample_id, result_type):
    data_paths = paths.loc[(isabl_sample_id, result_type)]
    assert data_paths.shape[0] == 1, data_paths
    assert 'result_filepath' in data_paths.columns, data_paths.columns
    path = data_paths['result_filepath'][0]
    assert os.path.exists(path), path
    return path

def _get_scdna_variantcalling(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.SCDNA-VARIANTCALLING.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[(paths["sample_category"] == "TUMOR")]
    paths.set_index(['isabl_sample_id', 'result_type'], inplace=True)

    result_types = ['snv_museq', 'snv_strelka', 'snv_mappability',
        'snv_cosmic_status', 'snv_dbsnp_status', 'snv_snpeff', 'snv_trinuc']
    dlp_paths = {result_type: _retrieve_path(paths, wildcards.sample, result_type)
        for result_type in result_types}
    print(dlp_paths)
    return dlp_paths

rule filter_pseudobulk_snv:
    input: 
        unpack(_get_scdna_variantcalling)
    output: 
        os.path.join(config['results_dir'], '{sample}.pseudobulk.csv')
    params:
        blacklist = '/juno/work/shah/reference/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt'
    log:
        os.path.join(config['log_dir'], '{sample}.pseudobulk.csv.log')
    singularity: 
        "docker://soymintc/clickpdvcf:latest"
    shell:
        'python scripts/filter_pseudobulk_snvs.py '
        '{input.snv_museq} ' 
        '{input.snv_strelka} ' 
        '{input.snv_mappability} ' 
        '{input.snv_cosmic_status} ' 
        '{input.snv_dbsnp_status} ' 
        '{input.snv_snpeff} ' 
        '{input.snv_trinuc} ' 
        '{params.blacklist} '
        '{output} '
        '&> {log}'
