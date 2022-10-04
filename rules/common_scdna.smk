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
    variants = os.path.join(config['results_dir'], '{sample}.filtered_snv.tsv')
    return {
        'counts': counts,
        'inputs': inputs,
        'clusters': clusters,
        'variants': variants,
    }

def _get_scdna_signals_rdatafile(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.SCDNA-SIGNALS.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[(paths["sample_category"] == "TUMOR")]
    paths.set_index(['isabl_sample_id', 'result_type'], inplace=True)
    path = _retrieve_path(paths, wildcards.sample, 'rdatafile')
    return path

def _get_scdna_variantcalling(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.SCDNA-VARIANTCALLING.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[(paths["sample_category"] == "TUMOR")]
    paths.set_index(['isabl_sample_id', 'result_type'], inplace=True)

    result_types = ['museq', 'snv_museq', 'snv_strelka', 'snv_mappability',
        'snv_cosmic_status', 'snv_dbsnp_status', 'snv_snpeff', 'snv_trinuc']
    dlp_paths = {result_type: _retrieve_path(paths, wildcards.sample, result_type)
        for result_type in result_types}
    return dlp_paths
