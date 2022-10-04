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

    result_types = ['museq', 'snv_museq', 'snv_strelka', 'snv_mappability',
        'snv_cosmic_status', 'snv_dbsnp_status', 'snv_snpeff', 'snv_trinuc']
    dlp_paths = {result_type: _retrieve_path(paths, wildcards.sample, result_type)
        for result_type in result_types}
    return dlp_paths

rule filter_pseudobulk_snv:
    input: 
        unpack(_get_scdna_variantcalling)
    output: 
        os.path.join(config['results_dir'], '{sample}.pseudobulk_snv.tsv')
    params:
        blacklist = '/juno/work/shah/reference/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt'
    log:
        os.path.join(config['log_dir'], '{sample}.pseudobulk_snv.tsv.log')
    singularity: 
        "docker://soymintc/clickpdvcf:latest"
    shell:
        'python scripts/filter_pseudobulk_snvs.py '
        '{input.museq} ' # museq vcf
        '{input.snv_museq} ' # museq csv 
        '{input.snv_strelka} ' 
        '{input.snv_mappability} ' 
        '{input.snv_cosmic_status} ' 
        '{input.snv_dbsnp_status} ' 
        '{input.snv_snpeff} ' 
        '{input.snv_trinuc} ' 
        '{params.blacklist} '
        '{output} '
        '&> {log}'

rule filtered_pseudobulk_snv_cell_counts:
	input:
        snvs=os.path.join(config['results_dir'], '{sample}.pseudobulk_snv.tsv'),
		counts=lambda wildcards: runinfo.paths[(wildcards.sample, 'SCDNA-SNVGENOTYPING', 'counts')],
		cells="{outdir}/{subdir}/sample_merge/{{sample}}/single_cell/filtered.tsv".format(
			outdir=output_dir,
			subdir=config['outputs']['out']
		),
	output:
		"{outdir}/{subdir}/filter_variants/{{sample}}/pseudobulk/filtered/snv_counts.tsv".format(
			outdir=output_dir,
			subdir=config['outputs']['out']
		),
	log:
		'{outdir}/{subdir}/filtered_pseudobulk_snv_cell_counts/{{sample}}.log'.format(
			outdir=output_dir,
			subdir=config['outputs']['log']
		),
	run:
		import numpy as np

		count_dtypes = {
			'chrom': 'category', 'ref': 'category', 'alt': 'category',
			'cell_id': 'category', 'sample_id': 'category',
			'library_id': 'category'
		}
		counts = pd.read_csv(input['counts'], dtype=count_dtypes)

		retained_cells = np.loadtxt(input['cells'], dtype=str)
		counts = counts[counts['cell_id'].isin(retained_cells)]
		counts = counts[counts['alt_counts'] > 0]
		counts = counts.copy()

		snv_dtypes = {
			'chrom': 'category', 'ref': 'category', 'alt': 'category',
			'tri_nucleotide_context': 'category'
		}
		snvs = pd.read_csv(input['snvs'], sep='\t', dtype=snv_dtypes)

		snvs = snvs.merge(counts)
		snvs.drop_duplicates(inplace=True)
		snvs.to_csv(output[0], sep='\t', index=False)
