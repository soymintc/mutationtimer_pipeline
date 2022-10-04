rule filter_pseudobulk_snv:
    input: 
        unpack(_get_scdna_variantcalling)
    output: 
        os.path.join(config['results_dir'], '{sample}.filtered_snv.tsv')
    params:
        blacklist = '/juno/work/shah/reference/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt'
    log:
        os.path.join(config['log_dir'], '{sample}.filtered_snv.tsv.log')
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

rule get_clusters: # cluster: clone_id
    input: 
        _get_scdna_signals_rdatafile # rdatafile
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
        '--variants {input.variants} ' # filtered snv
        '--output {output} '
        '&> {log}'

rule count_clusters_and_make_vcf: # make "clusters" input for MutationTimeR
    input:
        labelled = os.path.join(config['results_dir'], '{sample}.labelled_snv.tsv'),
        purityploidy = os.path.join(config['results_dir'], '{sample}.purity_ploidy.csv'),
    output:
        clusters = os.path.join(config['results_dir'], '{sample}.clusters.tsv'),
        vcf = os.path.join(config['results_dir'], '{sample}.vcf'),
    params:
        vcf_header = '/juno/work/shah/users/chois7/retreat/mttest/data/header.vcf',
    log:
        os.path.join(config['log_dir'], '{sample}.count_clusters_and_make_vcf.log'),
    singularity:
        "docker://soymintc/clickpdvcf"
    shell:
        'python scripts/make_clusters_and_vcf.py '
        '--labelled {input.labelled} '
        '--purityploidy {input.purityploidy} '
        '--vcf_header {params.vcf_header} '
        '--out_clusters {output.clusters} '
        '--out_vcf {output.vcf} '
        '&> {log}'
