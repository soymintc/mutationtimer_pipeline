import os
import yaml
import pandas as pd

def get_ploidy(purity_and_ploidy):
    purity, ploidy = open(purity_and_ploidy, 'r').read().strip().split(',')
    return ploidy

def get_purity(purity_and_ploidy):
    purity, ploidy = open(purity_and_ploidy, 'r').read().strip().split(',')
    return purity

rule mutationtimer:
    input:
        vcf=os.path.join(config['results_dir'], '{aliquot_id}.vcf'),
        cn=os.path.join(config['results_dir'], '{sample}.cna_bins_consensus_mutationtimer.tsv'),
        purity_and_ploidy=os.path.join(config['results_dir'], 
                                       '{aliquot_id}.purity,ploidy.txt'),
        clusters=(os.path.join(config['results_dir'], '{aliquot_id}.clusters.tsv')),
    output:
        pdf=os.path.join(config['results_dir'], '{aliquot_id}.pdf'),
        rdata=os.path.join(config['results_dir'], '{aliquot_id}.RData')
    params:
        purity = lambda wildcards, input: get_purity(input.purity_and_ploidy),
        ploidy = lambda wildcards, input: get_ploidy(input.purity_and_ploidy)
    log:
        os.path.join(config['log_dir'], '{aliquot_id}.pdf.log')
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/run_mutationtimer.R '
        '{input.vcf} {input.cn} '
        '{params.purity} ' # clonal_freq in MutationTimeR
        '{params.ploidy} ' # ploidy in MutationTimeR
        '{output.pdf} {output.rdata} '
        '&> {log}'
