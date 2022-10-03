import os
import yaml
import pandas as pd

def get_ploidy(purity_and_ploidy):
    for line in open(purity_and_ploidy, 'r'):
        if line.startswith('purity'): continue
        purity, ploidy = line.rstrip().split(',')
    return ploidy

def get_purity(purity_and_ploidy):
    for line in open(purity_and_ploidy, 'r'):
        if line.startswith('purity'): continue
        purity, ploidy = line.rstrip().split(',')
    return purity

rule mutationtimer:
    input:
        vcf=os.path.join(config['results_dir'], '{sample}.vcf'),
        cn=os.path.join(config['results_dir'], '{sample}.cn.tsv'),
        purity_and_ploidy=os.path.join(config['results_dir'], 
                                       '{sample}.purity_ploidy.csv'),
        clusters=(os.path.join(config['results_dir'], '{sample}.clusters.tsv')),
    output:
        pdf=os.path.join(config['results_dir'], '{sample}.pdf'),
        rdata=os.path.join(config['results_dir'], '{sample}.RData')
    params:
        purity = lambda wildcards, input: get_purity(input.purity_and_ploidy),
        ploidy = lambda wildcards, input: get_ploidy(input.purity_and_ploidy)
    log:
        os.path.join(config['log_dir'], '{sample}.pdf.log')
    resources:
        mem_mb = 4096
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/run_mutationtimer.R '
        '{input.vcf} {input.cn} {input.clusters} '
        '{params.purity} ' # clonal_freq in MutationTimeR
        '{params.ploidy} ' # ploidy in MutationTimeR
        '{output.pdf} {output.rdata} '
        '&> {log}'
