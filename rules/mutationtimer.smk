rule mutationtimer:
    input:
        vcf=os.path.join(config['results_dir'], '{sample}.formt.vcf'),
        cn=os.path.join(config['results_dir'], '{sample}.formt.cn.tsv')
    output:
        os.path.join(config['results_dir'], '{sample}.pdf')
    log:
        os.path.join(config['log_dir'], '{sample}.pdf.log')
    params:
        clonal_freq = 0.54 # TODO: soft code
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/run_mutationtimer.R '
        '{input.vcf} {input.cn} {params.clonal_freq} {output} '
        '&> {log}'
