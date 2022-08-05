import os 

rule all:
    input:
        os.path.join(config['results_dir'], config['sample_id'] + '.formt.vcf'),
        os.path.join(config['results_dir'], config['sample_id'] + '.formt.cn.tsv')

rule maf_to_vcf:
    input:
        maf = config['maf'],
        ref = '/juno/work/shah/mondrian/pipelines/mutect/reference/vep/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' # TODO: soft code
    output:
        os.path.join(config['results_dir'], config['sample_id'] + '.vcf')
    log:
        os.path.join(config['log_dir'], config['sample_id'] + '.vcf.log')
    params:
        outdir = config['results_dir']
    singularity: '/home/chois7/chois7/singularity/sif/var.sif' # TODO: soft code
    shell: 
        'maf2vcf.pl --input-maf {input.maf} '
        '--ref-fasta {input.ref} '
        '--output-dir {params.outdir} '
        '--output-vcf {output} '
        '&> {log}'

rule proc_vcf_for_mt:
    input:
        os.path.join(config['results_dir'], config['sample_id'] + '.vcf')
    output:
        os.path.join(config['results_dir'], config['sample_id'] + '.formt.vcf')
    log:
        os.path.join(config['log_dir'], config['sample_id'] + '.formt.vcf.log')
    shell:
        'python proc_consensus_vcf_for_mt.py '
        '--in_vcf {input} --out_vcf {output} &> {log}'

rule proc_remixtpp_for_mt:
    input:
        remixtpp = config['remixtpp']
    output:
        os.path.join(config['results_dir'], config['sample_id'] + '.formt.cn.tsv')
    log:
        os.path.join(config['log_dir'], config['sample_id'] + '.formt.cn.tsv.log')
    params:
        bin_size = int(5e7)
    shell:
        'python proc_remixtpp_for_mt.py '
        '--in_pp {input} --out_cn {output} --bin_size {params.bin_size} '
        '&> {log}'

