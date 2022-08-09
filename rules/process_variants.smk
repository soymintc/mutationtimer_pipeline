def get_maf_input_paths(wildcards):
    return {
        'consensus_somatic_maf': runinfo.paths[(wildcards.aliquot_id, 'WGS-SOMATICCALLING', 'consensus_somatic_maf')],
        'museq_paired_annotated': runinfo.paths[(wildcards.aliquot_id, 'WGS-SOMATICCALLING', 'museq_paired_annotated')],
        'strelka_snv_annotated': runinfo.paths[(wildcards.aliquot_id, 'WGS-SOMATICCALLING', 'strelka_snv_annotated')],
    }

def get_remixtpp_input_paths(wildcards):
    return {
        'remixt_cn': runinfo.paths[(wildcards.aliquot_id, 'WGS-REMIXT-POSTPROCESS', 'remixt_cn')],
    }

rule filter_maf:
    input: unpack(get_maf_input_paths)
    output: os.path.join(config['intermediate_dir'], '{aliquot_id}_filtered.maf'),
    log: os.path.join(config['log_dir'], '{aliquot_id}_filtered.maf.log'),
    singularity: "docker://soymintc/clickpdvcf:latest"
    shell: 
        'python scripts/filter_maf.py '
        '--in_maf {input.consensus_somatic_maf} --out_maf {output} '
        '&> {log}'

rule maf_to_vcf:
    input:
        maf = os.path.join(config['intermediate_dir'], '{aliquot_id}_filtered.maf'),
        ref = '/juno/work/shah/mondrian/pipelines/mutect/reference/vep/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' # TODO: soft code
    output:
        os.path.join(config['intermediate_dir'], '{aliquot_id}_filtered.vcf')
    log:
        os.path.join(config['log_dir'], '{aliquot_id}_filtered.vcf.log')
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
        os.path.join(config['intermediate_dir'], '{aliquot_id}_filtered.vcf')
    output:
        os.path.join(config['results_dir'], '{aliquot_id}.vcf')
    log:
        os.path.join(config['log_dir'], '{aliquot_id}.vcf.log')
    singularity: 
        "docker://soymintc/clickpdvcf:latest"
    shell:
        'python scripts/proc_consensus_vcf_for_mt.py '
        '--in_vcf {input} --out_vcf {output} &> {log}'

rule proc_remixtpp_for_mt:
    input:
        unpack(get_remixtpp_input_paths)
        #remixtpp = config['remixtpp']
    output:
        os.path.join(config['results_dir'], '{aliquot_id}.cn.tsv')
    log:
        os.path.join(config['log_dir'], '{aliquot_id}.cn.tsv.log')
    params:
        bin_size = int(5e7)
    singularity: 
        "docker://soymintc/clickpdvcf:latest"
    shell:
        'python scripts/proc_remixtpp_for_mt.py '
        '--in_pp {input.remixt_cn} --out_cn {output} --bin_size {params.bin_size} '
        '&> {log}'

