import os
import yaml
import pandas as pd

def extract_purity(aliquot_id, annot):
    assert 'tumor_aliquot_id' in annot.columns
    assert 'tumour_proportion' in annot.columns # purity | clonal frequency
    aliquot_and_purity = annot[['tumor_aliquot_id', 'tumour_proportion']].fillna('NA')
    aliquot_and_purity = aliquot_and_purity.sort_values(
        by=['tumour_proportion', 'tumor_aliquot_id'],
        ascending=False) # if duplicate aliquots with both tumour_proportion and NA, use tumour_proportion values instead of NA
    aliquot_to_purity = {x[0]: x[1] for x in aliquot_and_purity.values}
    purity = 'NA'
    if aliquot_id in aliquot_to_purity:
        purity = aliquot_to_purity[aliquot_id]
    purity = 0.5 if purity == 'NA' else purity # ad hoc replacement to 0.5 if NA
    return purity

def extract_ploidy(aliquot_id):
    remixtpp_path = '/juno/work/shah/users/chois7/tables/all.WGS-REMIXT-POSTPROCESS.tsv'
    assert os.path.exists(remixtpp_path)
    remixtpp = pd.read_table(remixtpp_path)
    meta = remixtpp[remixtpp['result_type']=='meta']
    meta = meta[meta['isabl_aliquot_id'] == aliquot_id]['result_filepath'].values
    assert len(meta) == 1
    assert os.path.exists(meta[0])
    meta = yaml.load(open(meta[0], 'r').read(), Loader=yaml.Loader)
    ploidy = meta['ploidy']
    return ploidy

rule parse_purity_and_ploidy:
    input:
        vcf=os.path.join(config['results_dir'], '{aliquot_id}.vcf') # not used
    output:
        os.path.join(config['results_dir'], '{aliquot_id}.purity,ploidy.txt')
    run:
        sheet_id = "1yZ0UYDm5JuY1FqBLImHjdHeBK6_6EfGpGwyb2X030NM" # TODO: include in config for info safety
        sheet_name = "main"
        url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/edit#gid=0&sheet={sheet_name}"
        url = url.replace('/edit#gid=', '/export?format=csv&gid=')
        df = pd.read_csv(url, header=0)
        df = df[df['tumor_aliquot_id'].notna()]

        annot = df[ df['tumor_aliquot_id'].str.startswith('SPECTRUM') ] # TODO: not only SPECTRUM...

        # get purity
        print("No problem before extract_purity")
        purity = extract_purity(wildcards.aliquot_id, annot)

        # get ploidy
        print("No problem before extract_ploidy")
        ploidy = extract_ploidy(wildcards.aliquot_id)

        with open(output[0], 'w') as outfile:
            outfile.write(f'{purity},{ploidy}\n')

def get_ploidy(purity_and_ploidy):
    purity, ploidy = open(purity_and_ploidy, 'r').read().strip().split(',')
    return ploidy

def get_purity(purity_and_ploidy):
    purity, ploidy = open(purity_and_ploidy, 'r').read().strip().split(',')
    return purity

rule mutationtimer:
    input:
        vcf=os.path.join(config['results_dir'], '{aliquot_id}.vcf'),
        cn=os.path.join(config['results_dir'], '{aliquot_id}.cn.tsv'),
        purity_and_ploidy=os.path.join(config['results_dir'], 
                                       '{aliquot_id}.purity,ploidy.txt')
    output:
        pdf=os.path.join(config['results_dir'], '{aliquot_id}.pdf'),
        rdata=os.path.join(config['results_dir'], '{aliquot_id}.RData')
    params:
        purity = lambda wildcards, input: get_purity(input.purity_and_ploidy),
        ploidy = lambda wildcards, input: get_ploidy(input.purity_and_ploidy)
    log:
        os.path.join(config['log_dir'], '{aliquot_id}.pdf.log')
    resources:
        mem_mb = 4096
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/run_mutationtimer.R '
        '{input.vcf} {input.cn} '
        '{params.purity} ' # clonal_freq in MutationTimeR
        '{params.ploidy} ' # ploidy in MutationTimeR
        '{output.pdf} {output.rdata} '
        '&> {log}'
