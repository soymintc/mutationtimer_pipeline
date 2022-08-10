import pandas as pd

rule parse_purity:
    input:
        vcf=os.path.join(config['results_dir'], '{aliquot_id}.vcf') # not used
    output:
        os.path.join(config['results_dir'], '{aliquot_id}.purity.txt')
    run:
        sheet_id = "1yZ0UYDm5JuY1FqBLImHjdHeBK6_6EfGpGwyb2X030NM" # TODO: include in config for info safety
        sheet_name = "main"
        url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/edit#gid=0&sheet={sheet_name}"
        url = url.replace('/edit#gid=', '/export?format=csv&gid=')
        df = pd.read_csv(url, header=0)
        df = df[df['tumor_aliquot_id'].notna()]
        spectrum = df[ df['tumor_aliquot_id'].str.startswith('SPECTRUM') ]
        aliquot_and_purity = spectrum[['tumor_aliquot_id', 'tumour_proportion']].fillna('NA')
        aliquot_and_purity = aliquot_and_purity.sort_values(
            by=['tumour_proportion', 'tumor_aliquot_id'], 
            ascending=False) # if duplicate aliquots with both tumour_proportion and NA, use tumour_proportion values instead of NA
        aliquot_to_purity = {x[0]: x[1] for x in aliquot_and_purity.values}
        purity = 'NA'
        if wildcards.aliquot_id in aliquot_to_purity:
            purity = aliquot_to_purity[wildcards.aliquot_id]
        purity = 0.5 if purity == 'NA' else purity # ad hoc replacement to 0.5 if NA
        with open(output[0], 'w') as outfile:
            outfile.write(f'{purity}\n')

def get_tumor_purity(purity_file):
    purity = open(purity_file, 'r').read().strip()
    return purity

rule mutationtimer:
    input:
        vcf=os.path.join(config['results_dir'], '{aliquot_id}.vcf'),
        cn=os.path.join(config['results_dir'], '{aliquot_id}.cn.tsv'),
        purity=os.path.join(config['results_dir'], '{aliquot_id}.purity.txt')
    output:
        pdf=os.path.join(config['results_dir'], '{aliquot_id}.pdf'),
        rdata=os.path.join(config['results_dir'], '{aliquot_id}.RData')
    params:
        clonal_freq = lambda wildcards, input: get_tumor_purity(input.purity)
    log:
        os.path.join(config['log_dir'], '{aliquot_id}.pdf.log')
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/run_mutationtimer.R '
        '{input.vcf} {input.cn} {params.clonal_freq} '
        '{output.pdf} {output.rdata} '
        '&> {log}'
