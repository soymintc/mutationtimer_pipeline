
rule all:
    input:
        os.path.join(config['results_dir'], f'{config["sample_id"]}.purity_ploidy.csv'),
        os.path.join(config['results_dir'], f'{config["sample_id"]}.pseudobulk_snv.tsv'),
        os.path.join(config['results_dir'], f'{config["sample_id"]}.labelled_snv.tsv'),
        os.path.join(config['results_dir'], f'{config["sample_id"]}.pdf'),

include: "rules/process_pseudobulk_variants.smk"
include: "rules/label_clusters.smk"
include: "rules/process_single_cell.smk"
include: "rules/mutationtimer_scdna.smk"
