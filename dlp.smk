
rule all:
    input:
        os.path.join(config['results_dir'], f'{config["sample_id"]}.pseudobulk.tsv'),
        os.path.join(config['results_dir'], f'{config["sample_id"]}.labelled.tsv'),
        #os.path.join(config['results_dir'], f'{config["sample_id"]}.pdf'),

include: "rules/process_pseudobulk_variants.smk"
include: "rules/label_clusters.smk"
#include: "rules/mutationtimer.smk"
