
rule all:
    input:
        os.path.join(config['results_dir'], f'{config["sample_id"]}.pseudobulk.csv'),
        #os.path.join(config['results_dir'], f'{config["sample_id"]}.pdf'),
include: "rules/process_pseudobulk_variants.smk"
#include: "rules/mutationtimer.smk"
