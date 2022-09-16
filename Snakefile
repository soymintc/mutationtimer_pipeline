
rule all:
    input:
        os.path.join(config['results_dir'], f'{config["aliquot_id"]}.pdf'),
        os.path.join(config['results_dir'], f'{config["aliquot_id"]}.RData')

include: "rules/common.smk"
include: "rules/process_variants.smk"
include: "rules/mutationtimer.smk"
