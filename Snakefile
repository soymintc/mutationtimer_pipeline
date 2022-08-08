
rule all:
    input:
        os.path.join(config['results_dir'], config['sample_id'] + '.pdf')

include: "rules/common.smk"
include: "rules/process_variants.smk"
include: "rules/mutationtimer.smk"
