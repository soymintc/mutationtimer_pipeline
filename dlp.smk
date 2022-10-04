import os 
import yaml
import pandas as pd

rule all:
    input:
        os.path.join(config['results_dir'], f'{config["sample_id"]}.pdf'),

include: "rules/common_scdna.smk"
include: "rules/process_pseudobulk_variants.smk"
include: "rules/process_pseudobulk_cna.smk"
include: "rules/mutationtimer_scdna.smk"
