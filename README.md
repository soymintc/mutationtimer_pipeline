# MutationTimeR pipeline
A wrapper pipeline for running MutationTimeR with Isabl consensus somatic maf and ReMixT post-processed csv.
Welcome to the `link_isabl` branch; adding isabl metadata to the pipeline. :construction_worker:

## Install
Cloning the directory is enough. However, you need to provide an accurate Isabl sample ID. See example in the next section.
```
git clone git@github.com:soymintc/mutationtimer_pipeline.git
```

## Example
At least running the following code will result in a successful run with `${out_dir}/results/${sample_id}.pdf` as an output.
```
out_dir=./analysis  # you need write privilege
sample_id=SPECTRUM-OV-081_S1_LEFT_OVARY  # matching Isabl sample ID

bash run_snakemake.sh $out_dir $sample_id
```

## Settings
1. Create an `$out_dir` for running `run_snakemake.sh`.
2. Make sure that the `$sample_id` input for `run_snakemake.sh` is exact with the Isabl sample ID included in the metadata file.
3. For first-time runs I recommend running `run_snakemake.sh` with `--dry-run`.
4. Let [Seongmin] (https://www.github.com/soymintc) know of any inconveniences :yum: and relax :beers:
