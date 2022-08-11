# MutationTimeR pipeline
A wrapper pipeline for running MutationTimeR with Isabl consensus somatic maf and ReMixT post-processed csv.

## Branches
:art: : feel free to add any branches, or to fork and add branches of your own!

## Install
Cloning the directory is enough. However, you need to provide an accurate Isabl aliquot ID. See example in the next section.
```
git clone git@github.com:soymintc/mutationtimer_pipeline.git
git checkout link_isabl  # main is a freeze of the most primitive working version!
```

## Example
At least running the following code will result in a successful run with `${out_dir}/results/${aliquot_id}.pdf` as an output.
```
out_dir=/path/to/out/dir  # you need write privilege
aliquot_id=SPECTRUM-OV-081_S1_LEFT_OVARY  # matching Isabl aliquot ID

bash run_snakemake.sh $out_dir $aliquot_id
```

## Settings
1. Create an `$out_dir` for running `run_snakemake.sh`.
2. Make sure that the `$aliquot_id` input for `run_snakemake.sh` is exact with the Isabl aliquot ID included in the metadata file.
3. For first-time runs I recommend running `run_snakemake.sh` with `--dry-run`.
4. Let [Seongmin](https://www.github.com/soymintc) know of any inconveniences :yum: and relax :beers:
