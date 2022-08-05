# MutationTimeR pipeline
A wrapper pipeline for running MutationTimeR with Isabl consensus somatic maf and RemiXT post-processed csv
`TODO`: add more docs

## Install
Cloning the directory is enough, but in this version you must provide the `maf` and `remixtpp` variables. See example in the next section.
```
git clone git@github.com:soymintc/mutationtimer_pipeline.git
```

## Example
At least running the following code will result in a successful run with `${out_dir}/results/${sample_id}.pdf` as an output
```
out_dir=./analyses
sample_id=OV-081
maf=/juno/work/shah/isabl_data_lake/analyses/82/06/18206/results/somatic/SPECTRUM-OV-081_S1_LEFT_OVARY/SPECTRUM-OV-081_S1_LEFT_OVARY_consensus_somatic.maf
remixtpp=/juno/work/shah/isabl_data_lake/analyses/37/74/23774/SPECTRUM-OV-081_S1_LEFT_OVARY_cn.csv
bash run_snakemake.sh $out_dir $sample_id $maf $remixtpp
```
