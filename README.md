# MutationTimeR pipeline
A wrapper pipeline for running MutationTimeR with Isabl consensus somatic maf and ReMixT post-processed csv.

## Branches
:art: : feel free to add any branches, or to fork and add branches of your own!

## Install
Cloning the directory is enough. However, you need to provide an accurate Isabl aliquot ID. See example in the next section.
```
git clone git@github.com:soymintc/mutationtimer_pipeline.git
```

## Lazy run example
At least running the following code will result in a successful run with `${out_dir}/results/${aliquot_id}.pdf` as an output.
```
out_dir=/path/to/out/dir  # you need write privilege
aliquot_id=SPECTRUM-OV-081_S1_LEFT_OVARY  # matching Isabl aliquot ID

bash run_snakemake.sh $out_dir $aliquot_id
```

## Input data description for MutationTimeR
You would want to customize however you want to preprocess data for MutationTimerR. Here's a description of how the input should be formatted for this case.<br>
Script `scripts/run_mutationtimer.R` takes four inputs: `vcf`, `cn`, `purity_and_ploidy`, and `clusters` files. `vcf` is advides to end in `.vcf` (not `.vcf.gz`), `cn` in `.cn.tsv`, `purity_and_ploidy` in `.csv`, and `clusters` in `.tsv`.

### `vcf` format
First of all, the vcf input should only include SNVs. It should also have only one genotype column, e.g. `TUMOR` as below. Here's an example visidata block of the `vcf` input:
```bash
 #CHROM │ POS      │ ID │ REF │ ALT │ QUAL │ FILTER │ INFO                          │ FORMAT   │ TUMOR     ║
 1      │ 840748   │ .  │ C   │ T   │ .    │ PASS   │ t_alt_count=2;t_ref_count=12  │ GT:AD:DP │ 0/1:2:14  ║
 1      │ 981866   │ .  │ C   │ A   │ .    │ PASS   │ t_alt_count=8;t_ref_count=35  │ GT:AD:DP │ 0/1:8:43  ║
 1      │ 1066117  │ .  │ G   │ C   │ .    │ PASS   │ t_alt_count=12;t_ref_count=8  │ GT:AD:DP │ 0/1:12:20 ║
```
Columns `#CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `FORMAT`, `TUMOR` are self explanatory. for `INFO`, MutationTimeR requires that you put in `t_alt_count` and `t_ref_count` for the tumor alt and ref counts, respectively. The formats for the `INFO` columns should be described in the header as well, so I advise that you keep the vcf meta-information (i.e. lines starting with `##`) and the header (i.e. line starting with `#CHROM`) as in `resources/meta-information.vcf` and `resources/header.vcf`, respectively.

### `cn` format
This is an arbitrary format that the `run_mutationtimer.R` wrapper for MutationTimeR takes in to process the tsv file as a `GRanges` object. Here's an example visidata block of the `cn` input:
```bash
 chromosome │ start     │ end       │ major_1 │ minor_1 ║
 1          │ 2000001   │ 12500000  │ 2       │ 2       ║
 1          │ 14000001  │ 18500000  │ 2       │ 2       ║
 1          │ 18500001  │ 108000000 │ 2       │ 1       ║
 1          │ 108000001 │ 108500000 │ 2       │ 2       ║
```
All columns are self explanatory except for `major_1` and `minor_1`, which are major and minor copy numbers for the block. It's also noteworthy that the blocks (`start` to `end`) in each line should not overlap, else you will get absurd results. Trust me.

### `purity_and_ploidy` format
This is an arbitrary format that has the following csv table structure (hooray another visidata example!):
```bash
 purity            │ ploidy            ║
 0.973384030418251 │ 3.213179615309869 ║
```
Values are self-explanatory. I'm not skipping explanation due to laziness. Really.<br>
Since MutationTimeR was tested out in WGS samples, you will have to summarize the purity and ploidy for a subclone or pseudobulk when applying it in scDNA results.

### `clusters` format [optional]
This is a table that will be read into a `data.frame`. It describes the information about "known" subclonal mutation clusters. Here's yet another visidata example for the format:
```bash
 cluster │ n_ssms │ proportion          ║
 0       │ 2635   │ 0.686895265171958   ║
 1       │ 603    │ 0.15719083297863023 ║
 2       │ 496    │ 0.1292979322676627  ║
```
The `cluster` column describes the identifiers for the subclones, where ID 0 is treated as the main subclone. `n_ssms` is the number of SNVs assigned to each cluster and is used as a prior downstream. `proportion` is a fraction of the purity for the subclone and does not have to sum up to 1.0, as depicted in https://github.com/gerstung-lab/MutationTimeR/blob/master/README.md#input-data.

## Settings
1. Create an `$out_dir` for running `run_snakemake.sh`.
2. Make sure that the `$aliquot_id` input for `run_snakemake.sh` is exact with the Isabl aliquot ID included in the metadata file.
3. For first-time runs I recommend running `run_snakemake.sh` with `--dry-run`.
4. Let [Seongmin](https://www.github.com/soymintc) know of any inconveniences :yum: and relax :beers:
