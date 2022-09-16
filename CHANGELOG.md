# Version history

## 0.0.1
- freeze of working version, but needs `maf` and `remixtpp` (remixt post-processed cn file) to be given as input

## 0.0.2
- changed `sample_id` to `aliquot_id`
- `run_mutationtimer.R` creates an RData file as well
- purity and ploidy parsing from Remixt Post-Processing google spreadsheet
- updated dockerfile to include scipy
- calculate `isWgd` parameter for MutationTimeR

## 0.0.3
- added a more stringent germline filter of `n_alt_count` == 0
