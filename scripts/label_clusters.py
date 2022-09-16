import re
import yaml
import argparse
import dask.dataframe as dd
import pandas as pd

def get_args():
    desc = 'Label pseudobulk snvs with SIGNALS clusters/clone_ids'
    p = argparse.ArgumentParser(description=desc)

    p.add_argument('--counts', help="SIGNALS counts file")
    p.add_argument('--inputs', help="SIGNALS input.yaml")
    p.add_argument('--clusters', help="clusters table from SIGNALS Rdata")
    p.add_argument('--variants', help="filtered variant tsv from VARIANTCALLING")

    p.add_argument('-o', '--output', help="output labelled tsv")

    return p.parse_args()

args = get_args()

# get clusters
clone_id = pd.read_csv(args.clusters)
clone_id = clone_id.set_index('cell_id')

# get cell_to_clone
inputs = yaml.load(open(args.inputs, 'r').read(), Loader=yaml.Loader)
samples = inputs['tumour_cells'].keys()
assert len(samples) == 1
sample = list(samples)[0]
keys = inputs['tumour_cells'][sample].keys()
assert len(keys) == 1
key = list(keys)[0]
cell_ids = list(inputs['tumour_cells'][sample][key].keys())
normal_cell_ids = [i for i in cell_ids if not i.startswith(sample)]
normal_cell_ids = sorted(normal_cell_ids, key = lambda x: re.search('R(\d+)-C(\d+)', x).groups())
tumor_cell_ids = [i for i in cell_ids if i.startswith(sample)]
tumor_cell_ids = sorted(tumor_cell_ids, key = lambda x: re.search('R(\d+)-C(\d+)', x).groups())

cell_to_clone = clone_id.to_dict()['clone_id']
for normal_cell_id in normal_cell_ids:
    cell_to_clone[normal_cell_id] = '0'

# process SNV-GENOTYPING file
varid = ['chrom', 'coord', 'ref', 'alt']

gt_path = '/juno/work/shah/users/chois7/retreat/mttest/test.csv.gz' #args.counts
gt = dd.read_csv(gt_path, dtype={'chrom':str, 'cell_id':str, 'coord':'int32',
                                 'ref':str, 'alt':str,
                                 'ref_counts':'int16', 'alt_counts':'int16'})
gt = gt[varid + ['ref_counts', 'alt_counts', 'cell_id']]

# process snv
snv_path = args.variants
snv = dd.read_csv(snv_path, sep='\t', dtype={'chrom':str, 'coord':'int32', 'ref':str, 'alt':str})

# get merged snv data
joint = snv.merge(gt, left_on=varid, right_on=varid, how="inner")
joint = joint.compute() # dd -> pd
joint['cluster'] = joint['cell_id'].map(cell_to_clone).fillna('NA')

cell_and_cluster = joint[varid + ['cell_id', 'cluster']]
cell_ids = cell_and_cluster.groupby(varid)['cell_id'].apply(lambda x: ','.join(x))
clusters = cell_and_cluster.groupby(varid)['cluster'].apply(lambda x: ','.join(x))

labelled = (cell_and_cluster[varid]
        .drop_duplicates()
        .merge(cell_ids, on=varid)
        .merge(clusters, on=varid))

labelled.to_csv(args.output, sep='\t', index=False)
