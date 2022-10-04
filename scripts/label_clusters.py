import re
import yaml
import argparse
import dask.dataframe as dd
import pandas as pd

INDEX = ['chrom', 'coord', 'ref', 'alt']

def get_args():
    desc = 'Label pseudobulk snvs with SIGNALS clusters/clone_ids'
    p = argparse.ArgumentParser(description=desc)

    p.add_argument('--counts', help="SIGNALS counts file")
    p.add_argument('--inputs', help="SIGNALS input.yaml")
    p.add_argument('--clusters', help="clusters table from SIGNALS Rdata")
    p.add_argument('--variants', help="filtered variant tsv from VARIANTCALLING")

    p.add_argument('-o', '--output', help="output labelled tsv")

    return p.parse_args()

def get_clone_id(clusters):
    # get clusters
    clone_id = pd.read_table(clusters)
    clone_id = clone_id.set_index('cell_id')
    return clone_id

def map_cell_to_clone(inputs, clone_id):
    inputs = yaml.load(open(inputs, 'r').read(), Loader=yaml.Loader)
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
    return cell_to_clone

def map_cluster_to_snvgenotyping(gt_path): # SCDNA-SNVGENOTYPING counts file
    gt = dd.read_csv(gt_path, dtype={'chrom':str, 'cell_id':str, 'coord':'int32',
                                     'ref':str, 'alt':str,
                                     'ref_counts':'int16', 'alt_counts':'int16'})
    gt['cluster'] = gt['cell_id'].map(cell_to_clone)
    gt['cluster'] = gt['cluster'].dropna()

    tumor_gt = gt.loc[gt['cluster'] != '0']
    normal_gt = gt.loc[gt['cluster'] == '0']

    tumor_groups = tumor_gt.groupby(['cluster', 'chrom', 'coord', 'ref', 'alt'])
    tumor_df = tumor_groups[['ref_counts', 'alt_counts']].sum()
    tumor_df = tumor_df.reset_index().compute()
    tumor_df = tumor_df.set_index(INDEX)

    normal_groups = normal_gt.groupby(['cluster', 'chrom', 'coord', 'ref', 'alt'])
    normal_df = normal_groups[['ref_counts', 'alt_counts']].sum()
    normal_df = normal_df.reset_index().compute()
    normal_df = normal_df.set_index(INDEX)

    tumor_df = tumor_df[tumor_df['alt_counts'] > 0]
    normal_df = normal_df[normal_df['alt_counts'] > 0]

    tumor_clust = tumor_df.groupby(INDEX)[['ref_counts', 'alt_counts']].sum()
    tumor_clust['cluster'] = tumor_df.groupby(INDEX)['cluster'].apply(lambda x: ','.join(x))
    normal_clust = normal_df.groupby(INDEX)[['ref_counts', 'alt_counts']].sum()
    normal_clust['cluster'] = normal_df.groupby(INDEX)['cluster'].apply(lambda x: ','.join(x))
    somatic_clust = tumor_clust[~(tumor_clust.index.isin(normal_clust.index))]

    return somatic_clust

def label_variantcalling_with_cluster(somatic_clust, snv, out_tsv_path):
    # join gt to snv, map cluster id to cell id
    snv = snv[INDEX].set_index(INDEX)
    labelled = snv.join(somatic_clust).dropna()

    return labelled

if __name__ == "__main__":
    args = get_args()

    # parse clone id / cluster id
    clone_id = get_clone_id(args.clusters)

    # get cell_to_clone
    cell_to_clone = map_cell_to_clone(args.inputs, clone_id)

    # process SCDNA-SNVGENOTYPING file
    somatic_clust = map_cluster_to_snvgenotyping(args.counts)

    # process filtered snv (from SCDNA-VARIANTCALLING proc pipeline)
    snv_path = args.variants
    snv = dd.read_csv(snv_path, sep='\t', 
            dtype={'chrom':str, 'coord':'int32', 'ref':str, 'alt':str})

    # join gt to snv, map cluster id to cell id
    labelled = label_variantcalling_with_cluster(somatic_clust, snv, args.output)
    labelled.to_csv(out_tsv_path, sep='\t', index=False)
