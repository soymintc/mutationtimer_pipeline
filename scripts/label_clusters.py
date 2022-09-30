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

def proc_snvgenotyping(gt_path): # SCDNA-SNVGENOTYPING counts file
    gt = dd.read_csv(gt_path, dtype={'chrom':str, 'cell_id':str, 'coord':'int32',
                                     'ref':str, 'alt':str,
                                     'ref_counts':'int16', 'alt_counts':'int16'})
    gt = gt[INDEX + ['ref_counts', 'alt_counts', 'cell_id']]
    gt = gt.loc[gt['alt_counts'] > 0, :]
    return gt

def map_cluster_to_variants(gt, snv, out_tsv_path):
    # join gt to snv, map cluster id to cell id
    joint = snv.merge(gt, left_on=INDEX, right_on=INDEX, how="inner")
    joint = joint.compute() # dask -> pandas
    joint['cluster'] = joint['cell_id'].map(cell_to_clone).fillna('NA')

    clust_info = joint[INDEX + ['cell_id', 'cluster']]
    cell_ids = clust_info.groupby(INDEX)['cell_id'].apply(lambda x: ','.join(x))
    clusters = clust_info.groupby(INDEX)['cluster'].apply(lambda x: ','.join(x))

    labelled = (clust_info[INDEX]
            .drop_duplicates()
            .merge(cell_ids, on=INDEX)
            .merge(clusters, on=INDEX))

    labelled.to_csv(out_tsv_path, sep='\t', index=False)
    return labelled

if __name__ == "__main__":
    args = get_args()

    # parse clone id / cluster id
    clone_id = get_clone_id(args.clusters)

    # get cell_to_clone
    cell_to_clone = map_cell_to_clone(args.inputs, clone_id)

    # process SCDNA-SNVGENOTYPING file
    gt = proc_snvgenotyping(args.counts)

    # process filtered snv (from SCDNA-VARIANTCALLING proc pipeline)
    snv_path = args.variants
    snv = dd.read_csv(snv_path, sep='\t', 
            dtype={'chrom':str, 'coord':'int32', 'ref':str, 'alt':str})

    # join gt to snv, map cluster id to cell id
    labelled = map_cluster_to_variants(gt, snv, args.output)
