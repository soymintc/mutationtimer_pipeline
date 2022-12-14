import click
import pandas as pd
from pandas.api.types import CategoricalDtype
import scipy.ndimage

def update_blocks_and_reset_prev(blocks, prev, row, 
        features=['chromosome', 'start', 'end', 'major_1', 'minor_1']):
    if prev['chromosome']: # not null
        blocks = pd.concat([blocks, prev.to_frame().T])
    prev = pd.Series({f: row[f] for f in features}, index=features)
    return blocks, prev


def get_blocks_from_remixt_pp(in_pp):
    #in_pp = '/juno/work/shah/isabl_data_lake/analyses/37/74/23774/SPECTRUM-OV-081_S1_LEFT_OVARY_cn.csv'
    remixt = pd.read_table(in_pp, dtype=str, low_memory=False)
    gaussian_filter_sigma = 30

    features = ['chromosome', 'start', 'end', 'major_1', 'minor_1']
    remixt = remixt[features] # reduce features for the sake of memory

    # smoothen major_1 and minor_1
    for f in ['major_1', 'minor_1']:
        remixt[f] = (scipy.ndimage.gaussian_filter(
                    remixt[f].astype(float), gaussian_filter_sigma
                )
                .round()
                .astype(int)
            )
    major_smaller_than_minor = (remixt['major_1'] < remixt['minor_1'])
    remixt.loc[ major_smaller_than_minor, 'major_1' ] = remixt.loc[ major_smaller_than_minor, 'minor_1' ]

    # merge segments with same major_1 and minor_1
    prev = pd.Series({f:None for f in features}, index=features) # init prev cn
    blocks = pd.DataFrame(columns=features) # init empty blocks

    for rix, row in remixt.iterrows():
        if prev['chromosome'] != row['chromosome']:
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        if (prev['major_1'] != row['major_1']) or (prev['minor_1'] != row['minor_1']):
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        if prev['end'] != row['start']:
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        # update prev block
        prev['end'] = row['end']

    blocks = pd.concat([blocks, prev.to_frame().T]) # add prev
    blocks.loc[:, ['start', 'end', 'major_1', 'minor_1']] = blocks.loc[:, ['start', 'end', 'major_1', 'minor_1']].astype(int)
    blocks['start'] = blocks['start'] + 1 # so segments don't overlap; w/o this !mixFlag[j] error occurs in MutationTimeR step
    return blocks


@click.command()
@click.option('--in_pp', required=True,
        help="Input RemiXT PostProcessing tsv")
@click.option('--out_cn', required=True,
        help="Output cn tsv path")
def proc_remixtpp_for_mt(in_pp, out_cn):
    cn_data = get_blocks_from_remixt_pp(in_pp)
    cn_data.to_csv(out_cn, sep='\t', index=False)


if __name__ == '__main__':
    proc_remixtpp_for_mt()
