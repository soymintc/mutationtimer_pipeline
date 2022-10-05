import click
import pandas as pd
from pandas.api.types import CategoricalDtype

def update_blocks_and_reset_prev(blocks, prev, row,
        features=['chromosome', 'start', 'end', 'major_1', 'minor_1']):
    if prev['chromosome']: # not null
        blocks = pd.concat([blocks, prev.to_frame().T])
    prev = pd.Series({f: row[f] for f in features}, index=features)
    return blocks, prev


def get_blocks_from_signals_cna(in_cna):
    #in_cna = '/juno/work/shah/isabl_data_lake/analyses/37/74/23774/SPECTRUM-OV-081_S1_LEFT_OVARY_cn.csv'
    merge_gap = 1000000
    dtype = {'chromosome':str, 'start':int, 'end':int, 'major_1': float, 'minor_1': float}
    cna = pd.read_table(in_cna, dtype=dtype, low_memory=False)
    
    cna = cna.dropna()
    cna['start'] -= 1 # current (start, end) example: (2500001, 3000000) -> (2500000, 3000000)

    features = ['chromosome', 'start', 'end', 'major_1', 'minor_1']
    cna = cna[features] # reduce features for the sake of memory

    major_smaller_than_minor = (cna['major_1'] < cna['minor_1'])
    cna.loc[ major_smaller_than_minor, 'major_1' ] = cna.loc[ major_smaller_than_minor, 'minor_1' ]

    # merge segments with same major_1 and minor_1
    prev = pd.Series({f:None for f in features}, index=features) # init prev cn
    blocks = pd.DataFrame(columns=features) # init empty blocks

    for rix, row in cna.iterrows():
        if prev['chromosome'] != row['chromosome']:
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        if (prev['major_1'] != row['major_1']) or (prev['minor_1'] != row['minor_1']):
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        if abs(prev['end'] - row['start']) > merge_gap:
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        # update prev block
        prev['end'] = row['end']

    blocks = pd.concat([blocks, prev.to_frame().T]) # add prev
    blocks.loc[:, ['start', 'end', 'major_1', 'minor_1']] = blocks.loc[:, ['start', 'end', 'major_1', 'minor_1']].astype(int)
    blocks['start'] = blocks['start'] + 1 # so segments don't overlap; w/o this !mixFlag[j] error occurs in MutationTimeR step
    return blocks


@click.command()
@click.option('--in_cna', required=True,
        help="Input SCDNA-SIGNALS aggregated CNA tsv")
@click.option('--out_cn', required=True,
        help="Output cn tsv path")
def proc_signals_cna_for_mt(in_cna, out_cn):
    cn_data = get_blocks_from_signals_cna(in_cna)
    cn_data.to_csv(out_cn, sep='\t', index=False)


if __name__ == '__main__':
    proc_signals_cna_for_mt()

