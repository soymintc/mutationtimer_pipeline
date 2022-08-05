import click
import pandas as pd
from pandas.api.types import CategoricalDtype

def update_blocks_and_reset_prev(blocks, prev, row, 
        features=['chromosome', 'start', 'end', 'major_1', 'minor_1']):
    if prev['chromosome']: # not null
        blocks = pd.concat([blocks, prev.to_frame().T])
    prev = pd.Series({f: row[f] for f in features}, index=features)
    return blocks, prev



def get_blocks_from_remixt_pp(in_pp):
    #in_pp = '/juno/work/shah/isabl_data_lake/analyses/37/74/23774/SPECTRUM-OV-081_S1_LEFT_OVARY_cn.csv'
    pp = pd.read_table(in_pp, dtype=str, low_memory=False)

    features = ['chromosome', 'start', 'end', 'major_1', 'minor_1']
    parsed_pp = pp[features]

    prev = pd.Series({f:None for f in features}, index=features) # init prev cn
    blocks = pd.DataFrame(columns=features) # init empty blocks

    for rix, row in parsed_pp.iterrows():
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
    return blocks


def get_cn_data_from_blocks(blocks, bin_size=int(5e7)):
    blocks['start_bin'] = blocks['start'] // bin_size
    blocks['end_bin'] = blocks['end'] // bin_size
    major_groups = blocks.groupby(['chromosome', 'start_bin'])['major_1'].mean().astype(int)
    minor_groups = blocks.groupby(['chromosome', 'start_bin'])['minor_1'].mean().astype(int)

    chroms, start_poss, end_poss, majors, minors = [], [], [], [], []
    for (chrom, pix), major in major_groups.iteritems():
        minor = minor_groups[(chrom, pix)]
        start_pos = bin_size * pix + 1
        end_pos = bin_size * (pix + 1)
        chroms.append(chrom)
        start_poss.append(start_pos)
        end_poss.append(end_pos)
        majors.append(major)
        minors.append(minor)
    cn_data = pd.DataFrame({
        'chromosome': chroms,
        'start': start_poss,
        'end': end_poss,
        'major_1': majors,
        'minor_1': minors,
    })

    chrom_list = CategoricalDtype(
        [str(c) for c in range(1, 22+1)] + ['X', 'Y'], 
        ordered=True
    )
    cn_data['chromosome'] = cn_data['chromosome'].astype(chrom_list)
    cn_data = cn_data.sort_values(by=['chromosome', 'start'])
    return cn_data


@click.command()
@click.option('--in_pp', required=True,
        help="Input RemiXT PostProcessing tsv")
@click.option('--out_cn', required=True,
        help="Output cn tsv path")
@click.option('--bin_size', default=int(5e7), type=int,
        show_default=True, help="CN bin size")
def proc_remixtpp_for_mt(in_pp, out_cn, bin_size):
    blocks = get_blocks_from_remixt_pp(in_pp)
    cn_data = get_cn_data_from_blocks(blocks, bin_size=bin_size)
    cn_data.to_csv(out_cn, sep='\t', index=False)


if __name__ == '__main__':
    proc_remixtpp_for_mt()
