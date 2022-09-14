from argparse import ArgumentParser
import numpy as np
import pandas as pd


def get_args():
    desc = 'Extract and filter pseudobulk SNV counts'
    p = ArgumentParser(description=desc)

    p.add_argument('museq', help='museq csv file')
    p.add_argument('strelka', help='strelka snv csv file')
    p.add_argument('mappability', help='mappability csv file')
    p.add_argument('cosmic', help='cosmic annotation csv file')
    p.add_argument('dbsnp', help='dbsnp annotation csv file')
    p.add_argument('snpeff', help='snpeff annotation csv file')
    p.add_argument('trinuc', help='trinucleotide context annotation csv file')
    p.add_argument('blacklist', help='blacklist region bed file')

    p.add_argument(
        '--remove-chr6p', action='store_true', help='remove all SNVs in chr6 p-arm'
    )

    p.add_argument('snvs', help='output snvs tsv file')

    return p.parse_args()


def in_any_region(pos, blacklist):
    idx = np.searchsorted(blacklist['start'], pos, side='right') - 1
    in_region = pos <= blacklist.loc[blacklist.index[idx], 'end'].values
    in_region &= idx >= 0
    return in_region


def filter_blacklist(snvs, blacklist):
    ok_idx = []
    for chrom in snvs['chrom'].unique():
        chrom_snvs = snvs[snvs['chrom'] == chrom]
        chrom_blacklist = blacklist[blacklist['chr'] == ('chr' + chrom)]
        try:
            is_low_mapp = in_any_region(chrom_snvs['coord'], chrom_blacklist)
        except IndexError:
            print(chrom_snvs['coord'])
            print(chrom_snvs['coord'].shape)
            print(chrom_blacklist)
            print(chrom_blacklist.shape)
        ok_idx.extend(chrom_snvs.index[~is_low_mapp])

    return snvs.loc[ok_idx]


# baed on scgenome
def filter_snvs(museq, strelka, mappability, cosmic, dbsnp, snpeff, trinuc,
        blacklist, remove_chr6p):

    museq = museq.rename(columns={'score': 'museq_score'})
    museq = museq[museq['museq_score'] > 0.9]

    print(len(museq))
    if remove_chr6p:
        is_chr6p = (museq['chrom'] == '6') & (museq['coord'] < 65e5)
        museq = museq[~is_chr6p]
    print(len(museq))

    strelka = strelka.rename(columns={'score': 'strelka_score'})
    strelka = strelka[strelka['strelka_score'] > 20.]

    mappability = mappability[mappability['mappability'] > 0.99]

    cosmic['is_cosmic'] = True
    cosmic = cosmic[[
        'chrom', 'coord', 'ref', 'alt', 'is_cosmic'
    ]].drop_duplicates()

    dbsnp['is_dbsnp'] = True
    dbsnp = dbsnp[[
        'chrom', 'coord', 'ref', 'alt', 'is_dbsnp'
    ]].drop_duplicates()

    snpeff = snpeff[[
        'chrom', 'coord', 'ref', 'alt', 'effect_impact'
    ]].drop_duplicates()
    snpeff['value'] = 1
    snpeff = snpeff.set_index([
        'chrom', 'coord', 'ref', 'alt', 'effect_impact'
    ])['value'].unstack(fill_value=0)
    snpeff = snpeff.rename(columns=str).reset_index()

    snvs = museq.merge(strelka)
    snvs = snvs.merge(mappability)
    snvs = snvs.merge(cosmic, how='left')
    snvs = snvs.merge(dbsnp, how='left')
    snvs = snvs.merge(snpeff, how='left')
    snvs = snvs.merge(trinuc, how='left')

    snvs['is_cosmic'] = snvs['is_cosmic'].fillna(False)
    snvs['is_dbsnp'] = snvs['is_dbsnp'].fillna(False)

    snvs = snvs[~snvs['is_dbsnp']]

    snvs = filter_blacklist(snvs, blacklist)
    snvs.drop_duplicates(inplace=True)

    return snvs

if __name__ == '__main__':
    argv = get_args()

    museq = pd.read_csv(argv.museq, dtype={'chrom': str})
    strelka = pd.read_csv(argv.strelka, dtype={'chrom': str})

    if len(museq) == 0 or len(strelka) == 0:
        snvs = pd.DataFrame(
            columns=[
                'chrom',
                'coord',
                'ref',
                'alt',
                'museq_score',
                'strelka_score',
                'mappability',
                'is_cosmic',
                'is_dbsnp',
                'HIGH',
                'LOW',
                'MODERATE',
                'MODIFIER',
                'tri_nucleotide_context'
            ],
            index=[]
        )
    else:
        mappability = pd.read_csv(argv.mappability, dtype={'chrom': str})
        cosmic = pd.read_csv(argv.cosmic, dtype={'chrom': str})
        dbsnp = pd.read_csv(argv.dbsnp, dtype={'chrom': str})
        snpeff = pd.read_csv(argv.snpeff, dtype={'chrom': str})
        trinuc = pd.read_csv(argv.trinuc, dtype={'chrom': str})
        blacklist = pd.read_csv(
            argv.blacklist, sep='\t', header=None,
            names=['chr', 'start', 'end', 'desc']
        )

        snvs = filter_snvs(
            museq, strelka, mappability, cosmic, dbsnp, snpeff, trinuc,
            blacklist, argv.remove_chr6p
        )

    snvs.to_csv(argv.snvs, sep='\t', index=False)