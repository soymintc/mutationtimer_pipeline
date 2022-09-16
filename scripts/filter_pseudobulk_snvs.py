import time
from collections import defaultdict
from argparse import ArgumentParser
import numpy as np
import pandas as pd


def get_args():
    desc = 'Extract and filter pseudobulk SNV counts'
    p = ArgumentParser(description=desc)

    p.add_argument('museq_vcf', help="museq vcf.gz file")
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
            print("[ERROR] chrom_snvs['coord']", chrom_snvs['coord'])
            print("[ERROR] chrom_snvs['coord'].shape", chrom_snvs['coord'].shape)
            print("[ERROR] chrom_blacklist", chrom_blacklist)
            print("[ERROR] chrom_blacklist.shape", chrom_blacklist.shape)
        ok_idx.extend(chrom_snvs.index[~is_low_mapp])

    return snvs.loc[ok_idx]


# baed on scgenome
def filter_snvs(museq, strelka, mappability, cosmic, dbsnp, snpeff, trinuc,
        blacklist, remove_chr6p):

    print(f"[LOG] museq.shape: {museq.shape}")
    museq = museq.rename(columns={'score': 'museq_score'})
    museq = museq[museq['museq_score'] > 0.9]
    print(f"[LOG] museq.shape after score filtering: {museq.shape}")

    print(f"[LOG] museq.shape: {museq.shape}")
    if remove_chr6p:
        is_chr6p = (museq['chrom'] == '6') & (museq['coord'] < 65e5)
        museq = museq[~is_chr6p]
    print(f"[LOG] museq.shape after chr6 filtering: {museq.shape}")

    print(f"[LOG] strelka.shape: {strelka.shape}")
    strelka = strelka.rename(columns={'score': 'strelka_score'})
    strelka = strelka[strelka['strelka_score'] > 20.]
    print(f"[LOG] strelka.shape after chr6 filtering: {strelka.shape}")

    print(f"[LOG] mappability.shape: {mappability.shape}")
    mappability = mappability[mappability['mappability'] > 0.99]
    print(f"[LOG] mappability.shape after mappability filtering: {mappability.shape}")

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

    print(f"[LOG] museq.shape: {museq.shape}")
    snvs = museq.merge(strelka)
    print(f"[LOG] museq.shape after strelka merge: {museq.shape}")
    snvs = snvs.merge(mappability)
    print(f"[LOG] museq.shape after mappability merge: {museq.shape}")
    snvs = snvs.merge(cosmic, how='left')
    print(f"[LOG] museq.shape after cosmic merge: {museq.shape}")
    snvs = snvs.merge(dbsnp, how='left')
    print(f"[LOG] museq.shape after dbsnp merge: {museq.shape}")
    snvs = snvs.merge(snpeff, how='left')
    print(f"[LOG] museq.shape after snpeff merge: {museq.shape}")
    snvs = snvs.merge(trinuc, how='left')
    print(f"[LOG] museq.shape after trinuc merge: {museq.shape}")

    snvs['is_cosmic'] = snvs['is_cosmic'].fillna(False)
    snvs['is_dbsnp'] = snvs['is_dbsnp'].fillna(False)

    snvs = snvs[~snvs['is_dbsnp']]
    print(f"[LOG] museq.shape after is_dbsnp filtering: {museq.shape}")

    snvs = filter_blacklist(snvs, blacklist)
    snvs.drop_duplicates(inplace=True)
    print(f"[LOG] museq.head() after drop_duplicates filtering: {museq.head()}")
    print(f"[LOG] museq.shape after drop_duplicates filtering: {museq.shape}")

    snvs.set_index(['chrom', 'coord', 'ref', 'alt'], inplace=True)

    return snvs

def process_vcf(museq_vcf):
    """Return (chrom, pos, ref, alt)-indexed (t_ref_count, t_alt_count, n_ref_count, n_alt_count) dict
    """
    vcf_header = ['chrom', 'coord', 'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'format', 'tumour', 'normal']
    vcf = pd.read_table(museq_vcf, comment='#', names=vcf_header, dtype=str, low_memory=False)
    vcf = vcf[['chrom', 'coord', 'ref', 'alt', 'format', 'tumour', 'normal']]
    vcfdata = defaultdict(list)
    for rix, row in vcf.iterrows():
        chrom = row['chrom']
        coord = int(row['coord'])
        ref = row['ref']
        alt = row['alt']
        vcfdata['chrom'].append(chrom)
        vcfdata['coord'].append(coord)
        vcfdata['ref'].append(ref)
        vcfdata['alt'].append(alt)
        ##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Count with reference to REF">
        ##FORMAT=<ID=AC,Number=1,Type=Integer,Description="Count with reference to ALT">
        ##FORMAT=<ID=NI,Number=1,Type=Integer,Description="Number of reads with an insertion">
        ##FORMAT=<ID=ND,Number=1,Type=Integer,Description="Number of reads with a Deletion">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth after filtering for Failed, duplicate and low quality reads">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Phred-scaled genotype likelihoods">
        formats = row['format'].split(':')
        t_gts = dict(zip(formats, row['tumour'].split(':')))
        n_gts = dict(zip(formats, row['normal'].split(':')))
        assert ',' not in t_gts['RC'], (ix, t_gts)
        assert ',' not in t_gts['AC'], (ix, t_gts)
        assert ',' not in n_gts['RC'], (ix, n_gts)
        assert ',' not in n_gts['AC'], (ix, n_gts)
        t_ref_count = int(t_gts['RC'])
        t_alt_count = int(t_gts['AC'])
        n_ref_count = int(n_gts['RC'])
        n_alt_count = int(n_gts['AC'])
        vcfdata['t_ref_count'].append(t_ref_count)
        vcfdata['t_alt_count'].append(t_alt_count)
        vcfdata['n_ref_count'].append(n_ref_count)
        vcfdata['n_alt_count'].append(n_alt_count)

    vcfdata = pd.DataFrame(vcfdata, columns=['chrom', 'coord', 'ref', 'alt', 
                                             't_ref_count', 't_alt_count',
                                             'n_ref_count', 'n_alt_count'])
    vcfdata.set_index(['chrom', 'coord', 'ref', 'alt'], inplace=True)
    return vcfdata


if __name__ == '__main__':
    argv = get_args()

    print(f"[LOG:{time.ctime()}] init")
    vcfdata = process_vcf(argv.museq_vcf)
    print(f"[LOG:{time.ctime()}] process_vcf complete")
    museq = pd.read_csv(argv.museq, dtype={'chrom': str})
    strelka = pd.read_csv(argv.strelka, dtype={'chrom': str})
    print(f"[LOG:{time.ctime()}] pd.read_csv complete")

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
                'tri_nucleotide_context',
                't_ref_count',
                't_alt_count',
                'n_ref_count',
                'n_alt_count',
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
            argv.blacklist, sep='\t', skiprows=1,
            names=['chr', 'start', 'end', 'width', 'desc']
        )
        print(f"[LOG:{time.ctime()}] pd.read_csv for supplementals complete")

        snvs = filter_snvs(
            museq, strelka, mappability, cosmic, dbsnp, snpeff, trinuc,
            blacklist, argv.remove_chr6p
        )
        snvs = snvs.join(vcfdata)
        print(f"[LOG:{time.ctime()}] filter_snvs complete")

    snvs.to_csv(argv.snvs, sep='\t')
