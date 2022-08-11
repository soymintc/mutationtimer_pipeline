import pandas as pd
import vcf
import click
import logging

@click.command()
@click.option('--in_maf', required=True,
        help="Input WGS-SOMATICCALLING maf")
@click.option('--out_maf', required=True,
        help="Output maf path")
def filter_maf(in_maf, out_maf):
    """ Filter maf using tempo-like filters 
    """
    logging.basicConfig(level=logging.DEBUG)
    dtypes = {
        'Chromosome': str,
    }

    maf = pd.read_csv(
        in_maf, # consensus somatic maf
        sep='\t', comment='#', dtype=dtypes, low_memory=False)
    logging.debug(f'maf shape just after read_csv: {maf.shape}')
    maf = maf.fillna(0)

    maf['chrom'] = maf['Chromosome']
    maf['coord'] = maf['Start_Position']
    maf['ref'] = maf['Reference_Allele']
    maf['alt'] = maf['Tumor_Seq_Allele2']

    # maf filters
    maf = maf[maf['t_depth'] > 20] # tumor_depth > 20
    logging.debug(f'maf shape after t_depth filter: {maf.shape}')
    maf = maf[maf['t_alt_count'] > 3] # Tumor_count > 3
    logging.debug(f'maf shape after t_alt_count filter: {maf.shape}')
    maf = maf[maf['n_depth'] > 10] # Normal_depth > 10
    logging.debug(f'maf shape after n_depth filter: {maf.shape}')
    maf = maf[maf['n_ref_count'] > 3] # Normal_count > 3
    logging.debug(f'maf shape after n_ref_count filter: {maf.shape}')
    maf = maf[maf['AF'] < 0.01] # AF-filter # Gnomad_allele_frequency[‘non_cancer_AF_popmax’]  < .01
    logging.debug(f'maf shape after AF filter: {maf.shape}')

    # Filter SNPs not predicted by both callers
    maf = maf.drop(['chrom', 'coord', 'ref', 'alt'], axis=1)

    maf.to_csv(out_maf, sep='\t', index=False)


if __name__ == '__main__':
    filter_maf()



