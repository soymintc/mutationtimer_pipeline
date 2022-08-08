import pandas as pd
import vcf
import fire


def filter_vcf(consensus_somatic_maf, output_maf):
    """ Filter maf using tempo-like filters 
    """
    dtypes = {
        'Chromosome': str,
    }

    maf = pd.read_csv(
        consensus_somatic_maf,
        sep='\t', comment='#', dtype=dtypes, low_memory=False)
    maf = maf.fillna(0)

    maf['chrom'] = maf['Chromosome']
    maf['coord'] = maf['Start_Position']
    maf['ref'] = maf['Reference_Allele']
    maf['alt'] = maf['Tumor_Seq_Allele2']

    # maf filters
    maf = maf[maf['t_depth'] > 20] # tumor_depth > 20
    maf = maf[maf['t_alt_count'] > 3] # Tumor_count > 3
    maf = maf[maf['n_depth'] > 10] # Normal_depth > 10
    maf = maf[maf['n_ref_count'] > 3] # Normal_count > 3
    maf = maf[maf['AF'] < 0.01] # AF-filter # Gnomad_allele_frequency[‘non_cancer_AF_popmax’]  < .01

    # Filter SNPs not predicted by both callers
    maf = maf.drop(['chrom', 'coord', 'ref', 'alt'], axis=1)

    maf.to_csv(output_maf, sep='\t', index=False)


if __name__ == '__main__':
    fire.Fire(filter_vcf)



