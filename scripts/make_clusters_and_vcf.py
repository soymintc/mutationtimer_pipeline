import argparse
import pandas as pd
from collections import Counter

def parse_args():
    description = "Create 'clusters' object input for MutationTimeR"
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--labelled', help="SNVs labelled with cluster info")
    p.add_argument('--purityploidy', help="Purity and ploidy info file")
    p.add_argument('--vcf_header', help="VCF header file")
    p.add_argument('--out_clusters', help="Output file with cluster proportions")
    p.add_argument('--out_vcf', help="Output VCF t_ref_count, t_alt_count in INFO")
    args = p.parse_args()
    return args

def get_purity_ploidy(purityploidy):
    pp = pd.read_csv(purityploidy) # purity,ploidy file
    purity = pp['purity'].values[0]
    ploidy = pp['ploidy'].values[0]
    return purity, ploidy

def get_labelled_table(labelled_table):
    lab = pd.read_table(labelled_table, 
        dtype={'chrom':str, 'coord':int, 'ref_counts':int, 'alt_counts':int}) # labelled
    print(f'[LOG] before filtering: {lab.shape}')
    lab = lab[lab['cluster'].notnull()]
    print(f'[LOG] after removing SNVs with no clusters mapped: {lab.shape}')
    return lab

def count_clusters(lab, purity):
    cluster_cnts = Counter()
    for rix, row in lab.iterrows():
        clusters = row['cluster'].split(',')
        cluster_cnts.update(Counter(clusters))

    cls = [] # cluster
    ssms = [] # n_ssms
    props = [] # proportion
    for ix, (cluster, cluster_cnt) in enumerate(cluster_cnts.most_common()):
        cls.append(ix)
        ssms.append(cluster_cnt)
    clustdf = pd.DataFrame({
        'cluster': cls,
        'n_ssms': ssms, 
    })
    clustdf['proportion'] = purity * (clustdf['n_ssms'] / clustdf['n_ssms'].sum())
    return clustdf

def write_vcf(lab, vcf_header, out_vcf):
    with open(out_vcf, 'w') as out:
        for line in open(vcf_header, 'r'):
            if line.startswith('##'): 
                out.write(line)
                continue
            elif line.startswith('#'):
                out.write(line)
                field = line.rstrip().split('\t')
                sample_ids = field[9:]
                continue

        for rix, row in lab.iterrows():
            chrom, pos, ref, alt = rix # ('X', 154022392, 'C', 'T')
            t_ref_count = row['ref_counts']
            t_alt_count = row['alt_counts']
            t_depth = t_ref_count + t_alt_count
            t_vaf = t_alt_count / t_depth
            is_hom = (t_vaf >= 0.9)
            gt = '1/1' if is_hom else '0/1'
            ID = '.'
            QUAL = '.'
            FILTER = 'PASS'
            INFO = f't_alt_count={t_alt_count};t_ref_count={t_ref_count}'
            FORMAT = 'GT:AD:DP'
            SAMPLE = f'{gt}:{t_alt_count}:{t_depth}'
            field = [chrom, pos, ID, ref, alt, QUAL, FILTER, INFO, FORMAT, SAMPLE]
            field = [str(x) for x in field]
            line = '\t'.join(field) + '\n'
            out.write(line)

if __name__ == '__main__':
    INDEX = ['chrom', 'coord', 'ref', 'alt']

    args = parse_args()

    purity, ploidy = get_purity_ploidy(args.purityploidy)

    lab = get_labelled_table(args.labelled)
    lab.set_index(INDEX, inplace=True)
    print(f'[LOG] lab.head(): \n{lab.head()}')
    clustdf = count_clusters(lab, purity)
    clustdf.to_csv(args.out_clusters, sep='\t', index=False)
    
    write_vcf(lab, args.vcf_header, args.out_vcf)
