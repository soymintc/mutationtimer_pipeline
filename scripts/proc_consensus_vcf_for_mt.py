import os
import vcf
import click
from collections import OrderedDict

@click.command()
@click.option('--in_vcf', required=True,
        help="Input vcf from consensus maf")
@click.option('--out_vcf', required=True,
        help="Output vcf path")
def proc_consensus_vcf_for_mt(in_vcf, out_vcf):
    #in_vcf = '/juno/work/shah/users/chois7/retreat/test/OV-081/results/OV-081.vcf'
    #out_vcf = '/juno/work/shah/users/chois7/retreat/test/OV-081/results/OV-081.addi.vcf' # tumor_ref_count, tumor_alt_count INFO added 
    assert os.path.exists(in_vcf)
    assert os.path.isdir(os.path.split(out_vcf)[0])

    vcf_reader = vcf.Reader(open(in_vcf, 'r'))

    # update vcf meta infos
    vcf_reader.infos['t_ref_counts'] = vcf.parser._Info(    
        't_ref_count', 1, 'Float', 'Number of reference reads',
        None, None)
    vcf_reader.infos['t_alt_counts'] = vcf.parser._Info(    
        't_alt_count', 1, 'Float', 'Number of alternate/variant reads',
        None, None)
    
    # init vcf writer after meta info update
    vcf_writer = vcf.Writer(open(out_vcf, 'w'), vcf_reader)

    # write records
    for record in vcf_reader:
        # get tumor_ref_count, tumor_alt_count
        flag_normal = False
        for sample in record.samples:
            if 'NORMAL' in sample.sample:
                flag_normal = True
                n_ref_count = sample['AD'][0]
                n_alt_count = sum(sample['AD'][1:])
                n_sample = sample
                continue
            if not flag_normal:
                t_ref_count = sample['AD'][0]
                t_alt_count = sum(sample['AD'][1:])
                t_sample = sample
                # break
                continue
        record.INFO = {'t_ref_count':t_ref_count, 't_alt_count':t_alt_count}
        vcf_writer.write_record(record)
        
    # close vcf writer
    vcf_writer.close() 


if __name__ == '__main__':
    proc_consensus_vcf_for_mt()
