#!/home/rk2643/miniforge3/envs/CMpop/bin/python

#Take the GWAS trend test results and get gene info


import pandas as pd
import numpy as np
import yaml
import sys
import janitor
import re



#read in GWAS trend
candid_col_names =['CHR', 'SNP','A1', 'A2', 'TEST', 'AFF', 'UNAFF', 'CHISQ', 'DF', 'P']
candid_df = pd.read_csv('GWAS_pfiltered_output.tsv', sep='\s+', header =0, names=candid_col_names)
candid_df[['CHR','SNP']]=candid_df[['CHR','SNP']].apply(pd.to_numeric)
candid_df =candid_df.sort_values(by=['CHR', 'SNP'])
#print(candid_df)

#print(f'Unique genes in candidates: {candid_df.drop_duplicates(subset=['gene'])}')
#print(candid_df.dtypes, candid_df)

#read in gtf reference
print('now gtf')
gtf_col_names =['acc_chr', 'source', 'feature', 'start', 'end', 'score', 'strand','frame', 'attribute']
gtf_site_df=pd.read_csv('gtf_reference.gtf', sep='\t', header=0, names=gtf_col_names)

def find_gene(attribute):
    return re.match(r'^gene_id "(.*?)"', attribute).group(1)

print(gtf_site_df['attribute'].loc[0])
print(re.match(r'^gene_id "(.*?)"',gtf_site_df['attribute'].loc[0]).group(1))

gtf_site_df['gene_id']=gtf_site_df['attribute'].apply(find_gene)
print(gtf_site_df)

gtf_site_df=gtf_site_df[~gtf_site_df['acc_chr'].astype(str).str.startswith('NW')]
gtf_site_df[['start','end']]=gtf_site_df[['start','end']].apply(pd.to_numeric)
gtf_site_df=gtf_site_df.sort_values(by=['acc_chr','start'])
#print(gtf_site_df['acc_chr'])


#fix gtf chr > normal chr
recode_col_names = ['CHR','acc_chr']
recode_dict={}
with open('chr_key.txt','r') as file:
    for line in file:
        line=line.rstrip()
        key_list=line.split(' ')
        if key_list[1] in recode_dict.keys():
            next
        else:
            recode_dict[key_list[1]]=key_list[0]    

gtf_site_df['acc_chr']=gtf_site_df['acc_chr'].replace(recode_dict).apply(pd.to_numeric)
#print(gtf_site_df)

trendsnps_to_gene_df = candid_df.conditional_join(gtf_site_df, 
                                          ('SNP', 'start', '>='),
                                          ('SNP', 'end', '<='),
                                          ('CHR', 'acc_chr', '=='),
                                          how='left')

print(trendsnps_to_gene_df['gene_id'])
trendsnps_to_gene_df.to_csv('GWAS_to_GTF_out.txt', index=False)

sys.exit()
merge_start = pd.merge(candid_df, gtf_site_df,how='inner',left_on=['chrom', 'start'], right_on=['chrom','start'])
#print(merge_start)

merge_end = pd.merge(candid_df, gtf_site_df, how='inner',left_on=['chrom', 'end'], right_on=['chrom','end'])
#print(merge_end)

uniq_cand_genes=pd.Series(candid_df['gene'].unique())

uniq_compare_start=pd.Series(merge_start['gene'].unique())

uniq_compare_end=pd.Series(merge_end['gene'].unique())


#print(f'original uniq: {uniq_cand_genes}, start uniq: {uniq_compare_start}, end uniq: {uniq_compare_end}')

matching_qtl_genes = pd.concat([uniq_compare_start, uniq_compare_end], ignore_index=True)
#print(matching_qtl_genes)



def check_range(x, code_range_df):
    #chromosome has to match
    index_check=code_range_df.chrom == x.chrom
    #print(f'cand {x.chrom}\n range: {code_range_df.chrom}')

    #fetch all matching chromosome positions from start/end table
    code_range_df=code_range_df.loc[index_check]
    #print(f'start/end positions for chrom_specified {code_range_df}')
    #print(f'start {code_range_df.start} end {code_range_df.end}')
    #print(f'snp: {x.snp}')
    #print(f'x: {x}')
    check_answer = (code_range_df.start <= x.snp) & (code_range_df.end >= x.snp)

    return (check_answer.any())

candid_df['output'] = candid_df.apply(lambda x: check_range(x, all_qtls), axis=1)

#print(candid_df)

#remove false/snps with no matching qtl, and multiples of a gene (based on multiple snps in one gene)
matching_snps_here = candid_df.loc[candid_df.output, :].drop_duplicates(subset=['gene'])
print(matching_snps_here)

non_matching_snps_here = candid_df[~candid_df["output"]].drop_duplicates(subset=['gene'])
print(non_matching_snps_here)

matching_snps_here.to_csv('GWAS_to_GTF_out.txt', index=False)

non_matching_snps_here.to_csv(f'nonmatching_GWAS_to_GTF_out.txt', index=False)
