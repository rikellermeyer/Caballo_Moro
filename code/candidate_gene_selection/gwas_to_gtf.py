#!/users/9/kell3262/miniforge3/envs/CMpop/bin/python

#Take the GWAS trend test results and get gene info


import pandas as pd
import numpy as np
import yaml
import sys
import janitor
import re
import yaml

####
#Load config settings
with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome = config["GENOME"]
genome_path = config["GENOME_PATH"]

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/candidate_gene_selection'

output_path_prefix = config["DATA_PATH"]
data_path = f'{output_path_prefix}/candidate_genes/GWAS'


#read in GWAS trend
candid_col_names =['CHR', 'SNP','A1', 'A2', 'TEST', 'AFF', 'UNAFF', 'CHISQ', 'DF', 'P']
candid_df = pd.read_csv(f'{data_path}/CA_trend.filtered.tsv', sep='\s+', header =0, names=candid_col_names)

#current SNP column format: `CHR:BP:`, change to `BP` only
candid_df['SNP'] = candid_df['SNP'].str.replace(r'\d*:(\d*):', r'\1', regex = True)
candid_df[['CHR','SNP']]=candid_df[['CHR','SNP']].apply(pd.to_numeric)
candid_df =candid_df.sort_values(by=['CHR', 'SNP'])

#print(candid_df.dtypes, candid_df)

#read in gtf reference
gtf_col_names =['acc_chr', 'source', 'feature', 'start', 'end', 'score', 'strand','frame', 'attribute']
gtf_site_df=pd.read_csv(f'{genome_path}/{genome}.gtf', sep='\t', comment = '#', names=gtf_col_names)


def find_gene(attribute):
    match_geneID = re.match(r'^gene_id "(.*?)"', str(attribute)).group(1)

    return match_geneID

#print(gtf_site_df['attribute'].loc[0])
#print(re.match(r'^gene_id "(.*?)"',gtf_site_df['attribute'].loc[0]).group(1))


gtf_site_df['gene_id']=gtf_site_df['attribute'].apply(find_gene)
#print(gtf_site_df)

gtf_site_df=gtf_site_df[~gtf_site_df['acc_chr'].astype(str).str.startswith('NW')]
gtf_site_df[['start','end']]=gtf_site_df[['start','end']].apply(pd.to_numeric)
gtf_site_df=gtf_site_df.sort_values(by=['acc_chr','start'])

#reduce gtf to genes only
gtf_site_df = gtf_site_df[gtf_site_df["feature"] == 'gene']
#print(gtf_site_df['acc_chr'])


#fix gtf chr > normal chr
recode_col_names = ['CHR','acc_chr']
recode_dict = pd.read_csv(f'{genome_path}/chr_key.txt', sep = '\t', names = recode_col_names)
 

gtf_site_df = gtf_site_df.merge(recode_dict, on = 'acc_chr', how = 'left')

gtf_export = gtf_site_df[['CHR', 'start', 'end', 'gene_id']]
gtf_export.to_csv(f'{genome_path}/gtf_new_chrom.txt', index = False)

exit()

print(f'candid_df, {candid_df.shape}:\n {candid_df.head(2)} \n {candid_df.tail(2)}')
print(f'gtf_site_df, {gtf_site_df.shape}:\n {gtf_site_df.head(2)} \n {gtf_site_df.tail(2)}')


trendsnps_to_gene_df = candid_df.conditional_join(gtf_site_df, 
                                          ('SNP', 'start', '>='),
                                          ('SNP', 'end', '<='),
                                          ('CHR', 'CHR', '=='),
                                          how='left')


trendsnps_to_gene_df = trendsnps_to_gene_df.drop(columns=[('right', 'CHR')])
trendsnps_to_gene_df.columns = trendsnps_to_gene_df.columns.get_level_values(1)

captured_genes_only = trendsnps_to_gene_df['gene_id'].dropna().unique()
print(captured_genes_only, len(captured_genes_only))
with open(f'{data_path}/GWAS_genes_only.txt', 'w') as file:
    file.write('\n'.join(captured_genes_only) +'\n')


trendsnps_to_gene_df.to_csv(f'{data_path}/GWAS_to_genes_withNA_total.txt', index=False)
trendsnps_to_gene_df.dropna().to_csv(f'{data_path}/GWAS_to_genes.txt', index = False)

sites_without_genes = trendsnps_to_gene_df.loc[trendsnps_to_gene_df['gene_id'].isna()]
sites_without_genes.to_csv(f'{data_path}/GWAS_to_no_genes.txt', index = False)



"""Not sure about this stuff 
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
"""