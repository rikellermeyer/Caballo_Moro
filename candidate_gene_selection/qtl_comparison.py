#!/home/rk2643/miniforge3/envs/CMpop/bin/python

#Take the list of eye term QTLs from Suzanne
#Compare to the genes captured in GWAS

import pandas as pd
import numpy as np
import yaml
import sys



#Load config settings
with open("../../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

data_path = config['DATA_PATH']
gwas_input_path = f'{data_path}/molino2/popgen/gwas/tests/goterm'

input_file = sys.argv[1]
output_file = sys.argv[2]

#read in candidate genes table from GWAS
candid_col_names =['chrom', 'snp','start', 'end', 'gene']
candid_df = pd.read_csv(input_file, sep='\s+', header =0, names=candid_col_names)
candid_df[['chrom','snp','start','end']]=candid_df[['chrom','snp','start','end']].apply(pd.to_numeric)
candid_df =candid_df.sort_values(by=['chrom', 'start'])
print(f'Unique genes in candidates: {candid_df.drop_duplicates(subset=['gene'])}')
#print(candid_df.dtypes, candid_df)

#read in QTL reference document
qtl_site_names =['id', 'chrom', 'start', 'end']
qtl_site_df=pd.read_csv('qtl_for_translate_input.csv', sep='\s+', header=0, names=qtl_site_names)
qtl_site_df=qtl_site_df.sort_values(by=['chrom','start'])
#print(qtl_site_df.dtypes, qtl_site_df)
qtl_compare = qtl_site_df[['chrom', 'start', 'end']]
#print(qtl_compare)

#Manually define molino-only QTLs
molino_specific_qtls = {'4':['4820000','5190000'], '5':['1280000', '1930000'], '24':['370000', '1850000']}
molino_df = pd.DataFrame.from_dict(molino_specific_qtls, orient='index', columns=['start', 'end']).reset_index(names='chrom')
molino_df[['chrom','start','end']]=molino_df[['chrom','start','end']].apply(pd.to_numeric)
#print(molino_df.dtypes)

all_qtls = pd.concat([qtl_compare,molino_df], ignore_index=True)
#print(all_qtls)

merge_start = pd.merge(candid_df,all_qtls,how='inner',left_on=['chrom', 'start'], right_on=['chrom','start'])
#print(merge_start)

merge_end = pd.merge(candid_df,all_qtls, how='inner',left_on=['chrom', 'end'], right_on=['chrom','end'])
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

matching_snps_here.to_csv(output_file, index=False)

non_matching_snps_here.to_csv(f'nonmatching_{output_file}', index=False)
