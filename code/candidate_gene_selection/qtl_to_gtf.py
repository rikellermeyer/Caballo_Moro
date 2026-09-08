#!/users/9/kell3262/miniforge3/envs/CMpop/bin/python

### Match QTL regions from Wiese et al., 2024 to gtf file to get number of genes in file

import pandas as pd
import janitor
import numpy as np
import yaml
import csv

####
#Load config settings
with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome = config["GENOME"]
genome_path = config["GENOME_PATH"]

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/candidate_gene_selection'

data_path = config["DATA_PATH"]
output_path_prefix = f'{data_path}/candidate_genes/QTLs'

def filter_qtl():
    #read in original QTL csv file
    raw_qtl = pd.read_csv(f'{output_path_prefix}/original_QTL_analysis/wiese_qtl_supp5.csv')
    #print(raw_qtl, raw_qtl.columns)

    #print(raw_qtl['Trait'].unique())

    eye_filter = raw_qtl[raw_qtl['Category'] == 'Eye'] 

    eye_filter.to_csv(f'{output_path_prefix}/original_QTL_analysis/wiese_eye_qtl_total.csv', index = False)

    #print(eye_filter, eye_filter['Trait'].unique())

    eye_filter = eye_filter.dropna(subset = ['Chromosome', 'Start', 'Stop'])

    #print(eye_filter)
    eye_filter.to_csv(f'{output_path_prefix}/original_QTL_analysis/wiese_eye_qtl_no_na.csv', index = False)

    #['trait', 'CHR', 'qtlStart', 'qtlStop']
    final_df = eye_filter[['Trait', 'Chromosome', 'Start', 'Stop']]

    final_df.to_csv(f'{output_path_prefix}/original_QTL_analysis/copywiese_eye_interval.txt', index = False, sep = ' ')


def janitor_join():

    gtf_col_names = ['CHR', 'geneStart', 'geneStop', 'geneID']
    gtf_site_df = pd.read_csv(gtf_file, skiprows=1, names = gtf_col_names) 


    qtl_col_names = ['trait', 'CHR', 'qtlStart', 'qtlStop']
    qtl_df = pd.read_csv(qtl_file, skiprows=1, sep = '\s', names = qtl_col_names)

    #remove MelEyes
    qtl_df = qtl_df[qtl_df['trait'] != 'MelEyes']

    qtl_to_gene_df = gtf_site_df.conditional_join(qtl_df, 
                                          ('geneStart', 'qtlStart', '>='),
                                          ('geneStop', 'qtlStop', '<='),
                                          ('CHR', 'CHR', '=='),
                                          how='inner')

    print(qtl_to_gene_df)

    qtl_to_gene_df = qtl_to_gene_df.drop(columns=[('right', 'CHR')]).dropna()
    qtl_to_gene_df.columns = qtl_to_gene_df.columns.get_level_values(1)

    qtl_to_gene_df.to_csv(f'{output_path_prefix}/QTL_gene_matches.txt', index = False)
    print(qtl_to_gene_df)


    uniq_output = qtl_to_gene_df.drop_duplicates(subset = 'geneID', keep = 'first')
    print(uniq_output)

    uniq_output['geneID'].to_csv(f'{output_path_prefix}/QTL_geneIDs_only.txt', index = False)



if __name__ == '__main__':
    qtl_file = f'{output_path_prefix}/original_QTL_analysis/wiese_eye_intervals.txt'
    gtf_file = f'{genome_path}/gtf_new_chrom.txt' # made from `gwas_to_gtf.py`

    #Filter the raw QTL file from Supplemental Table 5 - Wiese
    filter_qtl()

    #Match QTL regions to gene names
    janitor_join()