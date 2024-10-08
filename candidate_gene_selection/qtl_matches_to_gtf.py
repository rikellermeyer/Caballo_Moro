#!/home/rk2643/miniforge3/envs/CMpop/bin/python

import re
import pandas as pd
import janitor
import numpy as np

# matching file is the chromosome, start, end, and geneID from the combined GWAS output and SnpEff output
# geneID taken from Amex3.0_surface.gtf
matching_file = './all_candidate_variants.csv'
qtl_ensembl_match = './qtl_for_translate_input.csv'


match_names=['chrom','SNP','gene_start1','gene_end1', 'gene_id','true']
match_db = pd.read_csv(matching_file, names=match_names)
print(match_db)

ensembl_labels=['ensembl','chr','region_start','region_end']
ensembl_db = pd.read_csv(qtl_ensembl_match, names=ensembl_labels, sep=' ')
#print(ensembl_db)

#ensembl_and_match=pd.merge(left=match_db, right=ensembl_db, how='left', 
#                           left_on=['CHR', 'start', 'end'], 
#                           right_on=['CHR', 'start', 'end'])


ensembl_and_match = match_db.conditional_join(ensembl_db, 
                                          ('SNP', 'region_start', '>='),
                                          ('SNP', 'region_end', '<='),
                                          ('chrom', 'chr', '=='),
                                          how='left')

#print(ensembl_and_match.keys())
#print(ensembl_and_match[ensembl_and_match['left']['gene_id'].isin(['gja8b'])])


original_qtl_info='./qtl_eye_history.csv'

og_qtl_names=['Study','Cross','Trait','ensembl','gene_biotype','chrom','gene_strand','gene_start','gene_end','gene_length','gene_QTL_overlap','QTL_start','QTL_end','QTL_peak_or_CI','LOD','PVE']
og_qtl_db=pd.read_csv(original_qtl_info, names=og_qtl_names)

#print(og_qtl_db)
#How I want it oriented:
# chr, snp, gene id, all the qtl things...


merged_matches_and_og_qtl = pd.merge(left= ensembl_and_match, 
                                     left_on= ['chrom', 'ensembl'],
                                     right= og_qtl_db,
                                     right_on= ['chrom', 'ensembl'])

#print(merged_matches_and_og_qtl.iloc[10])

#this is for dropping duplicates
ready_for_printing=merged_matches_and_og_qtl.drop(['gene_start1','gene_end1','true','chr'], axis=1)
print(ready_for_printing.keys())
drop_these_dups = ready_for_printing[['chrom','SNP','ensembl']].drop_duplicates()

#['chrom', 'SNP', 'gene_id', 'ensembl', 'region_start', 'region_end',
#       'Study', 'Cross', 'Trait', 'gene_biotype', 'gene_strand', 'gene_start',
#       'gene_end', 'gene_length', 'gene_QTL_overlap', 'QTL_start', 'QTL_end',
#       'QTL_peak_or_CI', 'LOD', 'PVE']

recombine_no_dups= drop_these_dups.join(ready_for_printing, how='inner', lsuffix='_L',rsuffix='_R').drop(['chrom_L',
                                                                                          'SNP_L',
                                                                                          'ensembl_L'], axis=1)

print(recombine_no_dups)
print(recombine_no_dups.keys())

recombine_no_dups.to_csv('merged_all_qtl_output.tsv', sep='\t',index=False)