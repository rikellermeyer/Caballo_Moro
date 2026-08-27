#!/users/9/kell3262/miniforge3/envs/CMpop/bin/python

### Match QTL regions from Wiese et al., 2024 to gtf file to get number of genes in file

import pandas as pd
import janitor
import numpy as np
import yaml
import csv

####
#Load config settings
with open("../../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome = config["GENOME"]
genome_path = config["GENOME_PATH"]

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/candidate_gene_selection'

output_path_prefix = config["DATA_PATH"]
data_path = f'{output_path_prefix}/candidate_gene_selection'

### first fix chromosome name from gtf
def fix_chr_names(file_to_fix):
    key_names = ['chr', 'old_chrom']
    column_dtypes = {'chr': str, 'old_chrom': str}
    chr_match = pd.read_csv(chr_key, names = key_names, sep = '\t', dtype=column_dtypes)

    gtf_key_names = ['old_chrom', 'gene_start', 'gene_end', 'gene']
    gtf_original = pd.read_csv(file_to_fix, names = gtf_key_names, sep = '\t', dtype=column_dtypes)

    merged_chr_name_df = gtf_original.merge(chr_match, on = 'old_chrom', how = 'left')

    gtf_coord_ready = merged_chr_name_df.dropna().drop('old_chrom', axis = 1)

    gtf_list_of_tuples = list(gtf_coord_ready[['chr','gene_start','gene_end','gene']].itertuples(index=False, name=None))
    #genes = list(gene_df[["chromosome", "gene_start", "gene_end", "gene_name"]].itertuples(index=False, name=None))

    #print(gtf_list_of_tuples)

    return gtf_list_of_tuples


def simple_qtl_regions(gtf_coord_file, qtl_coord_file):
    genes = fix_chr_names(gtf_coord_file)

    # Load QTL regions
    qtl_regions = []
    with open(qtl_coord_file) as qtl_file:
        next(qtl_file) #skip the header
        for line in qtl_file:
            chrom, start, end = line.strip().split(" ")
            #print(chrom, '\t', start, '\t', end)
            #print(type(chrom), type(start), type(end))
            qtl_regions.append((chrom, int(start), int(end)))
    

    # Find genes within QTLs
    genes_in_qtls = []
    n=1
    for qtl_chrom, qtl_start, qtl_end in qtl_regions:
        #print(qtl_chrom, qtl_start, qtl_end)
        for gene_chrom, gene_start, gene_end, gene_name in genes:
            #print(gene_chrom, gene_start, gene_end, gene_name)
            if qtl_chrom == gene_chrom and gene_start >= qtl_start and gene_end <= qtl_end:
                #print(f'For {qtl_chrom}/{gene_chrom}:\n \
                #gene_start {gene_start}>= qtl_start {qtl_start} \n\
                #gene_end {gene_end} <= qtl_end {qtl_end}')
                genes_in_qtls.append((qtl_chrom, gene_start, gene_end, gene_name))
                if n==1:
                    print((qtl_chrom, gene_start, gene_end, gene_name))
                    n+=1
            
    #print(genes_in_qtls)

    # Output results
    #for gene in genes_in_qtls:
    #    print("\t".join(map(str, gene)))

    with open(f"{output_path_prefix}/genes_in_eye_qtls.txt", "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["chromosome", "gene_start", "gene_end", "gene_name"])  # Optional header
        writer.writerows(genes_in_qtls)


if __name__ == '__main__':

    qtl_coord_file = f'{data_path}/intermediate_files/qtl_eye_chr_start_stop.txt'
    gtf_coord_file = f'{genome_path}/genes_and_pos.{genome}.csv'
    chr_key = f'{genome_path}/chr_key.txt' # new \t old

    #fix_chr_names(gtf_coord_file)

    simple_qtl_regions(gtf_coord_file, qtl_coord_file)