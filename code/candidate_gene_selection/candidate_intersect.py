#!/users/9/kell3262/miniforge3/envs/CMpop/bin/python


############################################################################
### This script takes in "gene hits" from various avenues and conducts total
### and pairwise comparisons.
### Gene hit origins: 
###### genotype filtering - see `variant_calling/variant_filter.py` & `snpeff_annotation.py`
###### GWAS - see `popgen/gwas_fixandtest.py` > `filter_gwas.sh` > 
###### QTL - Wiese et al., 2024 > `qtl_to_gtf.py`
###### Zebrafish GO Terms - AmiGO2 & PANTHER
###### Selective sweeps - Moran et al., 2023
### Identifies shared and unique alleles
### Makes a venn diagram of shared alleles
############################################################################

import yaml
import pandas as pd
from itertools import combinations
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships
from venn import venn

####
#Load config settings
with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome = config["GENOME"]
genome_path = config["GENOME_PATH"]

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/candidate_gene_selection'

data_path = config["DATA_PATH"]
output_path_prefix = f'{data_path}/candidate_genes'

report_path = config["REPORT_PATH"]
report_output_prefix = f'{report_path}/candidate_genes_venn'

#First, process Rachel's sweeps file
def process_sweeps_file(raw_sweeps_file):
    raw_column_names = ['ensemble', 'geneID', 'Scaffold', 'Start', 'End',
       'Pachon', 'Tinaja',
       'Yerbaniz', 'Molino',
       'Vasquez', 'Jineo',
       'Escondido', 'Peroles',
       'Rascon', 'Mante',
       'Choy']
    raw_sweeps_df = pd.read_csv(raw_sweeps_file, skiprows = 1, header = None, names = raw_column_names)
    #print(raw_sweeps_df.columns)
    #print(raw_sweeps_df.iloc[0:3,])

    #Remove extra columns:
    # keep Gene_Symbol (rename to geneID)
    # Filter to only caves
    cave_list = ['ensemble','Scaffold','Start','End',
        'Peroles','Rascon', 'Mante', 'Choy']
    cave_sweeps = raw_sweeps_df.drop(cave_list, axis = 1).dropna(subset = 'geneID', inplace = False)
    #print(cave_sweeps.columns)
    #print(cave_sweeps.iloc[0:3,])

    # Filter to only genes with hard or soft sweeps in caves
    hard_or_soft_sweeps = cave_sweeps[cave_sweeps.apply(
        lambda row: row.astype(str).str.contains(
            'hard|soft', case=False, na=False).any(), axis=1)]

    #print(hard_or_soft_sweeps.iloc[0:3,0:16])

    ## Reformat to columns: geneID, hard, soft where hard and soft contain a list of the caves that match
    def extract_conditions(row):
        hard_cols = []
        soft_cols = []
        for col in hard_or_soft_sweeps.columns[1:]:  # Skip 'geneID'
            #print(col)
            val = row[col]
            if pd.notna(val):
                if 'hard' in val:
                    hard_cols.append(col)
                if 'soft' in val:
                    soft_cols.append(col)
        return pd.Series({
            'geneID': row['geneID'],
            'hard': ', '.join(hard_cols),
            'soft': ', '.join(soft_cols)
        })

    # Apply the function row-wise
    sweeps_by_cave = hard_or_soft_sweeps.apply(extract_conditions, axis=1)
    intermediate_sweeps_file = f'{output_path_prefix}/sweeps/Moran_et_al_2023_sweeps_in_at_least_one_cave.tsv'
    sweeps_by_cave.to_csv(intermediate_sweeps_file, sep = '\t', index=False)
    #print(sweeps_by_cave[sweeps_by_cave['geneID'].str.contains('gja8b', case = False, na = False)])


    sweeps_hard_geneID_only = sweeps_by_cave[sweeps_by_cave['hard'] != ''][['geneID']]
    sweeps_hard_geneID_only['origin']='Sweeps'
    #print(sweeps_hard_geneID_only)

    
    print(f'Sweeps file saved to: {sweeps_file}')
    sweeps_hard_geneID_only.to_csv(sweeps_file, index=False)    

    return sweeps_hard_geneID_only
   

def make_df():
    #read and name various inputs
    col_names = ['geneID', 'origin']

    geno_filter_df = pd.read_csv(geno_filter_file, comment = '#', names = col_names)
    geno_filter_df['origin'] = 'geno_filter'


    gwas_df = pd.read_csv(gwas_file, names = col_names)
    gwas_df['origin'] = 'GWAS'


    qtl_hits = pd.read_csv(qtl_file, names = col_names)
    qtl_hits['origin']='QTL'


    zf_col_names = ['db_gene_id','description', 'geneID']
    zf_go_terms_hits = pd.read_csv(zf_db_file, sep = '\t', names = zf_col_names)
    zf_go_terms_hits['origin']='GOterm'
    GO_terms = zf_go_terms_hits[['geneID', 'origin']]


    sweeps_hits = pd.read_csv(sweeps_file, header = 0)

    list_of_dfs = gwas_df, geno_filter_df, qtl_hits, GO_terms, sweeps_hits
    all_candidates = pd.concat([gwas_df, geno_filter_df, qtl_hits, GO_terms, sweeps_hits])
    all_candidates.to_csv(total_file, index = False)

    summary = []
    for df in list_of_dfs:
        name = df['origin'].loc[0]
        length_of_df = len(df['geneID'])
        summary.append(f'{name}: {length_of_df}')

    return summary

def shared_genes():
    all_candidates = pd.read_csv(total_file, header = 0)
    # Get unique values for each group as sets
    group_sets = all_candidates.groupby('origin')['geneID'].apply(lambda x: set(x.unique()))

    # Find the intersection of all sets
    shared_genes = set.intersection(*group_sets.values)
    #print("Shared genes across all groups:", shared_genes)
    #print(f'{len(shared_genes)} shared across all groups')
    with open(f'{output_path_prefix}/tempintersect_all.txt', 'w') as file:
        file.write('\n'.join(shared_genes))

    summary = list(shared_genes)
    summary.append(f'Total number: {len(shared_genes)}')
    return summary

def unique_genes():
    all_candidates = pd.read_csv(total_file, header = 0)
    group_sets = all_candidates.groupby('origin')['geneID'].apply(lambda x: set(x.unique()))

    unique_to_group={}

    for group_name, current_set in group_sets.items():
        # Create a set of all other group's values
        other_groups_values = set.union(*(s for name, s in group_sets.items() if name != group_name))
        
        # Find elements unique to the current group
        unique_to_group[group_name] = current_set - other_groups_values

    summary = []
    for group, genes in unique_to_group.items():
        summary.append(f"{group}: {len(genes)}")
        file_name = f'{output_path_prefix}/unique_to_{group}.txt'
        with open(file_name, 'w') as file:
            file.write('#'+group + '\n' + '\n'.join(genes))
    return summary


def pairwise_genes(graph_y_or_no):

    all_candidates = pd.read_csv(total_file, header = 0)
    groups = all_candidates["origin"].unique()
    pairwise = list(combinations(groups, 2))

    # Count shared genes
    shared_counts = []
    summary = []
    for g1, g2 in pairwise:
        genes_g1 = set(all_candidates[all_candidates["origin"] == g1]["geneID"])
        genes_g2 = set(all_candidates[all_candidates["origin"] == g2]["geneID"])
        
        shared = genes_g1 & genes_g2
        shared_counts.append((g1, g2, len(shared)))
        unique_to_g1 = genes_g1 - genes_g2
        unique_to_g2 = genes_g2 - genes_g1

        summary.append(f'{g1} vs. {g2}: {len(shared)}')

        with open(f'{output_path_prefix}/tempshared_{g1}v{g2}.txt', 'w') as file:
            file.write('\n'.join(shared))

        if graph_y_or_no == 'yes':
            plt.figure(figsize=(4, 4))
            venn2([genes_g1, genes_g2], set_labels=(g1, g2))
            plt.title(f"{g1} vs {g2} Shared Genes")
            plt.tight_layout()
            plt.savefig(f"{report_output_prefix}/{g1}v{g2}_genes.png")


    # Build dictionary of sets
    group_sets = {
        group: set(all_candidates[all_candidates["origin"] == group]["geneID"])
        for group in all_candidates["origin"].unique()
    }

    #print(group_sets)

    if graph_y_or_no == 'yes':
        # Plot 4-set Venn diagram
        venn(group_sets)
        plt.title("5-Group Gene Overlap")
        plt.savefig(f"{report_output_prefix}/5_group_venn.png", dpi=300, bbox_inches="tight")
        plt.savefig(f"{report_output_prefix}/5_group_venn.svg", dpi=300, bbox_inches="tight", format = "svg")
        plt.close()

    #print(summary)
    return summary


if __name__ == '__main__':
    #input files
    gwas_file = f'{output_path_prefix}/GWAS/GWAS_genes_only.txt'
    geno_filter_file = f'{output_path_prefix}/geno_filter/geno_filter_gene_id_only.txt'
    gwas_and_geno_file = f'{output_path_prefix}/gwas_and_geno_filter_hits.txt'
    #qtl_file = f'{output_path_prefix}/QTLs/genes_in_eye_qtl_list.txt'
    qtl_file = f'{output_path_prefix}/QTLs/QTL_geneIDs_only.txt'
    zf_db_file = f'{output_path_prefix}/GOTerm/ZF_DB_eye_genes.txt'
    raw_sweeps_file = f'{output_path_prefix}/sweeps/Moran_sweeps_supp.csv'

    ### Run once to process sweeps file ###
    #process_sweeps_file(raw_sweeps_file)

    #intermediate file names
    sweeps_file = f'{output_path_prefix}/sweeps/Moran_hard_and_soft_sweeps.csv'
    total_file = f'{output_path_prefix}/list_of_all_candidate_genes.csv'

    # A summary dictionary that lays out all the numbers!
    summary_dict = {}

    summary_dict['Number of genes per dataset:'] = make_df()

    summary_dict['Number of unique genes per group:'] = unique_genes()
    
    summary_dict['Number of shared genes, pairwise:'] = pairwise_genes('no')

    summary_dict['Genes shared in all datasets:'] = shared_genes()


    with open(f'{output_path_prefix}/tempIntersection_Summary.txt', 'w') as file:
        for key, value in summary_dict.items():
            file.write(f'###{key}\n{'\n'.join(value)}\n\n')
