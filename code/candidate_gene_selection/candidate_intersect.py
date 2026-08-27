#!/users/9/kell3262/miniforge3/envs/CMpop/bin/python


############################################################################
### This script takes in "gene hits" from various avenues and conducts total
### and pairwise comparisons.
### Gene hit origins: 
###### genotype filtering - see `variant_calling/variant_filter.py` & `snpeff_annotation.py`
###### GWAS - see `popgen/gwas_fixandtest.py`
###### QTL - Wiese et al., 2024
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


#First, process Rachel's sweeps file

def process_sweeps_file():
    raw_column_names = ['ensemble', 'geneID', 'Scaffold', 'Start', 'End',
       'Pachon', 'Tinaja',
       'Yerbaniz', 'Molino',
       'Vasquez', 'Jineo',
       'Escondido', 'Peroles',
       'Rascon', 'Mante',
       'Choy']
    raw_sweeps_df = pd.read_csv(sweeps_file, skiprows = 1, header = None, names = raw_column_names)
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

    sweeps_output_file = f'{output_path_prefix}/sweeps/Moran_hard_and_soft_sweeps.csv'
    print(f'Sweeps file saved to: {sweeps_output_file}')
    sweeps_hard_geneID_only.to_csv(sweeps_output_file, index=False)    

    return sweeps_hard_geneID_only
   

def make_df():
    gwas_col_names = ['geneID', 'origin']
    gwas_and_geno_hits = pd.read_csv(gwas_and_geno_file, sep= ' ', names = gwas_col_names)

    qtl_col_names = ['geneID']
    qtl_hits = pd.read_csv(qtl_file, sep = '\t', names = qtl_col_names)
    qtl_hits['origin']='QTL'

    zf_col_names = ['db_gene_id','description', 'geneID']
    zf_go_terms_hits = pd.read_csv(zf_db_file, sep = '\t', names = zf_col_names)
    zf_go_terms_hits['origin']='GOterm'

    sweeps_hits = process_sweeps_file()

    all_candidates_int = pd.concat([gwas_and_geno_hits, qtl_hits, zf_go_terms_hits, sweeps_hits])
    
    #print(all_candidates_int)

    all_candidates = all_candidates_int.iloc[:,0:2]

    all_candidates.to_csv(f'{output_path_prefix}/list_of_all_candidate_genes.csv', index = False)

    return all_candidates

def shared_genes():
    all_candidates = make_df()
    # Get unique values for each group as sets
    group_sets = all_candidates.groupby('origin')['geneID'].apply(lambda x: set(x.unique()))

    # Find the intersection of all sets
    shared_genes = set.intersection(*group_sets.values)
    print("Shared genes across all groups:", shared_genes)
    print(f'{len(shared_genes)} shared across all groups')

def unique_genes():
    all_candidates = make_df()
    group_sets = all_candidates.groupby('origin')['geneID'].apply(lambda x: set(x.unique()))

    unique_to_group={}

    for group_name, current_set in group_sets.items():
        # Create a set of all other group's values
        other_groups_values = set.union(*(s for name, s in group_sets.items() if name != group_name))
        
        # Find elements unique to the current group
        unique_to_group[group_name] = current_set - other_groups_values

    print("\nGenes unique to each group:")
    for group, genes in unique_to_group.items():
        print(f"Group {group}: {len(genes)}\n")

def pairwise_genes():
    all_candidates = make_df()
    groups = all_candidates["origin"].unique()
    pairwise = list(combinations(groups, 2))

    # Count shared genes
    shared_counts = []
    for g1, g2 in pairwise:
        genes_g1 = set(all_candidates[all_candidates["origin"] == g1]["geneID"])
        genes_g2 = set(all_candidates[all_candidates["origin"] == g2]["geneID"])
        
        shared = genes_g1 & genes_g2
        shared_counts.append((g1, g2, len(shared)))
        unique_to_g1 = genes_g1 - genes_g2
        unique_to_g2 = genes_g2 - genes_g1

        plt.figure(figsize=(4, 4))
        venn2([genes_g1, genes_g2], set_labels=(g1, g2))
        plt.title(f"{g1} vs {g2} Shared Genes")
        plt.tight_layout()
        plt.savefig(f"./candidate_gene_images/{g1}v{g2}_genes.png")


    # Build dictionary of sets
    group_sets = {
        group: set(all_candidates[all_candidates["origin"] == group]["geneID"])
        for group in all_candidates["origin"].unique()
    }

    # Plot 4-set Venn diagram
    venn(group_sets)
    plt.title("5-Group Gene Overlap")
    plt.savefig("candidate_gene_images/5_group_venn.png", dpi=300, bbox_inches="tight")
    plt.savefig("candidate_gene_images/5_group_venn.svg", dpi=300, bbox_inches="tight", format = "svg")
    plt.close()
    



if __name__ == '__main__':
    gwas_and_geno_file = f'{output_path_prefix}/gwas_and_geno_filter_hits.txt'
    qtl_file = f'{output_path_prefix}/QTLs/genes_in_eye_qtl_list.txt'
    zf_db_file = f'{output_path_prefix}/ZF_DB_eye_genes.txt'
    sweeps_file = f'{output_path_prefix}/sweeps/Moran_sweeps_supp.csv'

    process_sweeps_file()

    #make_df()

    #shared_genes()

    #unique_genes()

    #pairwise_genes()