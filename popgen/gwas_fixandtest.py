#!/home/rk2643/miniforge3/envs/CMpop/bin/python

# This takes the plink/bed files generated from `preplink.py`
# Fixes the extra files: .fam, .ped, .map
# NOTE that after running `ped_recode()` you have to run an awk script to do the replacing
# because it's super slow through python.
# See ``{output_path}/gwas/replace_with_new.sh`` for this script
# And runs some statistical gwas tests for visualization with R
# Pay attention to the gtf script. It does quite a bit to get the trend test usable for R

import logging
import yaml
import pandas as pd
import janitor


logging.basicConfig(level=logging.INFO, filename = '_plinkpca.log', filemode = 'w')

#Load config settings
with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

code_path = config["CODE_PATH"]
code_path = f'{code_path}/popgen/scripts/gwas'

data_path=config["DATA_PATH"]


input_path = f'{data_path}/molino2/popgen/gwas'

vcf_input_file = f'{input_path}/nomolino.Amex3.0_surface.vcf.gz'
gwas_input_file = f'{input_path}/gwas.Amex3.0_surface'

gwas_output_path = f'{data_path}/molino2/popgen/gwas'



#Handy function to pass bash commands to and generate batch scripts
def single_script(cmd_generator):
    sha_bang="#!/usr/bin/bash"

    prefix, command = cmd_generator
    
    batch_info=f'\n#SBATCH --job-name={prefix} \
        \n#SBATCH --cpus-per-task=4 \
        \n#SBATCH --time="24-00:00" \
        \n#SBATCH --mem=32G \
        \n#SBATCH --mail-user=rk2643@stowers.org \
        \n#SBATCH --mail-type=FAIL \
        \n#SBATCH --mail-type=END\
        \n#SBATCH --output=./slurmout/{prefix}.%A.out \
        \n#SBATCH --error=./slurmout/{prefix}.%A.err'

    file_name = f'{prefix}.sh'
    with open(f'{code_path}/{file_name}', 'w') as file:
        file.write(f'{sha_bang}\n{batch_info}\n\n{command}')
    print(f'making {file_name} script')
    print(f'script location: {code_path}/{file_name}')
    return(file_name)



######
#fam file format: A text file with no header line, and one line per sample with the following six fields:
#Family ID ('FID')
#Within-family ID ('IID'; cannot be '0') << sample as sample id apparently
#Within-family ID of father ('0' if father isn't in dataset)
#Within-family ID of mother ('0' if mother isn't in dataset)
#Sex code ('1' = male, '2' = female, '0' = unknown)
#Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

#phenotype = control/unaffected (1) or affected (2)
#in this case unaffected = eyed, surface, affected = eyeless
fam_pheno_key = {'C':['C','2'], 'E':['E','1'], 'S':['S','1']}

def fix_fam_file():
    #ped_df = ped_recode(genome)
    print('fixing fam files')
    fam_input = f'{gwas_input_file}.fam'
    print(fam_input)
    #fam_input= f'{output_path}/gwas/{genome}.pruned.fam'
    fam_in_df = pd.read_csv(fam_input, header = None, sep=' ')

    #replace sample names in "FamilyID" with actual FamilyID; pheno with actual pheno
    fam_replacement = {}
    pheno_replace = {}
    for index,row in fam_in_df.iterrows():
        for key in fam_pheno_key:
            if row[0].startswith(key):
                fam_replacement[row[0]]=key
                pheno_replace[row[0]]=fam_pheno_key[key][1]
    pheno_df = pd.DataFrame.from_dict(pheno_replace, orient='index')
    #print(fam_replacement)

    fam_in_df[5]= pheno_df[0].values
    fam_in_df[0]=fam_in_df[0].map(fam_replacement).fillna(fam_in_df[0])
    print(fam_in_df)

    fam_in_df.to_csv(fam_input, index=False, header=False, sep=' ')


#this makes ped and map files
#and filters out missing alleles
def make_ped_files():
   
    ped_output=f'{gwas_output_path}/gwas.Amex3.0_surface'
    #add in ped/map files, also filters out alleles with missing alleles (maf)
    plink_topedmap_cmd = f'plink \
        --vcf {vcf_input_file}\
        --out {ped_output}\
        --maf 0.05\
        --allow-extra-chr\
        --chr 1-25\
        --recode\
        --set-missing-var-ids @:#$1:$2'

    ped_script_cmd = ['makeped', f'{plink_topedmap_cmd}']
    single_script(ped_script_cmd)
    return(ped_output)


#This takes the existing .ped file and generates a new .ped file with the correct family id, sample id, and phenotype data
#Then use a custom awk script in the output path to replace whole columns with the new .ped values here
#Have to use the awk because the full size files are bananas big
def ped_recode():
    ped_input = f'{gwas_output_path}/gwas.Amex3.0_surface.ped'

    col_names = ['FamilyID', 'SampleID', 'Phenotype']
    #edit id's to match fam_pheno_key
    values_to_edit = pd.read_csv(ped_input, sep=' ', usecols= [0,1,5], names=col_names,  header=None)    
    new_values_dict={}
    new_index_list = []
    for index,row in values_to_edit.iterrows():
        fam_id = row["FamilyID"]
        samp_id = row["SampleID"]

        new_index_list.append(fam_id + '_' + samp_id)
    
        real_fam_id = fam_pheno_key[fam_id[0]][0]
        real_pheno_id = fam_pheno_key[fam_id[0]][1] 
        
        real_samp_id = ''
        if samp_id.startswith("Rnd"):
            real_samp_id = f'{fam_id}_{samp_id}'
            #print(index,fam_id,samp_id,real_samp_id)
        else:
            real_samp_id = samp_id

        new_values_dict[index]=[real_fam_id, real_samp_id, real_pheno_id]
    
    new_ped_df = pd.DataFrame.from_dict(new_values_dict, orient = "index")
    new_ped_df.columns = col_names
    new_ped_df.index = new_index_list

    print(values_to_edit, new_ped_df)
    new_ped_df.to_csv(f"{gwas_output_path}/newids.gwas.Amex3.0_surface.ped", index=False, header=False, sep=' ')

    #use a custom awk to replace each column. see `{output_path}/gwas/replace_with_new.sh`
    logging.info(f'old: {values_to_edit}, new: {new_ped_df}')


#tests done:
    # --assoc w/ ci
    # --model to make cochrane-armitage
    # --fst
    # --logistic
    # --adjust does fdr adjustment, Kyle did this as a supp. 


#in "adjust"- FDR adjusted, use FDR-BH (Benjamini-Hoschberg)
def test_assoc():
    #--assoc = basic association
    output_file = f'{gwas_output_path}/tests/assoc_test'

    assoc_cmd = f'echo "association test"\n\
        plink \
        --file {gwas_input_file}\
        --assoc \
        --adjust \
        --allow-extra-chr\
        --allow-no-sex\
        --out {output_file}'

    script_cmd = ['assoc', assoc_cmd]
    file_name = single_script(script_cmd)    
    return(file_name)

def trend_test():
    test_name = 'CA_trend'

    output_file = f'{gwas_output_path}/tests/{test_name}'

    assoc_cmd = f'echo "{test_name}"\n\
        plink \
        --file {gwas_input_file}\
        --model  \
        --allow-extra-chr\
        --allow-no-sex\
        --out {output_file}'

    script_cmd = [test_name, assoc_cmd]
    file_name = single_script(script_cmd)    
    return(f'{output_file}.model')


##### Get gene names from chr:bp > GTF
#add bp, filter to trend
#the output has some duplicated stuff, filter to distict in r
def fix_gtf():
    snp_anno_path = f'{gwas_output_path}/snp_annotation'
    """
    #pull in trend file
    trend_input_file = trend_test()
    trend_col_names = ['CHR', 'SNP', 'A1', 'A2', 'TEST', 'AFF', 'UNAFF', 'CHISQ', 'DF', 'P']
    trend_df = pd.read_csv(trend_input_file, delim_whitespace=True, header=0, names=trend_col_names)
    #skiprows= 67169300 #this is for testing because there's a lot

    #original snp format is 'chr:bp:'
    #fix the snp to be just bp
    #in plotting with R, just set 'SNP'=gene; 'BP'=SNP
    new_snp = trend_df['SNP'].str.extract(r'(.*):(.*):')
    trend_df['SNP']=new_snp[1]
    trend_df[['CHR','SNP']] = trend_df[['CHR','SNP']].apply(pd.to_numeric)

    #apply filters, filter to trend only, reduce to low pvalues
    #reduces to 10,092 sites
    trend_df = trend_df.loc[trend_df['TEST']=='TREND']
    trend_df = trend_df.loc[trend_df['P']!='NaN']
    trend_df = trend_df.loc[trend_df['P']<1e-6]
    print(f'trend_df post filter = {trend_df.head()}; {trend_df.shape}; {trend_df.dtypes}')
    trend_df.to_csv(f'{gwas_output_path}/tests/less_p_filtering_gwas/pe6_filt_trend.txt', sep ='\t', index=False)
    """
    trend_input_file = f'/n/projects/rk2643/caballo_moro_genomics/streamlined/reports/GWAS_pfiltered_output.tsv'
    trend_col_names = ['CHR', 'SNP', 'A1', 'A2', 'TEST', 'AFF', 'UNAFF', 'CHISQ', 'DF', 'P']
    trend_df = pd.read_csv(trend_input_file, delim_whitespace=True, header=0, names=trend_col_names)
    print(f'trend_df: {trend_df}')

    #pull in refseq names, start, stop and gene name from .gtf file (done with awk)
    chr_gene_names = ['real','start', 'end', 'gene']
    chr_gene_df=pd.read_csv(f'{snp_anno_path}/chr_bprange_gene.csv', names=chr_gene_names, header=None)
    #print(chr_gene_df.head())
    print(f'Imported chr:gene: {chr_gene_df.head()}')
    #print(f'head of chr-gene file:\n{chr_gene_df}')

    #this is a table with refseq chr name and common chr name
    rename_key_names = ['real', 'simple_ish']
    rename_key_df=pd.read_csv(f'{snp_anno_path}/rename_chr.txt', names=rename_key_names, header=None, sep=' ')
    print(rename_key_df.head())



    #map chr with gene + to common chr name, reduce to the 25 chr names
    #conditional join needs the values to be numeric because bp are matched to regions from gtf

    gene_with_simple_one = pd.merge(chr_gene_df, rename_key_df, on = 'real')
    gene_with_simple_df = gene_with_simple_one[['simple_ish','start', 'end', 'gene']]
    gene_with_simple_df.iloc[:,1:3].apply(pd.to_numeric)
    gene_with_simple_df=gene_with_simple_df[gene_with_simple_df['simple_ish'].between(1,25, inclusive='both')]


    print(f'added SNP to gene-chr db: {gene_with_simple_df.head()}')
    #print(f'head of snp,gene, chr:\n {gene_with_simple_df}')
    
    """
    for index, row in trend_snps_only.iterrows():
        if row['CHR'] == gene_with_simple_df['simple_ish']:
            print(gene_with_simple_df['simple_ish'])
            if row['SNP'].between(gene_with_simple_df['start'], gene_with_simple_df['end']):
                print(row, gene_with_simple_df['start'], gene_with_simple_df['end'], gene_with_simple_df['gene'])
            break
        break
    """

    trendsnps_to_gene_df = trend_df.conditional_join(gene_with_simple_df, 
                                          ('SNP', 'start', '>='),
                                          ('SNP', 'end', '<='),
                                          ('CHR', 'simple_ish', '=='),
                                          how='left')
    #print(f'final snp to gene file: {trendsnps_to_gene_df.keys()}')
    
    #print(f'head of db file: {trendsnps_to_gene_df}')
    #print(f'select columns: {trendsnps_to_gene_df[['CHR','SNP','P', 'start', 'end', 'gene']]}')
    #print(trendsnps_to_gene_df.loc[(trendsnps_to_gene_df['CHR']=='5')&(trendsnps_to_gene_df['SNP']=='45464322')])
    

    """pyjanitor:
    A.conditional_join(
       B (range), 
       ('A_value', 'B_low', '>='), 
       ('A_value', 'B_high', '<='), 
       how = 'right')"""

    print(trendsnps_to_gene_df)
    print(f'{gwas_output_path}/tests/GWAS_out.txt')
    trendsnps_to_gene_df.to_csv(f'{gwas_output_path}/tests/GWAS_out.txt', sep ='\t', index=False)


def fst_test():
    test_name = 'assoc_fst'

    output_file = f'{gwas_output_path}/tests/{test_name}'

    assoc_cmd = f'echo "{test_name}"\n\
        plink \
        --file {gwas_input_file}\
        --fst  \
        --family\
        --allow-extra-chr\
        --allow-no-sex\
        --out {output_file}'

    script_cmd = [test_name, assoc_cmd]
    file_name = single_script(script_cmd)    
    return(file_name)

#logistical regression test is only for covariates, ignore!
def logistic_test():
    test_name = 'logistic_covar'

    output_file = f'{gwas_output_path}/tests/{test_name}'

    assoc_cmd = f'echo "{test_name}"\n\
        plink \
        --file {gwas_input_file}\
        --fst  \
        --family\
        --ci 0.9999\
        --allow-extra-chr\
        --allow-no-sex\
        --out {output_file}'

    script_cmd = [test_name, assoc_cmd]
    file_name = single_script(script_cmd)    
    return(file_name)


def chr21_gtf():
    snp_anno_path = f'{gwas_output_path}/snp_annotation'

    #pull in trend file
    trend_input_file = f'{gwas_output_path}/tests/only_chr21/chr21_annot_only.txt'
    trend_col_names = ['CHR', 'SNP', 'A1', 'A2', 'TEST', 'AFF', 'UNAFF', 'CHISQ', 'DF', 'P', 'unk','start','stop','gene']
    trend_df = pd.read_csv(trend_input_file, sep = '\t', header=0, names=trend_col_names)
    print(trend_df.head())
    #skiprows= 67169300 #this is for testing because there's a lot

    #apply filters, filter to trend only, reduce to low pvalues
    #reduces to 10,092 sites

    trend_df = trend_df.loc[trend_df['P']<1e-8]
    trend_df[['CHR','SNP']] = trend_df[['CHR','SNP']].astype('int64')
    print(f'trend_df post filter = {trend_df.head()}')
    trend_df.to_csv(f'{gwas_output_path}/tests/only_chr21/ch21_annot_only.txt', sep ='\t', index=False)


    #pull in refseq names, start, stop and gene name from .gtf file (done with awk)
    chr_gene_names = ['real','position','start', 'end', 'gene_name']
    chr_gene_df=pd.read_csv(f'{snp_anno_path}/use_this_chr21.csv', delim_whitespace=True,names=chr_gene_names, header=0)
    print(chr_gene_df.head())
    chr_gene_df[['real','start','end']] = chr_gene_df[['real','start','end']].fillna(0.0).astype(int, errors='ignore')
    print(type(chr_gene_df['start'][1]))
    print(f'Imported chr:gene: {chr_gene_df.shape}')

    trendsnps_to_gene_df = trend_df.conditional_join(chr_gene_df, 
                                          ('SNP', 'start', '>='),
                                          ('SNP', 'end', '<='),
                                          ('gene', 'gene_name', '=='),
                                          how='left')
    #print(f'final snp to gene file: {trendsnps_to_gene_df.keys()}')
    
    #print(f'head of db file: {trendsnps_to_gene_df}')
    #print(f'select columns: {trendsnps_to_gene_df[['CHR','SNP','P', 'start', 'end', 'gene']]}')
    #print(trendsnps_to_gene_df.loc[(trendsnps_to_gene_df['CHR']=='5')&(trendsnps_to_gene_df['SNP']=='45464322')])
    

    """pyjanitor:
    A.conditional_join(
       B (range), 
       ('A_value', 'B_low', '>='), 
       ('A_value', 'B_high', '<='), 
       how = 'right')"""

    trendsnps_to_gene_df.to_csv(f'{gwas_output_path}/tests/only_chr21/filtered_annot_trend.txt', sep ='\t', index=False)

def go_filter_to_gene():
    snp_anno_path = f'{gwas_output_path}/snp_annotation'

    #pull in trend file
    trend_input_file = f'{gwas_output_path}/tests/goterm/final_GO_list.txt'
    trend_col_names = ['CHR', 'SNP', 'P', 'gene']
    trend_df = pd.read_csv(trend_input_file, sep='\s+', header=0, names=trend_col_names)
    print(trend_df.head())
    #skiprows= 67169300 #this is for testing because there's a lot

    #pull in refseq names, start, stop and gene name from .gtf file (done with awk)
    chr_gene_names = ['real','start', 'end', 'gene_name']
    chr_gene_df=pd.read_csv(f'{snp_anno_path}/chr_bprange_gene.csv',names=chr_gene_names, header=0)
    print(chr_gene_df.head())
    chr_gene_df[['real','start','end']] = chr_gene_df[['real','start','end']].fillna(0.0).astype(int, errors='ignore')
    print(type(chr_gene_df['start'][1]))
    print(f'Imported chr:gene: {chr_gene_df.shape}')

    trendsnps_to_gene_df = trend_df.conditional_join(chr_gene_df, 
                                          ('SNP', 'start', '>='),
                                          ('SNP', 'end', '<='),
                                          ('gene', 'gene_name', '=='),
                                          how='left')
    #print(f'final snp to gene file: {trendsnps_to_gene_df.keys()}')
    
    #print(f'head of db file: {trendsnps_to_gene_df}')
    #print(f'select columns: {trendsnps_to_gene_df[['CHR','SNP','P', 'start', 'end', 'gene']]}')
    #print(trendsnps_to_gene_df.loc[(trendsnps_to_gene_df['CHR']=='5')&(trendsnps_to_gene_df['SNP']=='45464322')])
    

    """pyjanitor:
    A.conditional_join(
       B (range), 
       ('A_value', 'B_low', '>='), 
       ('A_value', 'B_high', '<='), 
       how = 'right')"""

    trendsnps_to_gene_df.to_csv(f'{gwas_output_path}/tests/goterm/go_filtered_with_anot.tsv', sep ='\t', index=False)



if __name__ == '__main__':
    #print('fixing fam file')
    #fix_fam_file()

    #print('making ped and map files')
    #make_ped_files()

    #print('fixing ped file')
    #ped_recode()

    #really the fix below is too too big for memory to handle
    #filter firrrsssttt then annotate
    fix_gtf()
    #print('association test scripts')
    #test_cmd = f'sh {test_assoc()}&&\n sh {trend_test()} &&\n sh {fst_test()}'
    #single_script(['all_tests', test_cmd])

    #chr21_gtf()
    #go_filter_to_gene()