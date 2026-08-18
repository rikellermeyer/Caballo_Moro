#!/users/9/kell3262/miniforge3/envs/CMpop/bin/python

## Note, between CM variant calling/hard_filtering and here, we are using different conda env

## this does some prep work for vcf files for downstream PCA/GWAS analysis
# First, rename chromosomes in original variant vcf file 
# Then, either remove Molino for GWAS or keep it in for admixture/PCA
# Remove all chromosomes outside of 1-25
# For pca, pair down to 2 nocalls and does LD filtering
# Makes pca files and admixture files, as well as the cross-validation error graph for admixture
# NOTE: Admixture files always output to the script location & have to be moved.
# For gwas, remove all molinos and does LD filtering, then run `gwas_fixandtest.py`

import logging
import yaml
import glob
import os
import pandas as pd
import subprocess as sbp
import matplotlib.pyplot as plt


logging.basicConfig(level=logging.INFO, filename = '_plinkpca.log', filemode = 'w')

#Load config settings
with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome = config["GENOME"]
genome_path = config["GENOME_PATH"]
molino_sample_names=config["MO_ONLY_SAMPLES"]

original_vcf_path = config["CM_VCF_PATH"]

#This comes from the hard filter step
#original_vcf = f'{original_vcf_path}/{genome}.all.filtered.vcf.gz'
original_vcf = f'{original_vcf_path}/2025_Amex3.0_surface.variants.filtered.vcf.gz'

#scripts go here
code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/popgen/scripts'

report_path_prefix = config["REPORT_PATH"]

#outputs go here
data_path = config["DATA_PATH"]
common_output_path = f'{data_path}/popgen'

#ALL outputs, vcf included go to popgen, not variants
preplink_output_path = f'{common_output_path}/preplink'
pca_output_path = f'{common_output_path}/pca'
gwas_output_path = f'{common_output_path}/gwas'


#Handy function to pass bash commands to and generate batch scripts
def generate_scripts(cmd_generator):
    sha_bang="#!/usr/bin/bash"

    prefix, command, n = cmd_generator
    
    batch_info=f'\n#SBATCH --job-name={prefix} \
        \n#SBATCH --array=0-{n-1}%{n} \
        \n#SBATCH --cpus-per-task=4 \
        \n#SBATCH --time="2-00:00" \
        \n#SBATCH --mem=32G \
        \n#SBATCH --mail-user=kell3262@umn.edu \
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

# For when you don't have an array of scripts to run, but still want to use sbatch
def single_script(cmd_generator):
    sha_bang="#!/usr/bin/bash"

    prefix, command = cmd_generator

    file_name = f'{prefix}.sh'
    
    with open(f'{code_path}/{file_name}', 'w') as file:
        file.write(f'{sha_bang}\n\n{command}')
    print(f'making {file_name} script')
    print(f'script location: {code_path}/{file_name}')
    return(file_name)



# This is some preprocessing for the vcf files that is common to both pca and gwas
# First, rename chromosomes in original variant vcf file {genome}.biallelic.filtered.vcf.gz
# Then, remove M2 from analysis
# Remove all chromosomes outside of 1-25 (this is a plink problem because it's made for human chromosome names)
def common_preplink():
    n=1

    ### rename chromosomes
    input_file = original_vcf
    input_code = 'renamed_chr'
    #basically a key of refseqid:chr
    new_file_name= f'{original_vcf_path}/{genome}.renamechr.txt'

    renamed_output_file = f'{common_output_path}/preplink/{genome}.{input_code}.vcf.gz'

    #rename.txt = "old_name new_name\n"
    accession_to_chr = f'{genome_path}/chr_key.txt'
    col_names = ['Chromosome name','RefSeq seq accession']
    acc_to_chr_df = pd.read_csv(accession_to_chr, names = col_names, sep='\t')


    #need to swap the columns
    column_list = list(acc_to_chr_df.columns)
    x, y = column_list.index('Chromosome name'), column_list.index('RefSeq seq accession')
    column_list[y], column_list[x] = column_list[x], column_list[y]
    chr_acc_key_df = acc_to_chr_df[column_list]

    #gives everyone their common chromosome names    
    chr_acc_key_df['Chromosome name'] = chr_acc_key_df.index + 1
    #print(chr_acc_key_df)
    def replace_func(x):
        if x >= 110:
            return str(x)
        else:
            return x

    chr_acc_key_df['Chromosome name'] = chr_acc_key_df['Chromosome name'].apply(replace_func)

    logging.info(f'{chr_acc_key_df}')

    #ouput the "key" for chromosome to accession for bcftools annotate
    chr_acc_key_df.to_csv(new_file_name, index=False, sep=' ', na_rep='NA')
    
    rename_cmd = f'bcftools annotate \
        --rename-chrs {new_file_name} \
        --write-index\
        -o {renamed_output_file} {input_file}'
    #print(rename_cmd)
    
    script_cmd = ['rename_vcf', rename_cmd, n]
    generate_scripts(script_cmd)


    #pair down to 25 chr
    #input is the renamed_output_file from above
    chr_corrected_output = f'{preplink_output_path}/chr_corrected.{genome}.vcf.gz'

    correct_vcf_cmd = f'bcftools view \
        -o {chr_corrected_output} \
        --force-samples {renamed_output_file}\
        --write-index\
        --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25'
    
    final_command =f'echo "correcting vcf"\n{correct_vcf_cmd}'
    script_cmd= ['preplink', final_command, n]

    generate_scripts(script_cmd)
    return(chr_corrected_output)


def pca_plink():
    corrected_vcf_input = common_preplink()
    #make gatk index file, needed for SelectVariants
    index_cmd = f'gatk IndexFeatureFile \
        -I {corrected_vcf_input}'

    #filter to 2 nocall, no M2
    nocall2_output = f'{pca_output_path}/nocall2.{genome}.vcf.gz'
    nocall_filt_cmd = f'gatk SelectVariants\
        -V {corrected_vcf_input}\
        --max-nocall-number 2\
        --exclude-sample-name M2\
        -O {nocall2_output}'

    window_sizes = ["50", "100", "150", "300"]
    n=0
    # LD filtering and SNP naming
    for window in window_sizes:
        LD_output=f'{pca_output_path}/LD_filter_in.{window}win.{genome}.pruned.{window}win'
        plink_LD_filter = f'plink --indep-pairwise {window} 10 0.1 \
            --vcf {nocall2_output} \
            --out {LD_output}\
            --double-id\
            --allow-extra-chr\
            --set-missing-var-ids @:#'
    
        pruned_output=f'{pca_output_path}/pruned.{genome}.{window}.win'
        plink_tobed_cmd = f'plink --pca\
            --vcf {nocall2_output}\
            --out {pruned_output}\
            --extract {LD_output}.prune.in\
            --make-bed\
            --double-id\
            --allow-extra-chr\
            --set-missing-var-ids @:#'
        
        #keep an eye on admixture output path...
        # aka the logs will output where you want 
        # but the admixture files ##### ALWAYS output to the script directory #########
        admix_output=f'{pca_output_path}/admixture/admix.{genome}.{window}.win'
        admix_cmd=f"""for i in {{1..5}}; \n\
            do \n\
            ./admixture -j16 --cv {pruned_output}.bed $i > {admix_output}.log${{i}}.out \n\
            done"""

        by_window_cmd = f'echo "generating LD filter"\n\
            {plink_LD_filter}&&\n\
            echo "making bed files"\n\
            {plink_tobed_cmd}&&\n\
            echo "running admixture"\n\
            {admix_cmd}'
        single_script_cmd =[f'{window}_fullpca_{n}', by_window_cmd]
        single_script(single_script_cmd)
        n+=1

    #print('making batch pca files')
        
    pre_window_cmd = f'echo "Indexing vcf for gatk"\n\
        {index_cmd}&&\n\
        echo "filtering to nocalls"\n\
        {nocall_filt_cmd}'
    pre_window_script = ['pre_window_pca', pre_window_cmd]
    pre_window_file = single_script(pre_window_script)

    batch_window_cmd = f'sh *_fullpca_${{SLURM_ARRAY_TASK_ID}}.{genome}.sh'
    final_cmd=f'echo "preprocessing pca"\n \
        sh {pre_window_file}&&\n \
        echo "by-window plink commands"\n \
        {batch_window_cmd}'
    #print(batch_window_cmd, final_cmd)

    generate_scripts([f"pca_plink", final_cmd, n])

#Plots the cross-validation error for each K in admixture run
def plot_cv():
    cv_files = glob.glob(f'{pca_output_path}/admixture/*.out')
    #print(cv_files)

    cv_dict = {'window':[], 'k_value':[], 'cv':[]}
    for file in cv_files:
        with open(file, 'r') as read_file:
            for line in read_file:
                if line.startswith('CV'):
                    #print(line)
                    cv_value = float(line.rsplit(r': ')[1].rstrip())
                    cv_dict['cv'].append(cv_value)
                    window = int(os.path.basename(file).split(sep='.')[3])
                    cv_dict['window'].append(window)
                    k_value = os.path.basename(file).split(sep='.')[5][3]
                    cv_dict['k_value'].append(k_value)
    print(cv_dict)
    cv_df = pd.DataFrame(cv_dict).sort_values(by=['k_value', 'window'])
    #df.Col = df.Col.astype(float)

    print(cv_df)
    #df.sort_values(by=['col1', 'col2'])

    window_set = set(cv_df['window'])
    plt.figure()
    for window in window_set:
        selected_data = cv_df.loc[cv_df['window'] ==window]
        plt.plot(selected_data['k_value'], selected_data['cv'], label = window)

    plt.legend()
    plt.xlabel('K-value')
    plt.ylabel('Cross-validation error')
    #plt.axis([None, None, 0.40, 0.55])
    #plt.ylim(top=0.55)
    #axis = plt.gca()
    #axis.set_ylim([0.4,0.55])

    plt.show()
    plt.savefig(f'{report_path_prefix}/CMC_graphics_files/admixture_cross_validation_error', \
                dpi='figure', format=None)
    #plt.plot(k_list, cv_list)
    #plt.show()
    #plt.savefig('test', dpi='figure', format=None)


                        
# Makes plink files and preps them for GWAS tests
def gwas_plink():
    #remove all molinos
    corrected_vcf_input = common_preplink()

    #print(type(molino_sample_names))
    #ignore_molino = ','.join(molino_sample_names)
    ignore_molino = molino_sample_names.replace(" ", ",")
    print("Removing:", ignore_molino)

    no_molino_output = f'{gwas_output_path}/nomolino.{genome}.vcf.gz'
    remove_molino_cmd = f'bcftools view -s^{ignore_molino} \
        -o {no_molino_output} \
        --write-index\
        --force-samples {corrected_vcf_input}'
    
        
    #filter for LD
    LD_output = f'{gwas_output_path}/LD_filter.{genome}'
    LD_filter_cmd = f'plink \
        --vcf {no_molino_output}\
        --out {LD_output}\
        --double-id\
        --chr-set 25\
        --keep-allele-order\
        --indep-pairwise 50 10 0.1'

    #convert to bed/binaries, extractLD
    pruned_output = f'{gwas_output_path}/gwas.{genome}'
    plink_tobed_cmd = f'plink \
        --vcf {no_molino_output} \
        --out {pruned_output}\
        --extract {LD_output}.prune.in\
        --make-bed\
        --double-id\
        --chr-set 25\
        --keep-allele-order'
    
    final_cmd = f'echo "removing molino samples"\n \
        {remove_molino_cmd}&&\n \
        echo "making LD files"\n\
        {LD_filter_cmd}&&\n \
        echo "making bed files"\n \
        {plink_tobed_cmd}'
    n=1
    generate_scripts(['gwas_plink', final_cmd, n])



if __name__=='__main__':

    #pca_plink and gwas_plink will both call on common_preplink(), so only use when testing

    #print('generating common plink scripts')
    #common_preplink()

    #print('making pca scripts')
    #pca_plink()

    #plot_cv()

    print('making gwas scripts')
    gwas_plink()

    

