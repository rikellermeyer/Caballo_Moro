#!/home/rk2643/miniforge3/envs/gatk/bin/python

###
# This makes a series of scripts (and SBATCH scripts) to variant calling
# Previous step: variant_calling.py
# Subset into monomorphic/SNPs/INDELs > Label filter parameters > apply filter by parameters
# Subset = SelectVariants ; Label = VariantFiltration ; Apply filter = SelectVariants
# Next step: `genotype_filter.py` or PCA/Admixture
# Should be generalizable by editing the config.yaml file
### 

import yaml
import logging
import os
import glob
import sys
import re
from pathlib import Path


### Admin ###

with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


### File locations ###
data_path_base = config["DATA_PATH"]
data_path = f'{data_path_base}/molino2'
code_path = config["CODE_PATH"]
genome = config["GENOME"]
gatk_genome_path = config["GATK_GENOME_PATH"]
sample_info = config["MOLINO_SAMPLES"]

sha_bang = '#!/usr/bin/bash'
script_path = f'{code_path}/variant_calling/scripts/molino2/hard_filtering'
script_errors = f'{script_path}/error_logs' # idk if these are necessary or should go to slurmout/*.err
gatk_genome_format = f'{gatk_genome_path}/{genome}.fna' #for any gatk

input_all_vcf = f'{data_path}/variant_calling/joint_call/jointcall.vcf.gz'
output_filtered_vcf = f'{data_path}/variant_calling/variants' #subsets > variants/subsets/

### SBATCH parameters ###
def sbatch_header(script_name, array_length):
    sbatch_params = f'#SBATCH --array=0-{array_length-1}%{array_length}\
        \n#SBATCH --job-name={script_name} \
        \n#SBATCH --cpus-per-task=8 \
        \n#SBATCH --mem=32G \
        \n#SBATCH --time="24-00:00" \
        \n#SBATCH --error=./slurmout/{script_name}.%a.%A.err \
        \n#SBATCH --output=./slurmout/{script_name}.%a.%A.out \
        \n#SBATCH --mail-user=rk2643@stowers.org \
        \n#SBATCH --mail-type=FAIL \
        \n#SBATCH --mail-type=END'
    return(sbatch_params)

### Path handling ###
def check_paths(new_path): #always use new_path as a list, even if only one entry
    for file in new_path:
        if os.path.exists(file):
            continue
        else:
            os.makedirs(file, exist_ok=True)
            print(f"Making new paths: {new_path}")

### Dictionary that provides the filtering parameters for each subset ###
subset_params = {'mono':
                            ['--select-type-to-exclude INDEL \
                            --select-type-to-exclude MIXED \
                            --select-type-to-exclude MNP \
                            --select-type-to-exclude SYMBOLIC \
                            --select-type-to-include NO_VARIATION',
                                '-filter "QD < 2.0" --filter-name "QD2" \
                                -filter "FS > 60.0" --filter-name "FS60" \
                                -filter "MQ < 40.0" --filter-name "MQ40" '], 
                    'snps':
                            ['--select-type-to-include MNP\
                              --select-type-to-include SNP', 
                                '-filter "QD < 2.0" --filter-name "QD2" \
                                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                                -filter "SOR > 3.0" --filter-name "SOR3" \
                                -filter "FS > 60.0" --filter-name "FS60" \
                                -filter "MQ < 40.0" --filter-name "MQ40" \
                                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'], 
                    'indels':
                            ['--select-type-to-include MIXED\
                              --select-type-to-include INDEL', 
                                '-filter "QD < 2.0" --filter-name "QD2" \
                                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                                -filter "FS > 200.0" --filter-name "FS200" \
                                -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"']}

def make_file_names(genome, subset_type):
    subset_specific_path = f'{output_filtered_vcf}/{subset_type}'
    check_paths([subset_specific_path, script_path])

    subset_file = f'{subset_specific_path}/subset.{genome}.{subset_type}.vcf.gz'
    labeled_file = f'{subset_specific_path}/labeled.{genome}.{subset_type}.vcf.gz'
    filtered_file = f'{subset_specific_path}/filter_done.{genome}.{subset_type}.vcf.gz'

    return(subset_file, labeled_file, filtered_file)

# First parse out monomorphic (unchanged), snps, and indels
# Then provide the filter ifnormation ("label") by subset-specific "hard filters" aka quality and what not
# Then actually filter out the sites that don't pass the filter
def subset_script(genome):        
    n=0
    for subset_type in subset_params.keys():
        subset_file, labeled_file, filtered_file = make_file_names(genome, subset_type)
        
        print(subset_type, subset_params[subset_type], subset_params[subset_type][0])

        subset_cmd = f'gatk SelectVariants\
            -R {gatk_genome_format}\
            -V {input_all_vcf}\
            -O {subset_file}\
            {subset_params[subset_type][0]}'

        label_cmd = f'gatk VariantFiltration\
            -R {gatk_genome_format}\
            -V {subset_file}\
            -O {labeled_file}\
            {subset_params[subset_type][1]}'

        filter_cmd = f'gatk SelectVariants\
            -R {gatk_genome_format}\
            -V {labeled_file}\
            -O {filtered_file}\
            --exclude-filtered'


        with open(f'{script_path}/filter_{subset_type}_{n}.sh', 'w') as file:
            file.write(f'{sha_bang}\n{subset_cmd}\n{label_cmd}\n{filter_cmd}')

        n+=1

    #ran in 65 min, 1.2G
    with open(f'{script_path}/batch_all_filter.sh', 'w') as file:
        file.write(f'{sha_bang}\n\n{sbatch_header("hard_filter", n)}\nsh filter*_${{SLURM_ARRAY_TASK_ID}}.sh')
        
### Combine back, 1 w/ biallelic, 1 w/ all
def combine_back(genome):
    V = '-I'
    biallelic_inputs = []

    for subset_type in subset_params.keys():
        subset_file, labeled_file, filtered_file = make_file_names(genome, subset_type)
        if subset_type == 'mono':
            mono_input = filtered_file
        else:
            biallelic_inputs.append(f'{V} {filtered_file}')
    print(biallelic_inputs, mono_input)
    
    bi_input_file_list = (' '.join(list(biallelic_inputs)))
    biallelic_output = f'{output_filtered_vcf}/{genome}.biallelic.filtered.vcf.gz'
    biallelic_merge_cmd = f'gatk MergeVcfs\
        -R {gatk_genome_format}\
        {bi_input_file_list}\
        -O {biallelic_output}'

    all_filtered_output = f'{output_filtered_vcf}/{genome}.all.filtered.vcf.gz'
    all_filtered_merge = f'gatk MergeVcfs\
        -R {gatk_genome_format}\
        -I {biallelic_output}\
        -I {mono_input}\
        -O {all_filtered_output}'

    script_name = 'merge_back'
    with open(f"{script_path}/merge_filtered.sh", "w") as file:
        file.write(f'{sha_bang}\n{sbatch_header(script_name,1)}\n{biallelic_merge_cmd}&&\n{all_filtered_merge}')



if __name__=='__main__':
    print('Making filtering scripts')
    molino_sample_info = config["MOLINO_SAMPLES"] # this is a dictionary of seq_id:sample_name
    total_sample_info = ''

    subset_script(genome)
    combine_back(genome)