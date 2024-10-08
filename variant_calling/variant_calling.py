#!/home/rk2643/miniforge3/envs/gatk/bin/python

###
# This makes a series of scripts (and SBATCH scripts) to variant calling
# Previous step: alignment_scripts.py
# Split to mapped/unmapped > haplotype caller > combinegvcfs > joint calling 
# Next step: hard_filtering.py
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
logging.basicConfig(level=logging.INFO, filename = '_mkscript.log', filemode = 'w')

with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


### File locations ###
data_path_base = config["DATA_PATH"]
data_path = f'{data_path_base}/molino2'
code_path = config["CODE_PATH"]
genome = config["GENOME"]
gatk_genome_path = config["GATK_GENOME_PATH"]

sha_bang = '#!/usr/bin/bash'
script_path = f'{code_path}/variant_calling/scripts/molino2/variant_calling'
script_errors = f'{script_path}/error_logs' # idk if these are necessary or should go to slurmout/*.err
gatk_genome_format = f'{gatk_genome_path}/{genome}.fna' #for any gatk
sample_info = config["MOLINO_SAMPLES"] # this is a dictionary of seq_id:sample_name


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


#####################################################
##### Labeled by sample name >> HaplotypeCaller #####

### GATK scripts use reference genome formats based on the bwa NOT bwa-mem2
### Split to mapped and unmapped reads ###
##Obviously fix this one day to appropriately match seq_id to sample_name 
def split_to_mapped(sample_name, n):

    input_file = f'{data_path}/alignment/renamed/{sample_name}.marked.bam'
    output_path = f'{data_path}/alignment/mapped'    
    unmapped_error_file = f'{script_errors}/{sample_name}.unmapped.err'

    mapped_error_file = f'{script_errors}/{sample_name}.mapped.err'
    check_paths([output_path, script_path])
    unmapped_out_file = f'{output_path}/{sample_name}.unmapped.bam'

    mapped_out_file = f'{output_path}/{sample_name}.mapped.bam'

    unmapped_cmd = f'samtools view -b -f 4 {input_file} > {unmapped_out_file} 2> {unmapped_error_file}'
    mapped_cmd = f'samtools view -b -F 4 {input_file} > {mapped_out_file} 2> {mapped_error_file}'
    index_cmd = f'samtools index {mapped_out_file}'
    

    script_name = f'{script_path}/{sample_name}_mapping_{n}.sh'
    with open(script_name, 'w') as file:
        file.write(f'{sha_bang}\n{mapped_cmd};\n{unmapped_cmd}\n{index_cmd}')
    return(f'{sample_name}_mapping_{n}.sh')    

#HaplotypeCaller labels REF / ALT variants
def haplotype_caller(sample_name, n):
    input_file = f'{data_path}/alignment/mapped/{sample_name}.mapped.bam'
    temp_space = f'{data_path}/scratch'
    output_path = f'{data_path}/variant_calling/haplotype_caller'
    error_path = f'{script_errors}/{sample_name}.haplo.err'
    check_paths([output_path])

    bam_out_file = f'{output_path}/{sample_name}.bam'
    vcf_out_file = f'{output_path}/{sample_name}.g.vcf.gz'
    assembly_out_file = f'{output_path}/{sample_name}.assembly.region.tsv'

    #index_cmd = f'samtools index {input_file}' 

    haplotype_caller_cmd = f"gatk HaplotypeCaller\
        -R {gatk_genome_format}\
        -I {input_file}\
        -O {vcf_out_file}\
        -ERC GVCF\
        -bamout {bam_out_file}\
        --assembly-region-out {assembly_out_file}\
        --tmp-dir {temp_space}"
    
    script_name = f'{script_path}/{sample_name}_haplotype_{n}.sh'
    with open(script_name, 'w') as file:
        file.write(f'{sha_bang}\n{haplotype_caller_cmd}')
    return(f'{sample_name}_haplotype_{n}.sh')  


#This combines all the gvcfs for variant calling
#input_files_2 are the ones from the original CM samples
def combine_gvcfs(input_file_2):
    input_files_1 = glob.glob(f'{data_path}/variant_calling/haplotype_caller/*.g.vcf.gz')
    int_input_2 = f'{input_file_2}.int'
    all_input_files = input_files_1 + [int_input_2]
    V = ' -V '
    all_the_inputs = f"{V}{f'{V}'.join(all_input_files)}"
    #print(input_file_2, all_input_files, all_the_inputs)

    output_path = f'{data_path}/variant_calling/combined_gvcf'
    output_file = f'{output_path}/all.g.vcf.gz'
    error_file = f'{script_errors}/combine.gvcfs.err'
    check_paths([output_path])
    temp_space = f'{data_path}/scratch'

    combine_cmd = f'gatk CombineGVCFs\
        {all_the_inputs}\
        -R {gatk_genome_format}\
        -O {output_file}\
        --tmp-dir {temp_space} 2> {error_file}'
    
    return(combine_cmd)  

# use the gcvf with all samples to do joint calling
def joint_call():
    input_file = f'{data_path}/variant_calling/combined_gvcf/all.g.vcf.gz'
    output_path = f'{data_path}/variant_calling/joint_call'
    check_paths([output_path])
    output_file = f'{output_path}/jointcall.vcf.gz'
    error_file = f'{script_errors}/joint_call.err'

    

    joint_call_cmd = f'gatk GenotypeGVCFs \
        -R {gatk_genome_format}\
        -V {input_file}\
        -O {output_file}\
        --include-non-variant-sites 2>{error_file}'
    return(joint_call_cmd)

# Batch scripts for splitting mapped/unmapped > joint calling
def batch_mapped_to_joint_calling(sample_names, input_file_2):
    n=0
    #make the split to mapped and haplotype caller scripts
    for sample_name in sample_names:
        map_script = split_to_mapped(sample_name, n)
        haplo_script = haplotype_caller(sample_name, n)

        #make batch scripts for splitting & hapoltype calling
        with open(f'{script_path}/map_haplo_batch_{n}.sh','w') as file:
            file.write(f'{sha_bang}\nsh {map_script} &&\nsh {haplo_script}')
            
        n+=1

    #make the single combine gvcf and joint call scripts
    combine_cmd = combine_gvcfs(input_file_2)
    joint_call_cmd = joint_call()

    header = sbatch_header('hap_comb',n)

    batch_script_1 = f'{sha_bang}\n{header}\n\nsh map_haplo_batch_${{SLURM_ARRAY_TASK_ID}}.sh'

    batch_script_2 = f'{sha_bang}\n{header}\n\n{combine_cmd} && \n{joint_call_cmd}'
    
    with open(f'{script_path}/big_vc_batch_1.sh','w') as file:
        file.write(batch_script_1)
    
    with open(f'{script_path}/big_vc_batch_2.sh','w') as file:
        file.write(batch_script_2)

if __name__ == '__main__':
    print('Making variant calling scripts')

    ## For adhoc molino # 2 addition, combine from original run + first adhoc molino.
    input_file_2 = glob.glob(f'{data_path_base}/variant_calling/joint_call/jointcall.vcf.gz')[0]
    print(input_file_2)
    sample_names = list(sample_info.values())
    print(sample_names)
    batch_mapped_to_joint_calling(sample_names, input_file_2)

