#!/home/rk2643/miniforge3/envs/gatk/bin/python

## After filtering by genotype/intersection from `genotype_filter.py`
## Annotate the variants with snpEFF
## Continue on with candidate gene selection

### Admin ###
import yaml
import os
import glob

with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


### File locations ###
data_path = config["DATA_PATH"]
code_path = config["CODE_PATH"]
genome = config["GENOME"]
gatk_genome_path = config["GATK_GENOME_PATH"]
gatk_genome_format = f'{gatk_genome_path}/{genome}.fna' #for any gatk
sample_info = config["SAMPLES"]

sha_bang = '#!/usr/bin/bash'
script_path = f'{code_path}/variant_calling/scripts/snpeff'

### SBATCH parameters ###
def sbatch_header(script_name, array_length):
    if array_length >1:
        array_param = f'#SBATCH --array=0-{array_length-1}%{array_length}'
    else:
        array_param=''
    sbatch_params = f'{array_param}\
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

######
def make_snpeff_script(filter_type, nocalls):
    input_vcf_path = f'{data_path}/variant_calling/variants/geno_filter/{filter_type}/{genome}.{nocalls}nocalls/intersected'
    input_files = glob.glob(f'{input_vcf_path}/*.intersected.vcf.gz')
    input_file = ''.join(input_files)

    output_path = f'{input_vcf_path}/snpeff'
    output_file = f'{output_path}/{filter_type}.{nocalls}nocalls.snpeff.vcf'
    summary_file = f'{output_path}/{filter_type}.{nocalls}nocalls.summary.html'
    check_paths([output_path, script_path])

    snpeff_location = 'snpeff-5.2-0/snpEff.jar'
    snpeff_config_path = 'snpeff-5.2-0/snpEff.config'

    snpeff_cmd = f'java -jar {snpeff_location} eff\
        -c {snpeff_config_path} -v\
        -s {summary_file}\
        Astyanax_mexicanus-3.0.103\
        {input_file} > {output_file}'
    
    with open(f'{script_path}/run_snpeff.sh', 'w') as file:
        file.write(f'{sha_bang}\n{sbatch_header("snpeff",1)}\n{snpeff_cmd}')


if __name__ == '__main__':
    filter_type = 'eyeless_dominant'
    nocalls = '0'
    make_snpeff_script(filter_type, nocalls)
