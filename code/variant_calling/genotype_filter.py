#!/home/rk2643/miniforge3/envs/gatk/bin/python

## This makes scripts to filter by number of no calls
## Separates into phenotype (eyed, eyeless, surface)
## And then by genotype (Alt/Ref)
## For example: snps that are only recessive (ALT) in eyeless, ref in eyed/surface

import yaml
import os
import glob

### Admin ###

with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


### File locations ###
data_path = config["DATA_PATH"]
code_path = config["CODE_PATH"]
genome = config["GENOME"]
gatk_genome_path = config["GATK_GENOME_PATH"]
sample_info = config["SAMPLES"]
#could use config["POPULATIONS"]

sha_bang = '#!/usr/bin/bash'
script_path = f'{code_path}/variant_calling/scripts/geno_filter'
gatk_genome_format = f'{gatk_genome_path}/{genome}.fna' #for any gatk


#Filter by genotype with snps only
input_all_vcf = f'{data_path}/variant_calling/variants/snps/filter_done.{genome}.snps.vcf.gz'
output_filtered_vcf = f'{data_path}/variant_calling/variants/geno_filter'

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

def make_file_names(genome, subset_type):
    subset_specific_path = f'{output_filtered_vcf}/{subset_type}'
    check_paths([subset_specific_path])

###################
#general definitions
phenotype_keys = {'eyeless':'C+','eyed':'E+','surface':'S+'}


# split into population/phenotype groups
def split_by_phenotype():

    check_paths([script_path, output_filtered_vcf])
    check_paths([f'{script_path}/slurmout'])

    input_param = f'-V {input_all_vcf}'
    ref_param = f'-R {gatk_genome_format}'

    split_cmd_dict={}
    for pheno, value in phenotype_keys.items():
        output_pheno_vcf=f'-O {output_filtered_vcf}/{genome}.filtered.{pheno}.vcf.gz'
        sample_param = f'--sample-expressions {value}'
        split_cmd_dict[pheno]=[output_pheno_vcf,sample_param]
    n=0

    for key in split_cmd_dict:
        output_param=split_cmd_dict[key][0]
        sample_param=split_cmd_dict[key][1]
        split_cmd = f'gatk SelectVariants\
            {ref_param}\
            {input_param}\
            {output_param}\
            {sample_param}'
        
        with open(f'{script_path}/{genome}.{key}.split.{n}.sh', 'w') as file:
            file.write(f'{sha_bang}\n{split_cmd}')
        n+=1

    with open(f'{script_path}/{genome}.splitbypheno.sh', 'w') as file:
        file.write(f'{sha_bang}\n{sbatch_header('splitpheno', n)}\nsh *.split.${{SLURM_ARRAY_TASK_ID}}.sh')


### Genotype filters
def eyeless_recessive():
    
    surface_key = "-select 'vc.getHomRefCount() == 6'" #fixed/homozygous reference
    eyed_key = "-select 'vc.getHomVarCount() == 0'" #essentially no fixed alt
    eyeless_key = "-select 'vc.getHomRefCount() == 0 && vc.getHetCount() == 0'" #no hets, no fix reference aka hom. alt

    eyeless_recessive_dict={'surface':surface_key,'eyed':eyed_key,'eyeless':eyeless_key}

    return('eyeless_recessive', eyeless_recessive_dict)

def eyeless_dominant():
    #R/R
    surface_key = "-select 'vc.getHomRefCount() == 6'" #fixed/homozygous reference
    #R/R
    eyed_key = "-select 'vc.getHomRefCount() == 6'" #no alt
    #A/R or A/A
    eyeless_key = "-select 'vc.getHomRefCount() == 0'" #no fixed reference

    eyeless_dominant_dict={'surface':surface_key,'eyed':eyed_key,'eyeless':eyeless_key}

    return('eyeless_dominant', eyeless_dominant_dict)

#take the genotype-filtering-parameters and select the appropriate variants
#The parameters are different for each phenotype, so `geno_type` calls on the appropriate one
#i.e. eyeless "recessive" (homozygous ALT) or eyeless "dominant" (heterozygous ALT) vs. surface/eyed reference
def geno_split_script(geno_type, nocalls):
    geno_name, genotype_dictionary = geno_type

    output_vcf_path = f'{output_filtered_vcf}/{geno_name}/{genome}.{nocalls}nocalls'
    check_paths([output_vcf_path])
    
    n=0
    for pheno in phenotype_keys:
        input_vcf = f'{output_filtered_vcf}/{genome}.filtered.{pheno}.vcf.gz'
        output_vcf = f'{output_vcf_path}/{genome}.{nocalls}nocalls.{geno_name}.{pheno}.vcf.gz'
        select_param = genotype_dictionary[pheno]

        filter_cmd = f'gatk SelectVariants \
            -V {input_vcf}\
            --exclude-filtered\
            --max-nocall-number {nocalls}\
            {select_param}\
            -O {output_vcf}'
        print(filter_cmd)

        with open(f'{script_path}/{geno_name}_filter_{pheno}.{n}.sh', 'w') as file:
            file.write(f'{sha_bang}\n{filter_cmd}')
        
        n+=1

    script_name = f'{geno_name}'
    file_name = f'batch.{script_name}.sh'

    with open(f'{script_path}/{file_name}', 'w') as file:
        file.write(f'{sha_bang}\n{sbatch_header(script_name,n)}\nsh {geno_name}_filter_*${{SLURM_ARRAY_TASK_ID}}.sh')

    return(file_name)


# intersect/test for variants (filtered by genotype condition) that are common between phenotypes
# The intersect identifies the variants that match other variants.
# The output .vcf display the sites in a file that match the criteria (n=3) must be present in all files
def intersect(geno_type, nocalls):
    geno_name, genotype_dictionary = geno_type

    input_vcfs = glob.glob(f'{output_filtered_vcf}/{geno_name}/{genome}.{nocalls}nocalls/*.vcf.gz')
    print(input_vcfs)
    
    output_path = f'{output_filtered_vcf}/{geno_name}/{genome}.{nocalls}nocalls/intersected'
    inputs=' '.join(input_vcfs)

    intersect_cmd = f'bcftools isec -n=3 -O z -p {output_path} {inputs}'

    file_name='interesect.sh'

    with open(f'{script_path}/{file_name}', 'w') as file:
        file.write(f'{sha_bang}\n{sbatch_header('intersect',2)}\n{intersect_cmd}')

    return(file_name)

#bcftools isec produces a list of sites but not in a gatk friendly format
#This converts the chromosome:position into a gatk format
#And then pulls out the variants from the original (all, biallelic) vcf file
def extract_intervals(geno_type, nocalls):
    geno_name, genotype_dictionary = geno_type

    input_sites = f'{output_filtered_vcf}/{geno_name}/{genome}.{nocalls}nocalls/intersected/sites.txt'
    output_intervals = f'{output_filtered_vcf}/{geno_name}/{genome}.{nocalls}nocalls/intersected/sites.intervals'
    vcf_output= f'{output_filtered_vcf}/{geno_name}/{genome}.{nocalls}nocalls/intersected/{geno_name}.{nocalls}nocalls.intersected.vcf.gz'

    extract_cmd = f"cut -f1-2 {input_sites} | \
        sed -E 's/(.+)\t(.+)/\\1:\\2/' > {output_intervals}&&\
        \n gatk SelectVariants\
        -R {gatk_genome_format}\
        -V {input_all_vcf}\
        -O {vcf_output}\
        -L {output_intervals}"
    
    file_name = 'extract_intersect_intervals.sh'
    with open(f'{script_path}/{file_name}', 'w') as file:
        file.write(f'{sha_bang}\n{sbatch_header('extract',1)}\n{extract_cmd}')
    
    return(file_name)


if __name__ == '__main__':
    #split_by_phenotype()
    genotype_filter = eyeless_recessive()
    nocalls = '0'
    #split = geno_split_script(genotype_filter, nocalls)
    intersect_files = intersect(genotype_filter, nocalls)
    #extract = extract_intervals(genotype_filter, nocalls)
