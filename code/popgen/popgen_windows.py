#!/users/9/kell3262/miniforge3/envs/CMpop/bin/python

#This is a script to set up and run popgenWindows.py from https://github.com/simonhmartin/genomics_general#processing-vcf-files
#popgenWindows from above has to be in the current directory for this script, scripts are output to this directory

import yaml
import pandas as pd
from itertools import combinations
import csv

#Load config settings
with open("../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

data_path = config["DATA_PATH"]
code_path = config["CODE_PATH"]
genome = config["GENOME"]
sample_names = config["SAMPLES"]
population_code = {'C':'eyeless', 'E': 'eyed', 'S': 'surface'}

sha_bang = '#!/usr/bin/bash'
sbatch_params = f'#SBATCH --job-name=popgen \
        \n#SBATCH --cpus-per-task=8 \
        \n#SBATCH --mem=64G \
        \n#SBATCH --time="2-00:00" \
        \n#SBATCH --error=./slurmout/popgen.%a.%A.err \
        \n#SBATCH --output=./slurmout/popgen.%a.%A.out \
        \n#SBATCH --mail-user=kell3262@umn.edu \
        \n#SBATCH --mail-type=FAIL \
        \n#SBATCH --mail-type=END'

#location of CMC all.vcf files
vcf_path = config["CM_VCF_PATH"]

#smartin's genomics scripts
geno_script_path = f'{code_path}/popgen/smartin_genomics'

### STEPS ###
# install, add to path
# generate a population file "sample_name \t population_name \n"
# parse VCF with smartin_genomics/parseVCF.py which takes .all.vcf.gz file and makes a .geno file - see 
# .geno file is a "genotype matrix" specific to the genomics_general code from smartin
# gzip .geno file > .geno.gz
# sbatch script took about 4 hours

def format_setup():

    ## first make a --popsFile 
    # format is sample \t population \n
    print('make Pop file')
    samples_list = sample_names.split(' ')
    samp_dict = population_code
    sample_match_dict = {}
    for sample in samples_list:
        samp_code = sample[0]
        sample_match_dict[sample]=samp_dict[samp_code]
    
    make_into_txt = ''
    for key, value in sample_match_dict.items():
        make_into_txt+=f'{key}\t{value}\n'
    
    with open('samples_to_populations.txt','w') as file:
        file.write(make_into_txt)

    ## convert all.vcf.gz file to .geno files using parseVCF.py
    input_vcf = f'{vcf_path}/2025_{genome}.all.filtered.vcf.gz'
    output_geno = f'{data_path}/popgen/{genome}.all.geno.gz'

    # Index command
    bcf_index_cmd = [f'bcftools index {input_vcf}']
    print(bcf_index_cmd)

    #note that this needs an index (with bcftools index)
    geno_file_cmd = f'{geno_script_path}/VCF_processing/parseVCFs.py -i {input_vcf} --threads 16 \
            | bgzip > {output_geno}' #used to be >2 mkgeno.err, now outputs to slurm
    with open("./scripts/popgenwindows/make_geno_file.sh", "w") as file:
        file.write(f'{sha_bang}\n{sbatch_params}\
                \n{bcf_index_cmd}\n{geno_file_cmd}')

    return output_geno
    

#sbatch took about 4 hours
def run_popgen():
    input_file = format_setup()

    print('Run popgen windows')
    #input_file = f'{data_path}/popgen/{genome}.all.geno.gz'
    output_file = f'{data_path}/popgen/{genome}.50kwin.output.csv' #or does it have to be .gz?
    pop_file = f"{code_path}/popgen/samples_to_populations.txt"

    window_size = '-w 50000'
    min_good_sites = '-m 5000'

    populations = [key for key in population_code.values()]

    pop_list = (" ").join([f'-p {pop}' for pop in populations])

    pop_gen_window_cmd = f'{geno_script_path}/popgenWindows.py \
                        -g {input_file} \
                        {window_size}\
                        {min_good_sites}\
                        {pop_list}\
                        -f phased\
                        -T 16\
                        --popsFile {pop_file}\
                        -o {output_file}'

    with open(f"{code_path}/popgen/scripts/popgenwindows/pop_windows.sh", "w") as file:
        file.write(f'{sha_bang}\n{sbatch_params}\
                   \n{pop_gen_window_cmd}')
    return(output_file)


def calculate_stats():
    input_file = run_popgen()

    #format = scaffold,start,end,mid,sites,pi_eyeless,pi_eyed,pi_surface,dxy_eyeless_eyed,dxy_eyeless_surface,dxy_eyed_surface,Fst_eyeless_eyed,Fst_eyeless_surface,Fst_eyed_surface
    pop_df = pd.read_csv(input_file, header = 0)
    print(pop_df.columns)

    #first calculate mean pi per population
    mean_dict = {}
    pops = ['eyeless','eyed', 'surface']
    for pop in pops:
        pi_mean = pop_df[[f'pi_{pop}']].mean().item()
        mean_dict[f'mean_pi_{pop}'] = pi_mean

    #pairwise
    for pair in list(combinations(pops, 2)):
        print(pair)
        ready_pair = f'{pair[0]}_{pair[1]}'
        dxy_mean = pop_df[[f'dxy_{ready_pair}']].mean().item()
        mean_dict[f'mean_dxy_{ready_pair}'] = dxy_mean

        fst_mean = pop_df[[f'Fst_{ready_pair}']].mean().item()
        mean_dict[f'mean_fst_{ready_pair}'] = fst_mean

    print(mean_dict)

    with open(f'{data_path}/popgen/pop_gen_windows_means.csv', 'w') as file:
            writer = csv.writer(file)
            # Write each key-value pair as a row
            for key, value in mean_dict.items():
                writer.writerow([key, value])


if __name__ == '__main__':
    #options are: make_pop_file, index_vcf_cmd or make_geno_cmd
    #format_setup()

    #run_popgen()

    calculate_stats()

