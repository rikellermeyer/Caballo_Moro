#!/home/rk2643/miniforge3/envs/CMpop/bin/python

#This is a script to set up and run popgenWindows.py from https://github.com/simonhmartin/genomics_general#processing-vcf-files
#Use on cerebro, conda env CMpop
#popgenWindows from above has to be in the current directory for this script

import yaml
import subprocess as sbp

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
        \n#SBATCH --time="12-00:00" \
        \n#SBATCH --error=./slurmout/popgen.%a.%A.err \
        \n#SBATCH --output=./slurmout/popgen.%a.%A.out \
        \n#SBATCH --mail-user=rk2643@stowers.org \
        \n#SBATCH --mail-type=FAIL \
        \n#SBATCH --mail-type=END'

#location of CMC all.vcf files
vcf_path = config["CM_VCF_PATH"]

#smartin's genomics scripts
geno_script_path = f'{code_path}/popgen/smartin_genomics'

### STEPS ###
# install, add to path
# parse VCF with smartin_genomics/parseVCF.py which takes .all.vcf.gz file and makes a .geno file
# .geno file is a "genotype matrix" specific to the genomics_general code from smartin
# gzip .geno file > .geno.gz
# generate a population file "sample_name \t population_name \n"
# run popgenWindows.py

def format_setup(which_cmd):

    ## first make a --popsFile 
    # format is sample \t population \n
    if which_cmd == "make_pop_file":
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
    #use subprocess to do the vcf
    input_vcf = f'{vcf_path}/{genome}.all.filtered.vcf.gz'
    output_geno = f'{data_path}/popgen/{genome}.all.geno.gz'

    if which_cmd == "index_vcf_cmd":
        sbp_cmd = [f'bcftools index {input_vcf}']
        print(sbp_cmd)

    #note that this needs an index (with bcftools index)
    if which_cmd == "make_geno_cmd":
        sbp_cmd = f'{geno_script_path}/VCF_processing/parseVCFs.py -i {input_vcf} --threads 16 2> mkgeno.err\
               | bgzip > {output_geno}'
        with open("make_geno_file.sh", "w") as file:
            file.write(f'{sha_bang}\n{sbatch_params}\
                   \n{sbp_cmd}')
    
    if which_cmd == "index_vcf_cmd" or which_cmd == "make_geno_cmd":
        output_sbp = sbp.run(sbp_cmd, shell = True, check = True)
        print(output_sbp.stdout)
        if sbp.CalledProcessError == True :
            print(f'Errors: {output_sbp.stderr}')
    

def run_popgen():
    print('yes for the popgen')

    input_file = f'{data_path}/popgen/{genome}.all.geno.gz'
    output_file = f'{data_path}/popgen/{genome}.50kwin.output.csv.gz'
    pop_file = "samples_to_populations.txt"

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

    with open(f"{code_path}/popgen/scripts/pop_windows.sh", "w") as file:
        file.write(f'{sha_bang}\n{sbatch_params}\
                   \n{pop_gen_window_cmd}')
    return(output_file)

    #From Rachel
    #python popgenWindows.py -w 50000 -m 5000 -g 
    # OL_vs_NL_biSNPs_wInvar.geno.gz -o May2022_OL_vs_NL_50k_popgenWindows.output.csv.gz 
    # -f phased -T 24  -p NLC -p NLS -p OLC -p OLS --popsFile OL_NL_samples_tabs.txt\
    ### explanation
    #-w = window size sites in bases; -m = minimum good sites for a window; 
    # -g = input.geno.gz; -o= outputfile.csv.gz' -f= phased/haplo/diplo/pairs; 
    # T = threads; -p = population names; --popsFile = see above

# Average dxy, pi, fst across all sites, 
# scaffold,start,end,mid,sites,
#pi_eyeless,pi_eyed,pi_surface,
#dxy_eyeless_eyed,dxy_eyeless_surface,dxy_eyed_surface,
#Fst_eyeless_eyed,Fst_eyeless_surface,Fst_eyed_surface
def average_across_sites():
    csv_file = run_popgen()
    


if __name__ == '__main__':
    print('ok')
    #options are: make_pop_file, index_vcf_cmd or make_geno_cmd
    #format_setup('make_pop_file')
    #format_setup('index_vcf_cmd')
    #format_setup('make_geno_cmd')

    run_popgen()

