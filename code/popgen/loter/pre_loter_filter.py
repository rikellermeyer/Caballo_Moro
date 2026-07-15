#!/users/9/kell3262/miniforge3/envs/gatk/bin/python

### This takes an all sites vcf file and parses it for Loter with gatk
## Filters to SNP / populations
## and phases with beagle


import yaml


#Load config settings
with open("../../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome_path = config["GATK_GENOME_PATH"]
genome_file = f'{genome_path}/Amex3.0_surface.fna'

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/popgen/scripts/loter'

data_path = config["DATA_PATH"]

all_vcf_input = f'{data_path}/variant_calling/variants/pop_filter/CabMoroProject_ALLpopulations_ALLsites.vcf.gz'
vcf_int_output = f'{data_path}/variant_calling/variants/pop_filter' #10152025

sha_bang = "#!/usr/bin/bash"

def add_sbatch(script):
    batch_info=f'#!/usr/bin/bash\n \
    \n#SBATCH --job-name=pre_loter \
    \n#SBATCH --array=0-4%5 \
    \n#SBATCH --cpus-per-task=4 \
    \n#SBATCH --time="15:00:00" \
    \n#SBATCH --mem=64G \
    \n#SBATCH --mail-user=kell3262@umn.edu \
    \n#SBATCH --mail-type=FAIL \
    \n#SBATCH --mail-type=END\
    \n#SBATCH --output=./slurmout/pre_loter.{script}.%A.%a.out \
    \n#SBATCH --error=./slurmout/pre_loter.{script}.%A.%a.err'
    return(batch_info)




### Make files for loter
# First, select variants for population: eyed, eyeless, sf, vasquez, all
def pop_filter(pop_keys):
    print('making scripts to filter by population')

    n=0
    for pop, expr in pop_keys.items():
        print(pop, expr)
        sv_cmd = f'gatk SelectVariants\
        -R {genome_file}\
        -V {all_vcf_input}\
        -O {vcf_int_output}/Amex3.0_surface.filtered.snps.{pop}.vcf.gz \
        --sample-expressions {expr} \
        --select-type-to-include SNP' #.snps = 10152023 #running to get original loter 10152025 - incl snps here
        print(sv_cmd)

        total_cmd = f'{sha_bang}\n\n\
        {sv_cmd}'
        #print(total_cmd)
        
        with open(f'{code_path}/original_preloter_filter.{pop}.{n}.sh', 'w') as file:
            file.write(total_cmd)
        n+=1

    sbatch_script_cmd = f'{add_sbatch("pop_filter")}\n\n\
    sh original_preloter_filter.*.${{SLURM_ARRAY_TASK_ID}}.sh'
    print(sbatch_script_cmd)

    with open(f'{code_path}/original_popfilter_batch.sh', 'w') as file:
        file.write(sbatch_script_cmd)


#Loter uses SNPs... not whole genome. Filter to SNPs fml
def select_snps(phasing_pops):
    print('making snp filter scripts')

    n=0
    for pop in phasing_pops:
        input_file = f'{vcf_int_output}/Amex3.0_surface.filtered.{pop}.vcf.gz'
        output_file = f'{vcf_int_output}/Amex3.0_surface.filtered.snps.{pop}.vcf.gz'

        sv_cmd = f'gatk SelectVariants\
        -R {genome_file}\
        -V {input_file}\
        -O {output_file} \
        --select-type-to-include SNP'

        
        with open(f'{code_path}/snp_preloter_filter.{pop}.{n}.sh', 'w') as file:
            file.write(f'{sha_bang}\n\n{sv_cmd}')
        n+=1

    sbatch_script_cmd = f'{add_sbatch("snp_filter")}\n\n\
    sh snp_preloter_filter.*.${{SLURM_ARRAY_TASK_ID}}.sh'
    print(sbatch_script_cmd)

    with open(f'{code_path}/snpfilter_batch.sh', 'w') as file:
        file.write(sbatch_script_cmd)



def phasing_scripts(phasing_pops):
    print('making phased scripts')

    n=0
    for pop in phasing_pops:
        input_file = f'{vcf_int_output}/Amex3.0_surface.filtered.snps.{pop}.vcf.gz'
        output_file = f'{vcf_int_output}/Amex3.0_surface.phased.snps.{pop}.vcf.gz'
        beagle_cmd = f"java -Djava.io.tmpdir=/scratch.global/kell3262/ \
        -Xmx32g -jar beagle.27May24.118.jar nthreads=16 \
        gt={input_file} \
        out={output_file}"
        with open(f'{code_path}/original_phasing_script.{pop}.{n}.sh', 'w') as file:
            file.write(f'{sha_bang}\n\n{beagle_cmd}')
        n+=1

    sbatch_script_cmd = f'{add_sbatch("beagle_phasing")}\n\n \
    sh original_phasing_script.*.${{SLURM_ARRAY_TASK_ID}}.sh'

    print(sbatch_script_cmd)

    #1015205
    with open(f'{code_path}/original_phasing_batch.sh', 'w') as file:
        file.write(sbatch_script_cmd)


if __name__ == '__main__':
    #pop_keys = {'eyeless':'"C\\d+"','eyed':'E+','surface':'S+', 'vasquez':'V+', 'CM':'"C\\d+" --sample-expressions E+'}
    pop_keys = {'eyeless':'"C\\d+"','eyed':'E+','surface':'S+'} #10152025

    print('making pop_filtered files')
    pop_filter(pop_keys)
    #print('STOP! You may need to adjust the regular expressions for gatk SelectVariants \n\
    #For example: when using regular expressions (ie. eyeless C\d+), use quotes!\n\
    #To designate two populations with their own regex - use --sample-expressions TWICE')

    #phasing_pops = ["eyed", "eyeless", "surface", "vasquez", "CM"]
    phasing_pops = ["eyed", "eyeless", "surface"] #10152025

    
    #print('select out snps')
    #select_snps(phasing_pops)
    #print('NEXT TIME DO SNP FILTERING AND POP FILTERING AT THE SAME TIME')
    
    print('making beagle scripts')
    phasing_scripts(phasing_pops)
    print('USE LOTER ENV!!!!!!!')