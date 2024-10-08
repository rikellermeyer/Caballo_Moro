#!/home/rk2643/miniforge3/envs/gatk/bin/python

###
# This makes a series of scripts (and SBATCH scripts) to run alignment
# fastq files > align > mark duplicate > rename to biological sample name
# Should be generalizable by editing the config.yaml file
# pay special attention to sequence id (fastq file name) and sample name
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
genome_path = config["GENOME_PATH"]
gatk_genome_path = config["GATK_GENOME_PATH"]
fastq_dir = config["M_FASTQ_PATH"] # this is the adding molino one, edit for MOLNGs

sha_bang = '#!/usr/bin/bash'
script_path = f'{code_path}/alignment/scripts/molino2'
script_errors = f'{script_path}/error_logs' # idk if these are necessary or should go to slurmout/*.err
mem2_genome_path = f'{genome_path}/{genome}.fna' #for bwa-mem2
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


### parse data for reads/lanes/ect ###
### This is the kinda beta way to do this, next time use the config["SAMPLES"] dictionary
    ##Stowers: n_{lane}_{read}_{barcode}.fastq.gz  w/ barcode specific to the individual and is consistent across flowcells within one MOLNG

    ##UMN: ***_R1|2_001.fastq.qz (What is I1?/L1?)
    ## for molino running ~12/2023, SRR#_1|2.fastq.gz (1/2 is reads for each individual)

def read_group_info(fastq_dir):
    ##Detect file naming structure and provide seq_id:sample_name

    dir_structure = glob.glob(f'{fastq_dir}/*.fastq.gz')
    base_id = [os.path.basename(file) for file in dir_structure]

    seq_pairs = {}
    for name in base_id:
        find_seq_names = re.search(r'(.+)_(\d)', name) #right now this only works for SRR#_1.fastq.gz
        seq_id= find_seq_names[1]
    
        if seq_id not in seq_pairs:
            seq_pairs[seq_id]=[]
            seq_pairs[seq_id].append(name)
        else:
            seq_pairs[seq_id].append(name)
    return(seq_pairs)


### Alignment ###
def alignment(seq_id, n, x):
    output_path = f'{data_path}/alignment/bam'
    output_file= f'{output_path}/{seq_id}.sorted.bam'
    inter_out = f'{output_path}/{seq_id}.unsorted.bam'
    #sorted_output_file= f'{output_path}/{seq_id}.sorted.bam'
    check_paths([output_path, script_path, script_errors])

    script_name = f'{seq_id}_bwa.sh'
    
    sample_pairs_dict = read_group_info(fastq_dir)
    #each individual/sample ID gets an alignment script with associated read groups
    files = sample_pairs_dict[seq_id]
    read1=f'{fastq_dir}/{files[0]}'
    read2=f'{fastq_dir}/{files[1]}'

    read_group = f'@RG\\tID:{seq_id}\\tSM:M{x}\\tPL:illumina' # This works for generic formats but should be done with a custom script for complicated sets
    # See cmc/code/alignment/readgroups.py for Stowers read group generation

    bwa_cmd = f'bwa-mem2 mem -t 16 -M -R "{read_group}" {mem2_genome_path} {read1} {read2} | samtools view -b > {inter_out}' 
    sort_cmd = f'samtools sort -@ 8 -o {output_file} -T {seq_id} {inter_out}'
    ### STUPID PROBLEM: calling samtools view: -b indicates .bam output... not the output file name!! Need >
    
    bwa_cmds = [script_name, read_group]

    with open(f'{script_path}/{script_name}', 'w') as file:
        file.write(f'{sha_bang}\n{bwa_cmd} &&\n \
                    echo "bwa complete"\n \
                   {sort_cmd} &&\n \
                   echo "sort complete"')
                   #\n{sort_cmd}\n echo "sort complete"')
   
    logging.info(f'{seq_id}\n {read_group} \n {bwa_cmd}')

    return(bwa_cmds)
        


#Runs picards MarkDuplicates to mark and consolidate between barcodes/samples.
def mark_duplicates(seq_id):
    input_file = f'{data_path}/alignment/bam/{seq_id}.sorted.bam'
    output_path = f'{data_path}/alignment/marked_dups'
    temp_space = f'{data_path}/scratch'
    check_paths([output_path, temp_space])


    metrics_file = f'{output_path}/{seq_id}.metrics.txt'
    output_file = f'{output_path}/{seq_id}.marked.bam'
    error_file = f'{script_errors}/{seq_id}.markdups.err'

    #seems to be a picard problem with the cerebro module, unsure about conda install
    pre_markdups_run = f'touch {metrics_file} \n touch {output_file}'

    mark_dups_cmd = f'picard MarkDuplicates \
                    -Xmx20g\
                    I={input_file}\
                    O={output_file}\
                    METRICS_FILE={metrics_file}\
                    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 \
                    SORTING_COLLECTION_SIZE_RATIO=0.05\
                    ASSUME_SORTED=true\
                    CREATE_INDEX=TRUE\
                    TMP_DIR={temp_space}\
                    2> {error_file}'

    mark_dups_check = f"""if "grep -e 'fail' -e 'error' -e 'killed' -e 'FAIL' -e 'ERROR' {error_file}"; then \n echo 'check on {output_file}' > {error_file}.lookout.errors \n exit 1 \n fi"""
    
    with open(f'{script_path}/{seq_id}_markdups.sh', 'w') as file:
        file.write(f'{sha_bang}\n{pre_markdups_run}\n{mark_dups_cmd}')
    
    mark_dups_cmds = [f'{seq_id}_markdups.sh', mark_dups_check]

    return(mark_dups_cmds)


#Renames SRRs to sample name, ie M1/2/3
def rename_to_sample(seq_id, n, x):
    input_file_name = f'{data_path}/alignment/marked_dups/{seq_id}.marked.bam'
    output_path = f'{data_path}/alignment/renamed'
    check_paths([output_path])

    item_list = alignment(seq_id, n, x)
    RG =  item_list [1]
    sample_name = re.search(r'.+\\tSM:(.+)\\t.+', RG)[1]

    new_file_name = f'{output_path}/{sample_name}.marked.bam'

    mv_cmd = f'set -e \n cp {input_file_name} {new_file_name}'

    with open(f'{script_path}/{seq_id}_rename.sh', 'w') as file:
        file.write(f'{sha_bang}\n{mv_cmd}')

    return(f'{seq_id}_rename.sh')


#############################################################
# Total alignment script by individual: alignment/bwa > mark
#bwa: {SRR:[script_name, bwa_check, read_group]}
#touch the picard outputs, otherwise it'll fail on slurm
#mark_dups {SRR:[mark_dups_cmd, mark_dups_check]}
#SRR to sample name, check
def batch_align_script(fastq_dir):
    fastq_info = read_group_info(fastq_dir)
    seq_ids = fastq_info.keys()

    array_length = len(seq_ids)

    n=1
    x=4 #for molino RG naming
    for seq_id in seq_ids:
        bwa_cmd = alignment(seq_id, n, x)
        bwa_script = bwa_cmd[0]
        mark_dups_cmds = mark_duplicates(seq_id)
        mark_dups_script = mark_dups_cmds[0]

        rename_script=rename_to_sample(seq_id, n, x)

        with open(f'{script_path}/batch_align_{n-1}.sh', 'w') as file:
            file.write(f'{sha_bang}\n\nsh {bwa_script}&&\nsh {mark_dups_script}&&\nsh {rename_script}')
        
        n +=1
        x +=1


    batch_header = sbatch_header(f'big_align', array_length)
    with open(f'{script_path}/big_align_batch.sh', 'w') as big_file:
        big_file.write(f'{sha_bang}\n{batch_header}\n\nsh batch_align_${{SLURM_ARRAY_TASK_ID}}.sh')
    

if __name__ == '__main__':
    print("Making a script to go from fastq > marked, renamed bams")
    print(f"Using files in {fastq_dir}")
    batch_align_script(fastq_dir)
    

