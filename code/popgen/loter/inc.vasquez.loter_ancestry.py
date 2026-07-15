#!/users/9/kell3262/miniforge3/envs/new_loter/bin/python

### This script runs Loter, a tool to estimate ancestry hybridization with haplotype variants
# See this tutorial: https://github.com/bcm-uga/Loter/blob/master/python-package/Local_Ancestry_Example.md
####### Steps ######
### 1. Convert genotypes into a numpy array
###### Adding additional populations increases the memory usage quite a lot, so the genotypes have to be processed in chunks
### 2. Run Loter on the chunks
### 3. Merge the post-loter chunks back together
### 4. Plot either by individual chromosomes or across the whole genome


import yaml #pyyaml
import allel #pip install scikit-allel
import numpy as np #careful with version
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
import gc
import os 
import psutil
import multiprocessing
import glob

#sys.path.append('../Loter/python-package')

module_path = '../Loter/python-package'
print(f"Here's the module path: {module_path}")

if module_path not in sys.path:
    sys.path.append(module_path)
print('path append done?')

import loter.locanc.local_ancestry as lc #v1.0 NOT 1.0.1 plz
print('import successful\n\n')

#Load config settings
with open("../../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome = config["GENOME"]
genome_path = config["GATK_GENOME_PATH"]

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/popgen/scripts/loter'

output_path_prefix = config["DATA_PATH"]
vcf_file_path = f'{output_path_prefix}/variant_calling/variants/pop_filter'

output_path = f'{output_path_prefix}/popgen/loter'
tmp_path = f'/scratch.global/kell3262'



######################################
# See troupbleshooting notes at the bottom


def log_memory(tag=""):
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024 / 1024
    print(f"[MEMORY] {tag}: {mem_mb:.2f} MB")

    #USAGE: log_memory(f'Memory while reading in {vcf_file}')


# Make a numpy file of the phased genotypes, this gets complicated for bigger datasets (eyeless, CM, and extras)
# So we break that into chunks with multiprocessing below
def vcf2npy(population, tmp_file, output_file):

    print(f'Make sure you chunk first! see data folder')
    if population == 'eyeless' or population == 'CM' or population == 'extras' or population == 'othersurface':

        if os.path.exists(output_file):
            print(f'File already made, skip')
        
        else:
            print(f'gathering chunked files for {population}')

            file_list = glob.glob(f'{vcf_file_path}/vcf_chunks/{genome}.phased.snps.{population}.*.vcf.gz')
            
            offset = 0
            memmap_path = tmp_file
            
            for chunk_file in file_list:  # List of VCF chunk paths
                log_memory(f'Memory before chunk_file: {chunk_file}')
                offset, n = append_chunk_to_memmap(chunk_file, memmap_path, offset) 
                #print(n)
                
            final_matrix = np.memmap(memmap_path, dtype=np.uint8, mode='r+', shape=(2*n, offset))
            final_matrix.flush()
            np.save(output_file, np.asarray(final_matrix))
            print(f'Saved to {output_file}')

    else:
        vcf_file = f'{vcf_file_path}/{genome}.phased.snps.{population}.vcf.gz.vcf.gz'
        if os.path.exists(output_file):
            print(f'File already made, skip')

        else:

            print(f'Reading in {vcf_file}')
            callset = allel.read_vcf(vcf_file, fields=['calldata/GT'])
            #print(f'callset keys: {callset.keys()}')
            #print(callset)
            print(callset['calldata/GT'].shape)
            
            haplotypes_1 = callset['calldata/GT'][:,:,0]
            haplotypes_2 = callset['calldata/GT'][:,:,1]
            
            m, n = haplotypes_1.shape
            log_memory(f'Memory while reading in {vcf_file}')
            
            
            #mat_haplo = np.empty((2*n, m), dtype=np.uint8)
            # instead of reading the whole vcf into a matrix, just immediately write it to a temp np.memmap file
            mat_haplo = np.memmap(tmp_file, dtype = np.uint8, mode = 'w+', shape = (2*n, m))

            mat_haplo[::2] = haplotypes_1.T
            mat_haplo[1::2] = haplotypes_2.T

            np.save(output_file, np.asarray(mat_haplo))
            
            print(f'read {vcf_file} and wrote matrix to {output_file}\n')

####### Checking callset / npy output shape/metrics ##########
# eyed: callset['calldata/GT'] = (53108128, 17, 2) (sites, individuals, haplos)
# eyeless: 
# surface: (53108128, 6, 2)

# translate to numpy file: 

### Some of these vcf files get killed with no error - eyeless and CM
## See variant_calling/variants/pop_filter for a script to chunk out those vcfs by chromosome
## `chunk_phased_eyeless.sh`
## This makes appends each chunk to the same memmap file

## Confirming this worked for eyeless, see `pop_filter/vcf_chunks/variant_num_check.sh`
## pre chunked sites: 53108128
## total sites after chunking: 53108128
def append_chunk_to_memmap(vcf_chunk, memmap_path, offset):
    callset = allel.read_vcf(vcf_chunk, fields=['calldata/GT'])
    print(callset['calldata/GT'].shape)
    gt = callset['calldata/GT']
    h1 = gt[:,:,0]
    h2 = gt[:,:,1]


    m, n = h1.shape  # m = variants, n = samples

    # If first chunk, create memmap file
    if offset == 0:
        total_snps_estimate = 10000000  # set a high enough estimate
        memmap = np.memmap(memmap_path, dtype=np.uint8, mode='w+', shape=(2*n, total_snps_estimate))
    else:
        memmap = np.memmap(memmap_path, dtype=np.uint8, mode='r+', shape=(2*n, offset + m))

    log_memory(f'Memory while reading in {vcf_chunk}')

    # Write chunk
    memmap[::2, offset:offset+m] = h1.T
    memmap[1::2, offset:offset+m] = h2.T

    #print(n)

    return ((offset + m),n)  # return new offset


# use multiprocess because it fully discards the memory before moving to the next run
def make_npy_memmap(populations):
    print(f'Making matrices for each population in {populations}\n')

    for population in populations:
        log_memory(f'Memory before reading {population}')

        #vcf_file = f'{vcf_file_path}/{genome}.phased.snps.{population}.vcf.gz.vcf.gz'
        tmp_file = f'{tmp_path}/int_matrix.{pop_to_run}.{population}.tmp'
        output_file = f'{output_path}/matrix.{pop_to_run}.{population}.saved.npy'

        print(f'Starting analysis for {population}')

        run_process = multiprocessing.Process(target = vcf2npy, args = (population, tmp_file, output_file))
        run_process.start()
        run_process.join()

        log_memory(f'After clean up for {population} in loop')
        print(f'Finished processing for {population}\n\n')


def run_loter(pop_to_run, chunk_index):
    print(f'Gathering matrices for {pop_to_run}')

    #fetch ind. nmpy arrays
    npy_dict = {}
    for population in populations:
        print(f'Reading in nmpy matrix for {population} and adding to dictionary')

        matrix_output = np.load(f'{output_path}/matrix.{pop_to_run}.{population}.saved.npy', mmap_mode = 'r')
        npy_dict[population] = matrix_output

        #log_memory(f'Memory while reading in {population}')
        
        print(f'for {population}: {matrix_output.shape} ; dictionary:{len(npy_dict)}, {len(npy_dict[population])}')

    #Loter (below) is killed silently. Like due to Loter reading in the memmap matrices fully and taking up all the RAM
    # need to do it in chunks of snps

    # Determine SNP count and chunk size
    snp_count = npy_dict[list(npy_dict.keys())[0]].shape[1]
    print(f'Number of SNPs: {snp_count}')
    chunk_size = 1_000_000
    num_chunks = (snp_count + chunk_size - 1) // chunk_size  # ceiling division
    
    #subset = slice(0,100000)
    #print(f'Running loter on {subset} of SNPs')

    #Define ancestry and admixed group by what pop_to_run
    # Set up ancestry and admixed groups
    if pop_to_run == 'normal':
        ancestry_keys = ['eyeless', 'surface']
        admixed_key = 'eyed'
        print('Eyed are the admixed pop')

    if pop_to_run == 'vasquez':
        ancestry_keys = ['vasquez', 'surface']
        admixed_key = 'CM'
        print('CM is the admixed pop')

    if pop_to_run == 'molino':
        ancestry_keys = ['molino', 'surface']
        admixed_key = 'CM'
        print('CM are the admixed pop')

    if pop_to_run == 'othersurface':
        ancestry_keys = ['othersurface', 'eyeless']
        admixed_key = 'eyed'
        print('Eyed are the admixed pop with othersurface & eyeless')

    log_memory('Memory before running loter')
    print('starting lc.loter_smooth')


    # Sbatch feeds the chunk index number
    output_file = f'{output_path}/{pop_to_run}_chunks/{pop_to_run}_output_chunk_{chunk_index:03d}.txt'

    start = chunk_index * chunk_size
    end = min((chunk_index + 1) * chunk_size, snp_count)
    subset = slice(start, end)

    print(f'\nRunning chunk {chunk_index+1}/{num_chunks} — SNPs {start} to {end}')
    log_memory(f'Memory before running chunk {chunk_index+1}')

    ancestry_group = [npy_dict[key][:, subset] for key in ancestry_keys]
    admixed_group = npy_dict[admixed_key][:, subset]

    res_loter = lc.loter_smooth(l_H=ancestry_group, h_adm=admixed_group, num_threads=2)
    print(f'\n Running: lc.loter_smooth(l_H={ancestry_group}, h_adm={admixed_group}, num_threads=2)')
    log_memory(f'Memory before writing chunk {chunk_index+1} to {output_file}')

    # Append results to file
    with open(output_file, 'a') as file:
        np.savetxt(file, res_loter, fmt="%i")

    print(f'finished, see {output_file}')
    return(npy_dict)


#### Merge all loter chunks back together horizontally
def merge_outputs(pop_to_run):
    
    # Set your parameters
    ####IF YOU RUN INTO A PROBLEM SET THIS BACK TO 54 and FIGURE IT OUT#######
    num_chunks = 53
    chunk_prefix = f'{output_path}/{pop_to_run}_chunks/{pop_to_run}_output_chunk_'
    chunk_files = [f'{chunk_prefix}{i:03d}.txt' for i in range(num_chunks)]


    for file in chunk_files:
        print(file)
    
    # Open all chunk files
    chunk_streams = [open(f, 'r') for f in chunk_files]
    output_file = f'{output_path}/{pop_to_run}_loter_output_phased.txt'

    # num_ind is number of haplotypes
    if pop_to_run == 'normal':
        num_ind = 34
    
    if pop_to_run == 'vasquez':
        num_ind = 72

    # num_ind is 78 if you have those extras in
    if pop_to_run == 'molino':
        num_ind = 72

    if pop_to_run == 'othersurface':
        num_ind = 34

    with open(output_file, 'w') as out:
        for row_idx in range(num_ind):  # assuming 34 rows
            row_pieces = []
            for f in chunk_streams:
                line = f.readline().strip()
                row_pieces.append(line)
            merged_line = ' '.join(row_pieces)
            out.write(merged_line + '\n')

    # Close all chunk files
    for f in chunk_streams:
        f.close()

    print(f"Saved merged output to {output_file}")



#Plots a single chromosome at a time
def positions_for_plotting(pop_to_run, hybrid, chunk_index):

    ### Need to have chromosomes / positions associated with the loter plots - Do we????
    
    # Right now, order of operations is: use vcf file to fetch chrosome and position
    # Then use the loter output to plot
    # if we use vcf files to get length of each chromosome, can we just subset the loter by that and output 1 chr/plot?
    # In which case, use the vcf file to fetch number of positions
    # Split the loter input by that number (load in only the first X that matches Chr Y)
    # Then plot and export

    #list of vcf chunk/chr files for hybrid population
    #print(chunk_index)

    
    ###############undo this - commented out only for running "original" loter##########
    vcf_files = glob.glob(f'{vcf_file_path}/vcf_chunks/{genome}.phased.snps.eyed.NC_*.vcf.gz')

    
    
    
    #print(vcf_files)
    #print(len(vcf_files))
    print(f'Using chromosome-chunked files to pull out positions for plotting.\n \
        Starting on {hybrid} \n \
        {vcf_files[chunk_index]}')


    # A dictionary of {original chr name: [start, end]}
    positions_dict = {}

    start = 0
    for file in vcf_files:
        #print(file)
        callset = allel.read_vcf(file, fields=['variants/CHROM', 'variants/POS'])

        #print(callset)
            #'variants/CHROM': array(['NC_064412.1', 'NC_064412.1', 'NC_064412.1', ..., 'NC_064412.1',
            #'NC_064412.1', 'NC_064412.1'], dtype=object), 'variants/POS': array([   13268,    39433,    39442, ..., 58282245, 58282246, 58282268],
            # dtype=int32)
        
        chromosome = callset['variants/CHROM'][0]
        position_length = len(callset['variants/POS'])
        position_range = [start, start + position_length - 1]

        #print(f'{chromosome}: positions: {position_length, position_range}')
        positions_dict[chromosome] = position_range

        #print(positions_dict)
        start += position_length

    #print(positions_dict)

    # chunk_index = what chr to make a plot for
    fetch_pos_range = list(positions_dict.items())[chunk_index]
    og_chr_name = list(positions_dict.keys())[chunk_index]
    start = fetch_pos_range[1][0]
    end = fetch_pos_range[1][1]

    print(og_chr_name, chunk_index, fetch_pos_range, start, end)

    print(np.arange(start, end))
    # read in loter output array
    loter_output = np.loadtxt(f"{output_path}/{pop_to_run}_loter_output_phased.txt", usecols=np.arange(start, end))

    print(loter_output.shape)

    # Group by columns, use the mean
    group_size = 20000
    num_groups = loter_output.shape[1] // group_size
    aggregated_cols = loter_output[:, :num_groups * group_size].reshape(loter_output.shape[0], num_groups, group_size).mean(axis=2)

    plt.figure(figsize=(12, 6), dpi=100)
    plt.imshow(loter_output, aspect='auto', cmap='coolwarm', interpolation='nearest')
    plt.colorbar(label='Value')
    plt.title('Eyeless vs. Surface assignment')
    plt.xlabel('SNPs')
    plt.ylabel('Haplotypes')
    plt.tight_layout()
    plt.show()
    plt.savefig(f"./loter_images/loter_{pop_to_run}_chr{chunk_index+1}.png")

#uses vcf file to get the positions where chromosome splits are
def make_chr_positions(pop_to_run, hybrid):
    vcf_files = sorted(glob.glob(f'{vcf_file_path}/vcf_chunks/{genome}.phased.snps.eyed.NC_*.vcf.gz'))

    print(f'Using chromosome-chunked files to pull out positions for plotting.\nStarting on {hybrid}')

    # this is a file with chr \t pos 
    print(f'{vcf_file_path}/eyed_chrom_pos.tsv')

    #### Original loter 10152025

    key_dict = {}
    with open(f'{genome_path}/chr_key.txt', 'r') as key_file:
        for line in key_file:
            line = line.rstrip()
            new_chrom, chrom_old = line.split() 
            key_dict[chrom_old]=new_chrom
   # print(key_dict)


    # Use a phased vcf to fetch positions to match to chromosomes for plotting
    # eyed_chrom_pos generated from bcftools query -f
    chrom_pos_dict = {}
    ignore_these = []

    with open(f'{vcf_file_path}/eyed_chrom_pos.tsv', 'r') as pos_file:
        for line in pos_file:
            line = line.rstrip()
            chrom_old, pos = line.split()
            #print(chrom_old, pos)

            if chrom_old in chrom_pos_dict:
                chrom_pos_dict[chrom_old].append(int(pos))
                #print(f'Adding entry to existing {chrom_old} key: {chrom_pos_dict}')
            
            if chrom_old not in chrom_pos_dict and chrom_old.startswith("NC_"):
                chrom_pos_dict[chrom_old] = []
                chrom_pos_dict[chrom_old].append(int(pos))
                #print(f'Adding {chrom_old}: {chrom_pos_dict}')
                   
            if chrom_old.startswith("NW_"):
                ignore_these.append(chrom_old)


    print(f'Chrom:pos dictionary made. \n Ignoring: {ignore_these}')


    positions_dict = {}
    all_keys_for_count = 0
    all_values_for_count = 0
    start_value = 0
    n=0

    #fetch absolute position values 0 - end of 
    for key, value in chrom_pos_dict.items():

        #print(f'old_chrom: {key}')

        new_chrom = key_dict[key]
        #print(new_chrom)
        
        pos_length = len(value)
        
        end_value = start_value + pos_length - 1
        print(f'{new_chrom}\tstart: {start_value}\tend: {end_value}')

        positions_dict[new_chrom] = [start_value, end_value]
        #print(positions_dict)

        # insert checks
        all_keys_for_count+=1
        all_values_for_count+= len(value)
        #print(all_keys_for_count, all_values_for_count)

        start_value += pos_length

        #dict format: {'chromosome':[start position, end position]}
        
    #print(positions_dict)
    print(f'Double checking: number of chromosomes = {all_keys_for_count}, number of positions = {all_values_for_count}')

    return(positions_dict)

def convert_txt_to_npy_in_chunks(loter_input, numpy_save_file, chunk_size=1000):
    data = []
    with open(loter_input, 'r') as f:
        for i, line in enumerate(f):
            row = np.array(line.strip().split(), dtype='float32')
            data.append(row)
            if (i + 1) % chunk_size == 0:
                print(f"Processed {i + 1} lines")

    data = np.array(data, dtype='float32')
    np.save(numpy_save_file, data)


# plots all chromosomes at once
def plot_all_chromosomes(pop_to_run, hybrid):

    print('Gathering chromosome positions to match to loter output')
    positions_dict = make_chr_positions(pop_to_run, hybrid)

    print('converting loter output to numpy array')
    log_memory('memory before converting loter to numpy array')
    loter_input = f"{output_path}/{pop_to_run}_loter_output_phased.txt"
    numpy_save_file = f"{output_path}/{pop_to_run}_loter_output_phased.npy"

    if os.path.exists(numpy_save_file):
        print(f'{numpy_save_file} already exists! continue')
    
    else:
        convert_txt_to_npy_in_chunks(loter_input, numpy_save_file)

    gap_width = 50
    slices = []
    ###### MAKE SLICES ###

    # Sort chromosomes to maintain consistent order
    #sorted_chromosomes = sorted(positions_dict.keys())
    sorted_chromosomes = positions_dict.keys()
    print(sorted_chromosomes)

    log_memory(f'Memory before loading output')
    mmap_array = np.memmap(numpy_save_file, dtype='float32', mode='r')
    #raw_array = np.loadtxt(f"{output_path}/{pop_to_run}_loter_output_phased.txt", dtype='float32')
    #np.save(f"{output_path}/{pop_to_run}_loter_output_phased.npy", raw_array)
    log_memory(f'Memory after loading and saving output as numpy file')

    for i, chr_name in enumerate(sorted_chromosomes):
        start, end = positions_dict[chr_name]
        log_memory(f'Memory before gathering slice: {chr_name}: positions {start} - {end}')


        # Load the slice from loter output
        mmap_array = np.load(f"{output_path}/{pop_to_run}_loter_output_phased.npy", mmap_mode='r')
        loter_slice = mmap_array[:, start:end + 1]
        np.save(f"{tmp_path}/{chr_name}_slice.npy", loter_slice)


        # Add a gap after each chromosome except the last
        if i < len(sorted_chromosomes) - 1:
            gap = np.full((loter_slice.shape[0], gap_width), np.nan)
            np.save(f"{tmp_path}/{chr_name}_gap.npy", gap)
        

    log_memory(f'Memory before starting slice appends')

    fig, ax = plt.subplots(figsize=(20, 8), dpi=100)
    x_offset = 0
    chrom_boundaries = []

    for chr_name in sorted_chromosomes:
        slice_path = f"{tmp_path}/{chr_name}_slice.npy"
        loter_slice = np.load(slice_path, mmap_mode='r')

        log_memory(f'In the slice: {chr_name}')

        ax.imshow(
            loter_slice,
            aspect='auto',
            cmap='coolwarm',
            interpolation='nearest',
            extent=[x_offset, x_offset + loter_slice.shape[1], 0, loter_slice.shape[0]]
        )

        chrom_boundaries.append((chr_name, x_offset))
        x_offset += loter_slice.shape[1]

        gap_file = f"{tmp_path}/{chr_name}_gap.npy"
        if os.path.exists(gap_file):
            gap = np.load(gap_file, mmap_mode='r')
            ax.imshow(
                gap,
                aspect='auto',
                cmap='coolwarm',
                interpolation='nearest',
                extent=[x_offset, x_offset + gap.shape[1], 0, gap.shape[0]]
            )
            x_offset += gap.shape[1]

    # Optional: add vertical lines and labels
    for chr_name, boundary in chrom_boundaries:
        ax.axvline(x=boundary, color='black', linestyle='--', linewidth=0.5)
        ax.text(boundary + 5, -5, chr_name, rotation=90, verticalalignment='bottom', fontsize=6)

    ax.set_xlim(0, x_offset)
    ax.set_ylim(0, loter_slice.shape[0])
    ax.set_title('Eyeless vs. Surface assignment across all chromosomes')
    ax.set_xlabel('SNPs')
    ax.set_ylabel('Haplotypes')
    plt.colorbar(ax.images[0], label='Value')
    plt.tight_layout()
    plt.savefig(f"./loter_images/loter_{pop_to_run}_all_chromosomes.png")
    plt.savefig(f"./loter_images/loter_{pop_to_run}_all_chromosomes.svg", format="svg")
    plt.show()
    

def plot_original_basic(pop_to_run):
    print(f'Making plots for loter output: {pop_to_run}')
    loter_output = np.loadtxt(f"{output_path}/{pop_to_run}_loter_output_phased.txt")

    plt.figure(figsize=(12, 6), dpi=100)
    plt.imshow(loter_output, aspect='auto', cmap='coolwarm', interpolation='nearest')
    plt.colorbar(label='Value')
    plt.title('Eyeless vs. Surface assignment')
    plt.xlabel('SNPs')
    plt.ylabel('Haplotypes')
    plt.tight_layout()
    plt.show()
    plt.savefig(f"./loter_images/loter_{pop_to_run}.png")



if __name__ == '__main__':
    # Use 'normal' or 'vasquez' or 'molino' to run 
    pop_to_run = sys.argv[1]
    chunk_index = int(sys.argv[2])
    
    print(f'Starting process for {pop_to_run} loter run\n')

    if pop_to_run == 'normal':
        populations = ['eyed', 'eyeless', 'surface']
        hybrid = 'eyed'
    
    if pop_to_run == 'vasquez':
        populations = ['CM', 'vasquez', 'surface']
        hybrid = 'eyed'

    if pop_to_run == 'molino':
        populations = ['CM', 'molino', 'surface']
        hybrid = 'CM'

    if pop_to_run == 'eyed_only':
        populations = ['eyed', 'molino', 'surface']
        hybrid = 'eyed'

    if pop_to_run == 'eyeless_only':
        populations = ['eyeless', 'molino', 'surface']
        hybrid = 'eyeless'

    if pop_to_run == 'original':
        hybrid = 'eyed'
    
    if pop_to_run == 'othersurface':
        populations = ['othersurface', 'eyeless', 'eyed']
        hybrid = 'eyed'

    ### Don't need to sbatch this with chunk index, just give it one chunk
    #make_npy_memmap(populations)

    # This runs with sbatch & slurm array task ID
    #I think you can run these together
    #run_loter(pop_to_run, chunk_index)

    #Don't run these together - otherwise it tries the merge 54 times lol - last one should do it.
    # run on command line with any chunk listed `./inc.vasquez.loter_ancestry.py othersurface 0`
    #merge_outputs(pop_to_run)

    #positions_for_plotting(pop_to_run, hybrid, chunk_index)

    plot_all_chromosomes(pop_to_run, hybrid)

    #plot_original_basic(pop_to_run)

