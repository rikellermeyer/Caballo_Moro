#!/users/9/kell3262/miniforge3/envs/new_loter/bin/python

### This script takes the loter files generated in inc.vasquez.loter_ancestry.py
## And calculates the statistics: tract length summaries, haplotype summaries, and confirms calculations

import yaml
import allel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.image
from PIL import Image
import glob
import pandas as pd
import statistics
from scipy.stats import sem
import sys
from joblib import Parallel, delayed
import psutil
import os
import csv
from collections import defaultdict

with open("../../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/popgen/scripts/loter'

data_path_prefix = config["DATA_PATH"]

output_path=f'{data_path_prefix}/popgen/loter'

#sys.path.append('./Loter/python-package')
loter_script_path=f'./Loter'


def log_memory(tag=""):
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024 / 1024
    print(f"[MEMORY] {tag}: {mem_mb:.2f} MB")


# the loter output is a big text file
# with rows = haplotypes (eyed = 17 ind x2 = 34 ; CM = 38 ind x2 = 72)
# and columns = 0 (eyeless) or 1 (eyed) for each SNP (53108128)

### First, find tract lengths - continuous SNPs for 0's and 1's
def analyze_loter_tracts(pop_to_run):

    def row_generator(input_file):
        """Yields one row at a time as a NumPy array of ints."""
        with open(input_file, 'r') as file:
            for line in file:
                yield np.fromstring(line.strip(), sep=' ', dtype=int)

    def compute_tract_lengths(row):
        """Computes tract lengths for a single row of ancestry calls."""
        tracts = []

        # Defines the first tract identity as 0/1
        current_ancestry = row[0]
        length = 1

        for val in row[1:]:

            # if next value matches original value, add to the number of SNPs
            if val == current_ancestry:
                length += 1

            # if the next value switches, make the next track
            else:
                tracts.append((current_ancestry, length))
                current_ancestry = val
                length = 1
        tracts.append((current_ancestry, length))  # Final tract
        return tracts

    def process_loter_file(input_file, output_file):
        """Processes the Loter output file row-by-row and computes tract stats."""
        all_tracts = []
        ancestry_stats = defaultdict(lambda: {
            'count': 0,
            'total_length': 0,
            'max_length': 0
        })

        # For every row, make tract lengths (1:333, 0:434, 1:756, ...)
        for row_index, row in enumerate(row_generator(input_file)):
            print(f'Computing tracts for row {row_index}')
            tracts = compute_tract_lengths(row)
            all_tracts.append(tracts)

            # Update ancestry stats
            for ancestry, length in tracts:
                stats = ancestry_stats[ancestry]
                stats['count'] += 1
                stats['total_length'] += length
                stats['max_length'] = max(stats['max_length'], length)

            # Write tract lengths to file
            with open(output_file, 'a') as out:
                tract_str = ' '.join([f"{ancestry}:{length}" for ancestry, length in tracts])
                out.write(f"Ind{row_index+1}: {tract_str}\n")

        # Compute averages in total
        for ancestry, stats in ancestry_stats.items():
            stats['average_length'] = stats['total_length'] / stats['count']

        return all_tracts, ancestry_stats

    def write_summary_to_csv(stats, summary_csv):
        """Writes ancestry summary statistics to a CSV file."""
        with open(summary_csv, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Ancestry", "Total Tracts", "Total Length", "Average Length", "Max Tract Length"])
            for ancestry in sorted(stats):
                s = stats[ancestry]
                writer.writerow([
                    ancestry,
                    s['count'],
                    s['total_length'],
                    f"{s['average_length']:.2f}",
                    s['max_length']
                ])

    
    input_file=f'{output_path}/{pop_to_run}_loter_output_phased.txt'
    output_file = f"{output_path}/{pop_to_run}_loter_tract_lengths.txt"
    summary_csv = f"{output_path}/{pop_to_run}_loter_tract_summary.csv"

    log_memory(f'\nMemory before reading in {input_file}: ')

    tracts, stats = process_loter_file(input_file, output_file)
    log_memory(f'\nMemory after computing tract length and stats: ')

    write_summary_to_csv(stats, summary_csv)

    print(f"\nSummary statistics written to: {summary_csv}")


### In contrast to above, this generates stats by row and inputs them into `haplotype_summary.csv`
def haplo_stats(pop_to_run):
    def row_generator(input_file):
        """Yields one row at a time as a NumPy array of ints."""
        with open(input_file, 'r') as file:
            for line in file:
                yield np.fromstring(line.strip(), sep=' ', dtype=int)

    def compute_tract_lengths(row):
        """Computes tract lengths for a single row of ancestry calls."""
        tracts = []
        current_ancestry = row[0]
        length = 1
        for val in row[1:]:
            if val == current_ancestry:
                length += 1
            else:
                tracts.append((current_ancestry, length))
                current_ancestry = val
                length = 1
        tracts.append((current_ancestry, length))  # Final tract
        return tracts

    def process_loter_file(input_file, summary_csv):
        print(f'Writing to new summary file: {summary_csv}')

        # Prepare haplotype summary CSV
        with open(summary_csv, 'w', newline='') as hap_file:
            hap_writer = csv.writer(hap_file, delimiter='\t')

            #before re-formatting
            #hap_writer.writerow([
            #    'Haplotype_ID',
            #    'Num_0_Tracts', 'Max_0_Length', 'Total_0_Length', 'Mean_0_Length', 'Std_0_length',
            #    'Num_1_Tracts', 'Max_1_Length', 'Total_1_Length', 'Mean_1_Length', 'Std_1_length', 'Percent_of_eyeless_ancestry', 'Tadmix'
            #])

            hap_writer.writerow([
                'Haplotype_ID, Num_0_Tracts, Min_0_Length, Max_0_Length, Total_0_Lengths, Mean_0_Length\
                , Std_0_Lengths, Num_1_Tracts, Min_1_Length, Max_1_Length, Total_1_Lengths, Mean_1_Length, Std_1_Lengths, Percent_0_ancestry, Tadmix'
            ])

            #length of all 0 and 1 tracts per haplotype
            all_tract_len_0 = []
            all_tract_len_1 = []
            all_Tadmix = []
            means_tracts_0 = []
            means_tracts_1 = []

            eyedCM_tracts_only = []
            eyelessCM_tracts_only = []

            for row_index, row in enumerate(row_generator(input_file)):

                # For every row, fetch 0:length of 0, 1:length of 1 and do stats
                tracts = compute_tract_lengths(row)

                # Separate tract lengths by ancestry

                #make a new dictionary that automatically makes an empty list value when you add a new key
                # then add ancestry and tract length as an entry {0:[2455, 54393], 1:[4555, 87563]}
                tract_lengths_by_ancestry = defaultdict(list)
                for ancestry, length in tracts:
                    tract_lengths_by_ancestry[ancestry].append(length)
                
                #print(f'Number of tracts in Hap{row_index+1} dictionary: \
                #    {len(tract_lengths_by_ancestry)}')

                ### then, fetch all of the entries with key value 0
                # dict.get(keyfor0entries, [the list of values in question])
                # Stats for ancestry 0 - eyeless

                #list of all tract length
                lengths_0 = tract_lengths_by_ancestry.get(0, [])
                #print(lengths_0)
                #the number of tract lengths
                num_0 = len(lengths_0)
                #the longest tract length
                max_0 = max(lengths_0) if lengths_0 else 0
                #the shortest tract length
                min_0 = min(lengths_0) if lengths_0 else 0
                #the sum of all tract lengths
                total_0 = sum(lengths_0)
                #mean tract length 
                mean_0 = total_0 / num_0 if num_0 > 0 else 0
                #std of tract lengths
                std_0 = np.std(lengths_0)

                #add the list of tract lengths [213, 2333, ...] to a master list of tract lengths for 0
                all_tract_len_0.extend(lengths_0)
                means_tracts_0.append(mean_0)
                print(f'all_tract_len_0, tail: {all_tract_len_0[-3:]}; length:{len(all_tract_len_0)}')
                print(f'means_tracts_0, tail: {means_tracts_0[-3:]}; length:{len(means_tracts_0)}')
                #print(all_tract_len_0)

                # Stats for ancestry 1 - surface
                lengths_1 = tract_lengths_by_ancestry.get(1, [])
                num_1 = len(lengths_1)
                min_1 = min(lengths_1) if lengths_1 else 0
                max_1 = max(lengths_1) if lengths_1 else 0
                total_1 = sum(lengths_1)
                mean_1 = total_1 / num_1 if num_1 > 0 else 0
                std_1 = np.std(lengths_1)

                all_tract_len_1.extend(lengths_1)

                means_tracts_1.append(mean_1)


                # Proportion stats
                # This calculates the total length/number of zero vs. the total of all the lengths
                # aka the proportion of the snps that are 0's by haplotype
                percent_of_0s = (total_0 / (total_0 + total_1)) * 100 
                percent_of_1s = (total_1 / (total_0 + total_1)) * 100


                if eyeless_pop_index[0] <= row_index <= eyeless_pop_index[1]:
                    print(f'Eyeless row: {row_index}')
                    eyelessCM_tracts_only.append(percent_of_0s)

                if eyed_pop_index[0] <= row_index <= eyed_pop_index[1]:
                    print(f'Eyed row: {row_index}')
                    eyedCM_tracts_only.append(percent_of_0s)

                #Tadmix by haplotype:
                #The number of generations since the onset of admixture (T admix ) 
                # can be estimated using the following equation:
                # where LM is the average ancestry tract length from the minor
                # parent in Morgans and pB is the proportion of the genome derived 
                # from the major parent (the probability of recombining) (31–33).
                #Tadmix = 1/LM*PB
                # PB(eyeless) = 0.88
                # bp to Morgans 1.16 cM/Mb (0.0000000116 Morgan/bp) 
                # LM = surface average bp * 1.13e-8

                if total_0 > total_1:
                    #1's are the minor ancestor
                    LM = mean_1*0.0000000116
                    PB = percent_of_0s / 100
                    Tadmix = 1/(LM*PB)
            
                if total_1 > total_0:
                    #0's are the minor ancestor
                    LM = mean_0*0.0000000116
                    PB = percent_of_1s / 100
                    Tadmix = 1/(LM*PB)

                all_Tadmix.append(Tadmix)


                #'Haplotype_ID\tNum_0_Tracts\tMin_0_Length\tMax_0_Length\tTotal_0_Lengths\tMean_0_Length\
                #\tStd_0_Lengths\tNum_1_Tracts\tMin_1_Length\tMax_1_Length\tTotal_1_Lengths\tMean_1_Length\tStd_1_Lengths\tPercent_0_ancestr\tTadmix'
            
                write_list = [f"Hap{row_index+1},{num_0},{min_0},{max_0},{total_0},\
            {mean_0:.2f},{std_0:.2f},{num_1},{min_1},{max_1},{total_1},{mean_1:.2f},\
            {std_1:.2f},{percent_of_0s:.2f},{Tadmix:.2f}"]

                
                hap_writer.writerow(write_list)


                print(f'Summary written for Hap{row_index+1}')

            # this calculates the average tract length of zeros and 1s across all samples 
            # and the standard deviation of all the samples compared to the mean
            mean_length_of_all_0s = sum(all_tract_len_0) / len(all_tract_len_0)
            std_of_all_0s = np.std(means_tracts_0)

            mean_length_of_all_1s = sum(all_tract_len_1) / len(all_tract_len_1)
            print('surface tracts, total:', all_tract_len_1[0:10], sum(all_tract_len_1), len(all_tract_len_1))
            std_of_all_1s = np.std(means_tracts_1)
            print('std of sf tracts', std_of_all_1s)

            #this calculates the percentage of 0's compared to all the SNPs across all the haplotypes
            percent_of_0s_in_all_samples = (sum(all_tract_len_0) / (sum(all_tract_len_0) + sum(all_tract_len_1))) * 100
        
            percent_of_1s_in_all_samples = (sum(all_tract_len_1) / (sum(all_tract_len_0) + sum(all_tract_len_1))) * 100

            average_percent_0s_in_eyelessCM = (sum(eyelessCM_tracts_only) / len(eyelessCM_tracts_only)) 
            std_of_percent_0s_in_eyelessCM = np.std(eyelessCM_tracts_only)

            average_percent_0s_in_eyedCM = (sum(eyedCM_tracts_only) / len(eyedCM_tracts_only)) 
            std_of_percent_0s_in_eyedCM = np.std(eyedCM_tracts_only)

            ## Calculate average Tadmix (time since admixture)
            mean_Tadmix = sum(all_Tadmix) / len(all_Tadmix)
            Admix_SE = sem(all_Tadmix)
            print(f'Tadmix: {Tadmix} +- {Admix_SE}')


            #Write full summary stats to file
            if pop_to_run == 'normal' or pop_to_run == 'othersurface':
                hap_writer.writerow([f'\nMean eyeless tract length: {mean_length_of_all_0s:.2f} \n\
            SE of eyeless tract lengths: {std_of_all_0s:.2f} \n \
            Mean surface tract length: {mean_length_of_all_1s:.2f} \n \
            SE of surface tract lengths: {std_of_all_1s:.2f} \n \
            Total percent of 0s: {percent_of_0s_in_all_samples:.2f} \n \
            SE percent of 0s: {std_of_percent_0s_in_eyedCM} \n\
            Time since admixture (Tadmix): {mean_Tadmix:.2f} \n \
            Tadmix SEM: {Admix_SE:.2f} \n \n \
            Percent of eyeless tracts in eyed CM fish (17 fish): {average_percent_0s_in_eyedCM:.2f} \n \t \
            SD of percent 0s in eyed fish: {std_of_percent_0s_in_eyedCM:.2f}'])
            else:
                hap_writer.writerow([f'\nMean eyeless tract length: {mean_length_of_all_0s:.2f} \n\
                SE of eyeless tract lengths: {std_of_all_0s:.2f} \n \
                Mean surface tract length: {mean_length_of_all_1s:.2f} \n \
                SE of surface tract lengths: {std_of_all_1s:.2f} \n \
                Total percent of 0s: {percent_of_0s_in_all_samples:.2f} \n \
                Time since admixture (Tadmix): {mean_Tadmix:.2f} \n \
                Tadmix SEM: {Admix_SE:.2f} \n \n \
                Percent of eyeless tracts in eyeless CM fish (19 fish): {average_percent_0s_in_eyelessCM:.2f} \n \t \
                SD of percent 0s in eyeless fish: {std_of_percent_0s_in_eyelessCM:.2f}\n \
                Percent of eyeless tracts in eyed CM fish (17 fish): {average_percent_0s_in_eyedCM:.2f} \n \t \
                SD of percent 0s in eyed fish: {std_of_percent_0s_in_eyedCM:.2f}'])

            
    #Run all the functions!
    input_file=f'{output_path}/{pop_to_run}_loter_output_phased.txt'
    summary_csv = f"{output_path}/{pop_to_run}_haplotype_summary.csv"

    log_memory(f'\nMemory before reading in {input_file}: ')

    process_loter_file(input_file, summary_csv)
    log_memory(f'\nMemory after computing tract length and stats: ')

    print(f"\nSummary statistics written to: {summary_csv}")
   

def test_output_match(pop_to_run):


    def row_generator(file_path):
        with open(file_path, 'r') as file:
            for line in file:
                yield np.fromstring(line.strip(), sep=' ', dtype=int)

    def compute_tract_lengths(row):
        tracts = []
        current_ancestry = row[0]
        length = 1
        for val in row[1:]:
            if val == current_ancestry:
                length += 1
            else:
                tracts.append((current_ancestry, length))
                current_ancestry = val
                length = 1
        tracts.append((current_ancestry, length))
        return tracts

    def match_loter_to_tracts(input_file, num_individuals=5):
        for idx, row in enumerate(row_generator(input_file)):
            if idx >= num_individuals:
                break
            tracts = compute_tract_lengths(row)
            print(f"\nIndividual {idx+1}")
            print("Loter ancestry calls (first 100):")
            print(row[:1000])
            print("Tract breakdown:")
            print(tracts[:10])  # Show first 10 tracts

            # Optional: reconstruct ancestry from tracts and compare
            reconstructed = []
            for ancestry, length in tracts:
                reconstructed.extend([ancestry] * length)
            match = np.array_equal(row, np.array(reconstructed))
            print(f"✅ Match with original: {match}")

    input_file_loter=f'{output_path}/{pop_to_run}_loter_output_phased.txt'
    match_loter_to_tracts(input_file_loter, num_individuals=3)



if __name__=='__main__':

    pop_to_run = sys.argv[1]
    print(pop_to_run)
    #print(f'Analyzing loter output for {pop_to_run}')
    #analyze_loter_tracts(pop_to_run)

    #print(f'Writing summary for {pop_to_run}')
    #tract_lengths.to_csv(output_file, index=False)

    if pop_to_run == 'vasquez' or pop_to_run == 'molino':
        end_of_eyeless = (19*2)-1
        eyeless_pop_index = [0,end_of_eyeless] # 0 to eyeless
        eyed_pop_index = [end_of_eyeless + 1, end_of_eyeless + (17*2) + 1]
        print(eyeless_pop_index)
    
    
    if pop_to_run == 'normal' or pop_to_run == 'othersurface':
        eyeless_pop_index = [0, 34]
        eyed_pop_index = [0, 34]

    print(eyeless_pop_index, eyed_pop_index)

    print(f'Making haplotype summaries for {pop_to_run}')
    haplo_stats(pop_to_run)


    print(f'Checking loter to summary for {pop_to_run}')
    test_output_match(pop_to_run)



    # Sad? Struggling? Feel like you've done this before but it's not here? 
    # Check out `old/inc.vasquez.summarystats.py` to see all the trouble shooting steps that got cut!
