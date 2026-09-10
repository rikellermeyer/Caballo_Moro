#!/usr/bin/bash
#SBATCH --cpus-per-task=8 
#SBATCH --mem=32G 
#SBATCH --time="2-00:00" 

####each of these were run, just not all at once

#samtools faidx Amex3.0_surface.fna

#bwa-mem2 index Amex3.0_surface.fna

picard CreateSequenceDictionary R=Amex3.0_surface.fna O=Amex3.0_surface.dict
