#!/usr/bin/bash

#SBATCH --job-name=all_tests         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/all_tests.%A.out         
#SBATCH --error=./slurmout/all_tests.%A.err

sh assoc.sh&&
 sh CA_trend.sh &&
 sh assoc_fst.sh