#!/usr/bin/bash

#SBATCH --job-name=batch_hard_filter         
#SBATCH --array=0-2%3         
#SBATCH --cpus-per-task=4         
#SBATCH --time="24-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=rk2643@stowers.org         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/batch_hard_filter.%A.out         
#SBATCH --error=./slurmout/batch_hard_filter.%A.err

sh filter*_${SLURM_ARRAY_TASK_ID}.sh