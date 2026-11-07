#!/usr/bin/bash

#SBATCH --job-name=pca_plink         
#SBATCH --array=0-3%4         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/pca_plink.%A.out         
#SBATCH --error=./slurmout/pca_plink.%A.err

echo "preprocessing pca"
         #sh pre_window_pca.sh&&
         echo "by-window plink commands"
         sh *_fullpca_${SLURM_ARRAY_TASK_ID}.sh
