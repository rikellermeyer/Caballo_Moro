#!/usr/bin/bash

#SBATCH --job-name=makeped         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/makeped.%A.out         
#SBATCH --error=./slurmout/makeped.%A.err

plink         --vcf /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/nomolino.Amex3.0_surface.vcf.gz        --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/gwas.Amex3.0_surface        --maf 0.05        --allow-extra-chr        --chr 1-25        --recode        --set-missing-var-ids @:#$1:$2