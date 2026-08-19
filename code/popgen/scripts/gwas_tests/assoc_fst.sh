#!/usr/bin/bash

#SBATCH --job-name=assoc_fst         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/assoc_fst.%A.out         
#SBATCH --error=./slurmout/assoc_fst.%A.err

echo "assoc_fst"
        plink         --file /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/gwas.Amex3.0_surface        --fst          --family        --allow-extra-chr        --allow-no-sex        --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/tests/assoc_fst