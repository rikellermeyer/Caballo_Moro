#!/usr/bin/bash

#SBATCH --job-name=CA_trend         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/CA_trend.%A.out         
#SBATCH --error=./slurmout/CA_trend.%A.err

echo "CA_trend"
        plink         --file /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/gwas.Amex3.0_surface        --model          --allow-extra-chr        --allow-no-sex        --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/tests/CA_trend