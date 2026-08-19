#!/usr/bin/bash

#SBATCH --job-name=gwas_plink         
#SBATCH --array=0-0%1         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/gwas_plink.%A.out         
#SBATCH --error=./slurmout/gwas_plink.%A.err

echo "removing molino samples"
         bcftools view -s^M1,M2,M3,M4,M5,M6,M7         -o /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/nomolino.Amex3.0_surface.vcf.gz         --write-index        --force-samples /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/preplink/chr_corrected.Amex3.0_surface.vcf.gz&&
         echo "making LD files"
        plink         --vcf /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/nomolino.Amex3.0_surface.vcf.gz        --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/LD_filter.Amex3.0_surface        --double-id        --chr-set 25        --keep-allele-order        --indep-pairwise 50 10 0.1&&
         echo "making bed files"
         plink         --vcf /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/nomolino.Amex3.0_surface.vcf.gz         --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/gwas.Amex3.0_surface        --extract /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/LD_filter.Amex3.0_surface.prune.in        --make-bed        --double-id        --chr-set 25        --keep-allele-order