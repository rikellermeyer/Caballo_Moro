#!/usr/bin/bash

#SBATCH --job-name=rename_vcf         
#SBATCH --array=0-0%1         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/rename_vcf.%A.out         
#SBATCH --error=./slurmout/rename_vcf.%A.err

bcftools annotate         --rename-chrs /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/Amex3.0_surface.renamechr.txt         --write-index        -o /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/preplink/Amex3.0_surface.renamed_chr.vcf.gz /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/2025_Amex3.0_surface.variants.filtered.vcf.gz