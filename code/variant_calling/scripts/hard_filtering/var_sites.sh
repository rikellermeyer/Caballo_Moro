#!/usr/bin/bash

#SBATCH --job-name=var_sites             
#SBATCH --cpus-per-task=8             
#SBATCH --mem=16G             
#SBATCH --time="2-00:00"             
#SBATCH --error=./slurmout/var_sites.%a.%A.err             
#SBATCH --output=./slurmout/var_sites.%a.%A.out             
#SBATCH --mail-user=kell3262@umn.edu             
#SBATCH --mail-type=FAIL             
#SBATCH --mail-type=END

gatk IndexFeatureFile        -I /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/2025_Amex3.0_surface.all.filtered.vcf.gz
gatk SelectVariants         -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna        -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/2025_Amex3.0_surface.all.filtered.vcf.gz        --exclude-non-variants        -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/Amex3.0_surface.variants.filtered.vcf.gz