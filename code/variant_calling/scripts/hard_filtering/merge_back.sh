#!/usr/bin/bash

#SBATCH --job-name=merge_back             
#SBATCH --cpus-per-task=8             
#SBATCH --mem=16G             
#SBATCH --time="2-00:00"             
#SBATCH --error=./slurmout/merge_back.%a.%A.err             
#SBATCH --output=./slurmout/merge_back.%a.%A.out             
#SBATCH --mail-user=kell3262@umn.edu             
#SBATCH --mail-type=FAIL             
#SBATCH --mail-type=END


gatk MergeVcfs        -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna        -I /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/snps/filter_done.snps.vcf.gz -I /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/indels/filter_done.indels.vcf.gz        -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/Amex3.0_surface.biallelic.filtered.vcf.gz&&
gatk MergeVcfs        -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna        -I /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/Amex3.0_surface.biallelic.filtered.vcf.gz        -I /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/mono/filter_done.mono.vcf.gz        -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/Amex3.0_surface.all.filtered.vcf.gz