#!/usr/bin/bash

#SBATCH --job-name=filter_mono_0             
#SBATCH --cpus-per-task=8             
#SBATCH --mem=16G             
#SBATCH --time="2-00:00"             
#SBATCH --error=./slurmout/filter_mono_0.%a.%A.err             
#SBATCH --output=./slurmout/filter_mono_0.%a.%A.out             
#SBATCH --mail-user=kell3262@umn.edu             
#SBATCH --mail-type=FAIL             
#SBATCH --mail-type=END

gatk SelectVariants            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/joint_call/jointcall.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/mono/subset.mono.vcf.gz            --select-type-to-exclude INDEL                             --select-type-to-exclude MIXED                             --select-type-to-exclude MNP                             --select-type-to-exclude SYMBOLIC                             --select-type-to-include NO_VARIATION
gatk VariantFiltration            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/mono/subset.mono.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/mono/labeled.mono.vcf.gz            -filter "QD < 2.0" --filter-name "QD2"                                 -filter "FS > 60.0" --filter-name "FS60"                                 -filter "MQ < 40.0" --filter-name "MQ40" 
gatk SelectVariants            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/mono/labeled.mono.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/mono/filter_done.mono.vcf.gz            --exclude-filtered