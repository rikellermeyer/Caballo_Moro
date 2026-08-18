#!/usr/bin/bash

#SBATCH --job-name=filter_indels_2             
#SBATCH --cpus-per-task=8             
#SBATCH --mem=16G             
#SBATCH --time="2-00:00"             
#SBATCH --error=./slurmout/filter_indels_2.%a.%A.err             
#SBATCH --output=./slurmout/filter_indels_2.%a.%A.out             
#SBATCH --mail-user=kell3262@umn.edu             
#SBATCH --mail-type=FAIL             
#SBATCH --mail-type=END

gatk SelectVariants            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/joint_call/jointcall.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/indels/subset.indels.vcf.gz            --select-type-to-include MIXED                              --select-type-to-include INDEL
gatk VariantFiltration            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/indels/subset.indels.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/indels/labeled.indels.vcf.gz            -filter "QD < 2.0" --filter-name "QD2"                                 -filter "QUAL < 30.0" --filter-name "QUAL30"                                 -filter "FS > 200.0" --filter-name "FS200"                                 -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
gatk SelectVariants            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/indels/labeled.indels.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/indels/filter_done.indels.vcf.gz            --exclude-filtered