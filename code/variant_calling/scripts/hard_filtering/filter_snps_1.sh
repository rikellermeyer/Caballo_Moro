#!/usr/bin/bash

#SBATCH --job-name=filter_snps_1             
#SBATCH --cpus-per-task=8             
#SBATCH --mem=16G             
#SBATCH --time="2-00:00"             
#SBATCH --error=./slurmout/filter_snps_1.%a.%A.err             
#SBATCH --output=./slurmout/filter_snps_1.%a.%A.out             
#SBATCH --mail-user=kell3262@umn.edu             
#SBATCH --mail-type=FAIL             
#SBATCH --mail-type=END

gatk SelectVariants            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/joint_call/jointcall.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/snps/subset.snps.vcf.gz            --select-type-to-include MNP                              --select-type-to-include SNP
gatk VariantFiltration            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/snps/subset.snps.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/snps/labeled.snps.vcf.gz            -filter "QD < 2.0" --filter-name "QD2"                                 -filter "QUAL < 30.0" --filter-name "QUAL30"                                 -filter "SOR > 3.0" --filter-name "SOR3"                                 -filter "FS > 60.0" --filter-name "FS60"                                 -filter "MQ < 40.0" --filter-name "MQ40"                                 -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5"                                 -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
gatk SelectVariants            -R /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/genome/Amex3.0_surface.fna            -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/snps/labeled.snps.vcf.gz            -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/variant_calling/variants/snps/filter_done.snps.vcf.gz            --exclude-filtered