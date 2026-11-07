#!/usr/bin/bash

echo "Indexing vcf for gatk"
        gatk IndexFeatureFile         -I /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/preplink/chr_corrected.Amex3.0_surface.vcf.gz&&
        echo "filtering to nocalls"
        gatk SelectVariants        -V /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/preplink/chr_corrected.Amex3.0_surface.vcf.gz        --max-nocall-number 2        --exclude-sample-name M2        -O /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/nocall2.Amex3.0_surface.vcf.gz