#!/usr/bin/bash

echo "generating LD filter"
            plink --indep-pairwise 150 10 0.1             --vcf /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/nocall2.Amex3.0_surface.vcf.gz             --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/LD_filter_in.150win.Amex3.0_surface.pruned.150win            --double-id            --allow-extra-chr            --set-missing-var-ids @:#&&
            echo "making bed files"
            plink --pca            --vcf /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/nocall2.Amex3.0_surface.vcf.gz            --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/pruned.Amex3.0_surface.150.win            --extract /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/LD_filter_in.150win.Amex3.0_surface.pruned.150win.prune.in            --make-bed            --double-id            --allow-extra-chr            --set-missing-var-ids @:#&&
            echo "running admixture"
            for i in {1..5}; 
            do 
            ./admixture -j16 --cv /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/pruned.Amex3.0_surface.150.win.bed $i > /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/admixture/admix.Amex3.0_surface.150.win.log${i}.out 
            done