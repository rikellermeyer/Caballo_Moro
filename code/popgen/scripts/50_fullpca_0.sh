#!/usr/bin/bash

echo "generating LD filter"
            plink --indep-pairwise 50 10 0.1             --vcf /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/nocall2.Amex3.0_surface.vcf.gz             --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/LD_filter_in.50win.Amex3.0_surface.pruned.50win            --double-id            --allow-extra-chr            --set-missing-var-ids @:#&&
            echo "making bed files"
            plink --pca            --vcf /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/nocall2.Amex3.0_surface.vcf.gz            --out /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/pruned.Amex3.0_surface.50.win            --extract /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/LD_filter_in.50win.Amex3.0_surface.pruned.50win.prune.in            --make-bed            --double-id            --allow-extra-chr            --set-missing-var-ids @:#&&
            echo "running admixture"
            for i in {1..5}; 
            do 
            ./admixture -j16 --cv /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/pruned.Amex3.0_surface.50.win.bed $i > /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/pca/admixture/admix.Amex3.0_surface.50.win.log${i}.out 
            done