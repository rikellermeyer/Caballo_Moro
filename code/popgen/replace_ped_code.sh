#!/usr/bin/bash

ped_file="/projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/gwas.Amex3.0_surface.ped"

recode_file="/projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/newids.gwas.Amex3.0_surface.ped"

new_ped_file="/projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/gwas.recoded.Amex3.0_surface.ped"

awk '
BEGIN { OFS="\t" }
FNR==NR {
    c1[FNR] = $1
    c2[FNR] = $2
    c6[FNR] = $3
    next
}
{
    $1 = c1[FNR]
    $2 = c2[FNR]
    $6 = c6[FNR]
    print
}
' $recode_file $ped_file > $new_ped_file

mv $ped_file /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/badids_gwas.Amex3.0_surface.ped
mv $new_ped_file $ped_file