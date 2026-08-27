#!/usr/bin/bash

input_file="/projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/gwas/tests/CA_trend.model"
output_file="/projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/candidate_genes/CA_trend.filtered.tsv"

awk 'NR==1 || ($5=="TREND" && $10!="NA" && $10<=1e-8)' $input_file > $output_file