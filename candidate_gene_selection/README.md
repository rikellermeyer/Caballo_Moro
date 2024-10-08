# Candidate gene selection

# Data for candidate gene comparisons
`ZF_DB_eye_genes.txt`: Database of zebrafish (D. rerio) eye genes pulled from Panther. 
`eye_qtls_from_Wiese2024.csv`: Eye QTLs were pulled from [Wiese et al., 2024](https://doi.org/10.1093/jhered/esae040)

# Annotating and comparing candidate genes
`gwas_to_gtf.py`: Annotating GWAS hits 
`qtl_comparison.py`: Takes genes identified by GWAS and SnpEff and compares them to known eye QTL regions.i
`qtl_matches_to_gtf.py`: Annotating QTL matches
