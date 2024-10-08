Popgen analysis
1. `preplink.py`: takes filtered, biallelic (SNPs and INDELs) vcf and runs plink: generates a PCA and Admixture 
2. `gwas_fixandtest.py`: Runs a series of GWAS trend tests with PLINKv1.90b6.21. Visualized with `../reports/CMC_graphics.rmd` 
3. `popgen_windows.py`: Does the genetic architecture compuations (Dxy, pi, etc). Uses smartin_genomics/.
4. `loter/`: does everything needed for local ancestry inference. Has it's own README!
