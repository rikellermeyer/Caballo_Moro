# Visualizing Caballo Moro Data

# CMC_graphics.Rmd
- Chunk: `popgenstats` is unused but informative Fst + Dxy plotting along chromosomes
- Chunk: `pca_dataframe` and `plot_pca`:
    1. Preps pca files
    2. Does a normal PCA plot (Figure 1)
    3. Calculates Euclidean distance, kmeans clustering, etc for supplemental figures
- Chunk: `load_admixture` and `plot_admixture`:
    1. Preps admixture files
    2. Makes admixture plots
    Note: Cross-validation error is calculated with Python in `../code/popgen/preplink.py`
- Various ways to plot GWAS. The best though is just from `CMC_GWAS.Rmd` and remove the QTL if you want only GWAS

# CMC_GWAS.Rmd
This is a later stage R markdown file that makes graphics for:
- GWAS + QTL at the bottom - chunk: `CA_and_QTL`, graphic: `CA_trend_and_QTL.svg`
- A plot for each chromosome with GWAS + QTL + regions with selective sweeps - chunk: `Chrom_plots_sep`, graphics: `chr_gwas_qtl_plots/updown_Chr_*`