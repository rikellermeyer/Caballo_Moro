# Local ancestry inference with Loter

## `pre_loter_filter.py`
- Takes in the original, all sites vcf file
- Splits by population and filters to snps
- Makes beagle scripts for phasing (Phasing done by Beagle 5.4)

## `loter_ancestry.py` runs Loterv1.0
- Takes in phased vcfs that are split by population
- From `variant_calling/variants/Amex3.0_surface.phased.{population}.vcf.gz.vcf.gz
- Loter outputs a matrix that has assigned SNPs to eyeless (0) or surface (1) ancestry
- This output is split by haplotypes (because phased), with two rows per individual
- Plots loter output as a heatmap

## `loter_summarystats.py` is a custom script to analyze loter output
- Parses the loter output matrix
- Calculates the tract lengths of 0's (eyeless) and 1's (surface)
- Then averages across haplotypes (not individuals!)
- And calculates Tadmix
- This script also has the depreciated code for running (collapse_ancestry)[https://github.com/armartin]

### For simpler analysis
See: `old/loter_ancestry.py` and `old/loter_output_analysis.py`