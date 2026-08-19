#!/usr/bin/bash

#SBATCH --job-name=preplink         
#SBATCH --array=0-0%1         
#SBATCH --cpus-per-task=4         
#SBATCH --time="2-00:00"         
#SBATCH --mem=32G         
#SBATCH --mail-user=kell3262@umn.edu         
#SBATCH --mail-type=FAIL         
#SBATCH --mail-type=END        
#SBATCH --output=./slurmout/preplink.%A.out         
#SBATCH --error=./slurmout/preplink.%A.err

echo "correcting vcf"
bcftools view         -o /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/preplink/chr_corrected.Amex3.0_surface.vcf.gz         --force-samples /projects/standard/mcgaughs/shared/kell3262/CaballoMoroOfficial/Caballo_Moro/data/popgen/preplink/Amex3.0_surface.renamed_chr.vcf.gz        --write-index        --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25