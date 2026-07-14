#!/usr/bin/bash

#SBATCH --mail-user=rk2643@stowers.org \
#SBATCH --cpus-per-task=4 \
#SBATCH --mail-type=FAIL \
#SBATCH --mail-type=END\
#SBATCH --output=./slurmout/beagle.%A.out \
#SBATCH --error=./slurmout/beagle.%A.err
#SBATCH --mem=64G \

#use conda loter
#beagle5.1
#Amex3.0_surface.filtered.eyed.vcf.gz

populations=(eyed eyeless surface)

for population in "${populations[@]}" ; do
    echo "$population"
    java -jar -Djava.io.tmpdir=./temp/ -Xmx32g beagle.27May24.118.jar nthreads=16 gt=/n/projects/rk2643/caballo_moro_genomics/streamlined/data/variant_calling/variants/geno_filter/Amex3.0_surface.filtered.$population.vcf.gz out=/n/projects/rk2643/caballo_moro_genomics/streamlined/data/variant_calling/variants/Amex3.0_surface.phased.$population.vcf.gz
done