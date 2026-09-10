#!/usr/bin/bash

#SBATCH --mail-user=kell3262@umn.edu \
#SBATCH --cpus-per-task=4 \
#SBATCH --time=8:00:00 \
#SBATCH --mail-type=FAIL \
#SBATCH --mail-type=END\
#SBATCH --output=./slurmout/phased_wc.slurmout.%A.%a.out \
#SBATCH --error=./slurmout/phased_wc.slurmout.%A.%a.err
#SBATCH --mem=64G \

#53108128
#for i in {1..22} X Y MT
#do
#  bcftools view ${vcf_in} --regions ${i} -o ${vcf_out_stem}_${i}.vcf.gz -Oz
#done

input_vcf=Amex3.0_surface.phased.snps.eyed.vcf.gz.vcf.gz
output_vcf_stem=./vcf_chunks/Amex3.0_surface.phased.snps.eyed
chrs=list_of_chromosomes.txt

for i in $(cat list_of_chromosomes.txt)
do 
	bcftools view ${input_vcf} --regions ${i} -o ${output_vcf_stem}.${i}.vcf.gz -Oz
done
