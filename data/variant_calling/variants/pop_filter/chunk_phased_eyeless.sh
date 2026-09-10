#!/usr/bin/bash

#53108128
#for i in {1..22} X Y MT
#do
#  bcftools view ${vcf_in} --regions ${i} -o ${vcf_out_stem}_${i}.vcf.gz -Oz
#done

input_vcf=Amex3.0_surface.phased.snps.eyeless.vcf.gz.vcf.gz
output_vcf_stem=./vcf_chunks/Amex3.0_surface.phased.snps.eyeless
chrs=list_of_chromosomes.txt

for i in $(cat list_of_chromosomes.txt)
do 
	bcftools view ${input_vcf} --regions ${i} -o ${output_vcf_stem}.${i}.vcf.gz -Oz
done
