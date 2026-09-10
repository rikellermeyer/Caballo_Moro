#!/usr/bin/bash

list_of_files=Amex3.0_surface.*.eyeless.*.vcf.gz
output_file=all_eyeless_chunks_wc.txt

echo $list_of_files

for file in $list_of_files
do
	wc_output=$(bcftools view $file | grep -v '^#' | wc)
	echo -e "${file}\t${wc_output}" >> "$output_file"	
done

total_of_variants=$(awk '{s+=$2} END {print s}' $output_file)

echo -e "Total number of variants in all files: ${total_of_variants}"
