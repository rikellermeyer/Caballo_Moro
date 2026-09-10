#!/usr/bin/bash

awk '
{
    chr=$1
    bp=$4

    if (!(chr in min) || bp < min[chr]) min[chr]=bp
    if (!(chr in max) || bp > max[chr]) max[chr]=bp
}
END {
    for (chr in min)
        print chr, min[chr], max[chr]
}' gwas.Amex3.0_surface.bim | sort -n
