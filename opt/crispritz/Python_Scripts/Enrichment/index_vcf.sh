#!/usr/bin/env bash
# $1 is vcf file (Eg. variants_hg38/ALL.chr1.vcf.gz)
# $2 is chrname (Eg. chr1)
echo "Parsing: " $2
bcftools index $1 2>/dev/null
