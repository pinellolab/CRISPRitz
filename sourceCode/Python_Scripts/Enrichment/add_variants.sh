#!/usr/bin/env bash
# $1 is vcf file (Eg. variants_hg38/ALL.chr1.vcf.gz)
# $2 is chrname with file extension(Eg. chr1.fa)
# $3 is chrname (Eg. chr1)

bcftools consensus -I -M N -f $2 $1 > $3'.enriched.fa'
