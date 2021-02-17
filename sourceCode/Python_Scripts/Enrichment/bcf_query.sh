#!/usr/bin/env bash

#INPUT
#$1 vcf file
#$2 chr name

echo "Parsing: " $2
#query VCF file and exctracting position, ref and alt variant, AF and GT*SAMPLEs, filtering out any line with FILTER != PASS

bcftools view -h $1 | fgrep -v '##' | cut -f1,6,7,8 --complement | sed 's/FORMAT/AF/' > $2.alt
bcftools query -i 'FILTER="PASS"' -f '%POS\t%ID\t%REF\t%ALT\t%AF[\t%GT]\n' $1 | sed 's/\//|/g' >> $2.alt
