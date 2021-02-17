#!/usr/bin/env bash

echo "Parsing: " $2
#query VCF file and exctracting position, ref and alt variant and GT*SAMPLEs, filtering out any line with FILTER != PASS
# bcftools query  -f '%POS %REF\t%ALT[\t%SAMPLES]\n' $1 > $2.alt
bcftools query -i 'FILTER="PASS"' -f '%POS\t%REF\t%ALT[\t%GT]\n' $1 > $2.alt