#!/usr/bin/env bash

echo "Parsing: " $2
#query VCF file and exctracting position, ref and alt variant and quality, then filter out any line with QUAL != PASS
bcftools query -f '%POS %REF %ALT %QUAL\n' $1 | fgrep 'PASS' > $2.alt
