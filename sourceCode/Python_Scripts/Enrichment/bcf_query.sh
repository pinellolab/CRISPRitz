#!/usr/bin/env bash

echo "Parsing: " $2
bcftools query -f '%POS %REF %ALT %QUAL\n' $1 > $2.alt
