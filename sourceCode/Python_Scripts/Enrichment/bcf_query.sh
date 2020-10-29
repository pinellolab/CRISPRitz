#!/usr/bin/env bash

echo "Parsing: " $2
bcftools query -f '%POS %REF %AL T %QUAL\n' $1 > $2.alt
