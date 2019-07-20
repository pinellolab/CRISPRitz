#!/usr/bin/env bash
VCFFILES=$1*
for f in $VCFFILES
do
	IFS='.' read -r -a array <<< "$f"
	echo "Parsing: " ${array[1]}
	bcftools query -f '%POS %REF %ALT %QUAL\n' $f > ${array[1]}.alt
done
