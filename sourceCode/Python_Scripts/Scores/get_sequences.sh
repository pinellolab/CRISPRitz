#!/bin/bash
#awk -v len="23"  '/\+/ {print $4"\t"($5-4)"\t"($5+len+3)"\t"$1"_"$2"_"$5"_"$7"_"$8"\t""1""\t""+"} /\t\-\t/ {print $4"\t"($5-3)"\t"($5+len+4)"\t"$1"_"$2"_"$5"_"$7"_"$8"\t""1""\t""-"}' sourceCode/Python_Scripts/CFD/big.targets.txt.order.txt > big.targets.txt.bed 

#Da modificare se non voglio tenere tutti quei dati nel 4 campo

if [[ $4 == "True" ]]; then
awk -v len="$2"  '/\+/ {print $4"\t"($5-4)"\t"($5+len+3)"\t""forw""\t""1""\t""+"} 
/\t\-\t/ {print $4"\t"($5-3)"\t"($5+len+4)"\t""rev""\t""1""\t""-"}' $1 > $3.bed 
else
awk -v len="$2"  '/\+/ {print $4"\t"($5-3)"\t"($5+len+4)"\t""forw""\t""1""\t""+"} 
/\t\-\t/ {print $4"\t"($5-4)"\t"($5+len+3)"\t""rev""\t""1""\t""-"}' $1 > $3.bed 
fi