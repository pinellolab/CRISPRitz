# $1 is genome directory
# $2 is bedfile.bed
mkdir tmp_seq
for chr in $1/*.fa; do
    chr_name=$(basename $chr)
    bedtools getfasta -fi $chr -bed $2 -fo tmp_seq/$chr_name'_seq.txt' 2>/dev/null
done
find ./tmp_seq -type f -empty -delete

for extr in tmp_seq/*_seq.txt; do
    sed -i '1~2d' $extr
done