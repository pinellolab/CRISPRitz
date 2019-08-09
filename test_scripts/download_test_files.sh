#!/bin/sh

#CREATE TEST DIRECTORY
mkdir test_crispritz
cd test_crispritz/

echo "DOWNLOADING ANNOTATIONS FILES"
#ANNOTATIONS FILE PATH AND ANNOTATIONS BED DIRECTORY DOWNLOAD
wget -q https://www.dropbox.com/s/ykw3dlxixhhs2rs/annotations_path.txt?dl=0
mv annotations_path.txt?dl=0 annotations_path.txt
wget -q https://www.dropbox.com/s/g5crgtmafd8ifw6/hg19_annotations.zip?dl=0
mv hg19_annotations.zip?dl=0 hg19_annotations.zip
unzip -q hg19_annotations.zip
rm hg19_annotations.zip
echo "ANNOTATIONS FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING GUIDES FILES"
#GUIDE DIRECTORY DOWNLOAD
wget -q https://www.dropbox.com/s/8th8yofj3l0q1gq/guides.zip?dl=0
mv guides.zip?dl=0 guides.zip
unzip -q guides.zip
rm guides.zip
echo "GUIDES FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING PAM FILES"
#PAM DIRECTORY DOWNLOAD
wget -q https://www.dropbox.com/s/z2d2s27ij3a7mky/pam.zip?dl=0
mv pam.zip?dl=0 pam.zip
unzip -q pam.zip
rm pam.zip
echo "PAM FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING CHR22 (TIME DEPENDING ON YOUR INTERNET CONNECTION)"
#CHR22 DOWNLOAD AND gunzip --quiet -q
mkdir hg19_ref
cd hg19_ref/
wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz
gunzip --quiet chr22.fa.gz
cd ..
echo "CHR22 DOWNLOADED AND EXTRACTED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING CHR22_VCF (TIME DEPENDING ON YOUR INTERNET CONNECTION)"
#CHR22 VCF DOWNLOAD
mkdir hg19_1000genomeproject_vcf
cd hg19_1000genomeproject_vcf/
wget -q http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
cd ..
echo "CHR22_VCF DOWNLOADED AND EXTRACTED"
echo "-------------------------------------------------------------------------"

echo "EVERYTHING DOWNLOADED AND READY TO USE"
