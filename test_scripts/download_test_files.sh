#!/bin/sh

#CREATE TEST DIRECTORY
mkdir test_crispritz
cd test_crispritz/

echo "DOWNLOADING ANNOTATIONS FILES"
#ANNOTATIONS FILE PATH AND ANNOTATIONS BED DIRECTORY DOWNLOAD
curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/hg38Annotation.zip --output hg38Annotation.zip
unzip hg38Annotation.zip
echo "ANNOTATIONS FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING GUIDES FILES"
#GUIDE DIRECTORY DOWNLOAD
curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/EMX1.sgRNA.txt --output EMX1.sgRNA.txt
echo "GUIDES FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING PAM FILES"
#PAM DIRECTORY DOWNLOAD
curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/20bp-NGG-SpCas9.txt --output 20bp-NGG-SpCas9.txt
echo "PAM FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING CHR22 (TIME DEPENDING ON YOUR INTERNET CONNECTION)"
#CHR22 DOWNLOAD AND gunzip --quiet -q
mkdir hg38_ref
cd hg38_ref/
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz --output chr22.fa.gz
gunzip --quiet chr22.fa.gz
cd ..
echo "CHR22 DOWNLOADED AND EXTRACTED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING CHR22_VCF (TIME DEPENDING ON YOUR INTERNET CONNECTION)"
#CHR22 VCF DOWNLOAD
mkdir hg38_1000genomeproject_vcf
cd hg38_1000genomeproject_vcf/
curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz --output ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
cd ..
echo "CHR22_VCF DOWNLOADED AND EXTRACTED"
echo "-------------------------------------------------------------------------"

echo "EVERYTHING DOWNLOADED AND READY TO USE"
