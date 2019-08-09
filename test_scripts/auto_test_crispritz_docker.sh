#!/bin/sh

#CREATE TEST DIRECTORY
mkdir test_crispritz
cd test_crispritz/

echo "DOWNLOADING ANNOTATIONS FILES"
#ANNOTATIONS FILE PATH AND ANNOTATIONS BED DIRECTORY DOWNLOAD
wget -q https://www.dropbox.com/s/cri6jqveh85u3do/annotations_path.txt?dl=0
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

echo ""
echo "STARTING THE AUTO TEST SCRIPT"
echo "THIS SCRIPT WILL NOW TEST ALL THE CRISPRitz FUNCTION TO CHECK THE INSTALLATION"
echo ""
echo "TESTING ADD-VARIANTS (IF THE DOCKER IS NOT AVAILABLE ON YOUR SYSTEM, WILL BE DOWNLOADED AND INSTALLED AUTOMATICALLY)"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py add-variants hg19_1000genomeproject_vcf/ hg19_ref/ > output.redirect.out
echo -e "ADD-VARIANTS \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING INDEX-GENOME"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py index-genome hg19_ref hg19_ref/ pam/pamNGG.txt -bMax 2 > output.redirect.out
echo -e "INDEX-GENOME \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING SEARCH WITH ONLY MISMATCHES"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search hg19_ref/ pam/pamNGG.txt guides/EMX1.txt emx1.hg19 -mm 4 -t -scores hg19_ref/ > output.redirect.out
echo -e "SEARCH WITH ONLY MISMATCHES \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING SEARCH WITH MISMATCHES AND BULGES"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search genome_library/NGG_2_hg19_ref/ pam/pamNGG.txt guides/EMX1.txt emx1.hg19.bulges -index -mm 4 -bDNA 1 -bRNA 1 -t > output.redirect.out
echo -e "SEARCH WITH MISMATCHES AND BULGES \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING ANNOTATE-RESULTS"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py annotate-results guides/EMX1.txt emx1.hg19.targets.txt annotations_path.txt emx1.hg19.annotated > output.redirect.out
echo -e "ANNOTATE-RESULTS \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING GENERATE-REPORT"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py generate-report GAGTCCGAGCAGAAGAAGAANNN -mm 4 -profile emx1.hg19.profile.xls -extprofile emx1.hg19.extended_profile.xls -exons emx1.hg19.annotated.ExonsCount.txt -introns emx1.hg19.annotated.IntronsCount.txt -dnase emx1.hg19.annotated.DNAseCount.txt -ctcf emx1.hg19.annotated.CTCFCount.txt -promoters emx1.hg19.annotated.PromotersCount.txt -gecko > output.redirect.out
echo -e "GENERATE-REPORT \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

rm output.redirect.out
echo ""
echo -e "EVERY TEST \e[32mPASSED\e[0m!!! ENJOY CRISPRitz"