# CRISPRitz

CRISPRitz is a software package containing 5 different tools dedicate to perform analysis on CRISPR/Cas off-target analysis.
The aim of the software is to help the use performing tediuos and long tasks, such as, the search for off-targets on a genome, or the assessment of guides, in an easy and fast way.

![tool_overview](https://user-images.githubusercontent.com/32717860/52447053-faa3b200-2b2f-11e9-96fb-e3761a9232c5.png)


With this aim in mind we create this package, containing the following tools:

- add-variants: Created to permit the user to encode genomic variants using the IUPAC code, outputin a genome with added variant in IUPAC notation form.
- index-genome: Created to permit the user to find the potential targets in a genome based on the PAM and saving them in a set of .bin files, avoiding the necessity to perform this recursive task on every search with the same PAM
- search: The main function of the package, created to help the user perform off-target search on a genome(with variants) or genome index to help create safe and secure guides to use in a wet-lab CRISPR/Cas experiment.
- annotate-results: Created to permit the user to annotate regions founded during search with functional annotations (promoter, chromatin accessibility, insulator, etc)
- generate-report: Created to permit the user to generate a graphical report with annotated and overall mismatch and bulge profile for a given guide, to help the assessment of guides and to have an easy view over guides and their behave.

# EXAMPLE OF USAGE

**SHOW THE HELP**
```
- python3 crispritz.py
help:
        crispritz add-variants <vcfFilesDirectory> <genomeDirectory>
        crispritz index-genome <name_genome> <genomeDirectory> <pamFile>
        crispritz search <genomeDirectory> <pamFile> <guideFile> <resultFile> {-db,-bf} (choose the type of search, db with index-genome, bf with fasta files) -mm <mm_num> [-bRNA <bRNA_num> | -bDNA <bDNA_num>] [-th <num_thread>] {-r,-p,-t} (write only off-targets results,write only profiles, write both)
        crispritz annotate-results <guidesFile> <resultsFile> <outputFile> -exons <exonsbedFile> -introns <intronsbedFile> -ctcf <ctcfbedFile> -dnase <dnasebedFile> -promoters <promotersbedFile>
        crispritz generate-report <guide> -mm <mm_num or range mm_min-mm_max> -profile <guideProfile> -extprofile <guideExtendedProfile> -exons <exonsCountFile> -introns <intronsCountFile> -ctcf <CTCFCountFile> -dnase <DNAseCountFile> -promoters <promotersCountFile> [-gecko (to use gecko pre-computed profile)]
```

**CREATE A VARIANT GENOME**
```
- python3 crispritz.py add-variants vcf_files/ chroms_hg19/

**INPUT**
- vcf_files directory, the directory of the folder containing vcf files to process (.vcf or .vcf.gz allowed)
- reference_genome directory, the directory of the folder containing genome files (.fasta allowed)

**OUTPUT**
- Set of chromosomes with SNPs, copy of original chromosomes with variants SNPs
- Set of chromosomes with INDELs, copy of original chromosomes with variants INDELs
```

**CREATE A GENOME INDEX**
```
- python3 crispritz.py index-genome HG_19 chroms_hg19/ pam/pamNGG.txt

**INPUT**
- name of index genome, the name of the folder that will contain the index-genome created by the tool
- reference_genome directory, the directory of the folder containing genome files (.fasta allowed)
- PAM, text file containing the PAM (the PAM file must be written in this form NNNNNNNNNNNNNNNNNNNNNGG 3, the PAM must be long as future guides, in the example 20Ns, and the selected PAM, in the example NGG, the '3' means that the PAM sequence is represented by the last three letters.

**OUTPUT**
- index genome folder, folder containing the .bin files (one per input chromosome), every file contains all the candidate target founded on the chrosome for the selected PAM in input
```

**SEARCH ON A GENOME INDEX (SEARCH WITH MISMATCHES AND BULGES)**
```
- python3 crispritz.py search genome_library/NGG_HG_19/ pam/pamNGG.txt guides/1guides.txt result.ngg_hg19 -db -mm 4 -bDNA 1 -bRNA 1

**INPUT**
- genome directory, the directory of the folder containing genome files (.bin allowed)
- PAM, text file containing the PAM (the PAM file must be written in this form, NNNNNNNNNNNNNNNNNNNNNGG 3, the PAM must be long as future guides, in the example 20Ns, plus the selected PAM, in the example NGG, the '3' means that the PAM sequence is represented by the last three letters. REMEMBER THAT LENGTH OF INPUTTED GUIDES AND PAM MUST BE THE SAME.
- Guide/s, text file containing a list of guides (guides must be written in this form, TCACCCAGGCTGGAATACAGNNN, a sequence containing the actual guide you want to search, TCACCCAGGCTGGAATACAG, plus a number of Ns equal to the PAM length used to create the genome, in the exaple NNN. REMEMBER THAT LENGTH OF INPUTTED GUIDES AND PAM MUST BE THE SAME.
- name of result file, the name selected to nominate the output files (in the example 'result.ngg_hg19')
- db, flag set to use the search method constructed to search on index-genome
- mm #, flag set to select the mismatch threshold, in the example 4
- bDNA #, flag set to select the DNA bulge threshold, in the example 1
- bRNA, flag set to select the RNA bulge threshold, in the example 1

**OUTPUT**
- list of targets and off-targets, a text file containing all the targets and off-targets found on genome with genomic coordinates
- set of profile files, set of files containing detailed profile of every inputed guide, with count of bp/mismatches and a grouped representation of targets and off-targets founded per every mismatch count
```

**SEARCH ON A GENOME WITHOUT INDEXING (SEARCH WITH MISMATCHES ONLY)**
```
- python3 crispritz.py search enriched_genome/ pam/pamNGG.txt guides/1guides.txt result.hg19 -bf -mm 4 -t

**INPUT**
- directory, the directory of the folder containing genome files (.fasta allowed)
- PAM, text file containing the PAM (the PAM file must be written in this form, NNNNNNNNNNNNNNNNNNNNNGG 3, the PAM must be long as future guides, in the example 20Ns, plus the selected PAM, in the example NGG, the '3' means that the PAM sequence is represented by the last three letters.
- guide/s, text file containing a list of guides (guides must be written in this form, TCACCCAGGCTGGAATACAGNNN, a sequence containing the actual guide you want to search, TCACCCAGGCTGGAATACAG, plus a number of Ns equal to the PAM length used to create the genome, in the exaple NNN. REMEMBER THAT LENGTH OF INPUTTED GUIDES AND PAM MUST BE THE SAME.
- name of result file, the name selected to nominate the output files (in the example 'result.hg19')
- bf, flag set to use the search method constructed to search on genome without index
- mm #, flag set to select the mismatch threshold, in the example 4
- r or -p or -t, select which type of output you want, -r output only file containing list of targets and off-targets, -p output only files containing the profile(overall mismatch profile), -t output the complete set of files, so targets list and profile files (in the example '-t' selected)

**OUTPUT**
- list of targets and off-targets, a text file containing all the targets and off-targets found on genome with genomic coordinates
- set of profile files, set of files containing detailed profile of every inputed guide, with count of bp/mismatches and a grouped representation of targets and off-targets founded per every mismatch count
```

**ANNOTATE RESULTS**
```
- python3 crispritz.py index-genome HG_19 chroms_hg19/ pam/pamNGG.txt

**INPUT**
- name of index genome, the name of the folder that will contain the index-genome created by the tool
- reference_genome directory, the directory of the folder containing genome files (.fasta allowed)
- PAM, text file containing the PAM (the PAM file must be written in this form NNNNNNNNNNNNNNNNNNNNNGG 3, the PAM must be long as future guides, in the example 20Ns, and the selected PAM, in the example NGG, the '3' means that the PAM sequence is represented by the last three letters.

**OUTPUT**
- index genome folder, folder containing the .bin files (one per input chromosome), every file contains all the candidate target founded on the chrosome for the selected PAM in input
```








SOFTWARE REQUIREMENTS:
- PYTHON 3.*
- GCC 5.4 or later
- bcftools 1.2 or later
