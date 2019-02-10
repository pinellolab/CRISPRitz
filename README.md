# CRISPRitz

CRISPRitz is a software package containing 5 different tools dedicated to perform predictive analysis on CRISPR/Cas experiments.
The aim of the software is to help the user, to perform tasks, such as, the search for off-targets on a genome, or the assessment of guides, in an easy and fast way.

![tool_overview](https://user-images.githubusercontent.com/32717860/52447053-faa3b200-2b2f-11e9-96fb-e3761a9232c5.png)


With this aim in mind we created this package, containing the following tools:

- add-variants: Tool to encode genomic variants using the IUPAC code, the ouput is a genome with added variant in IUPAC notation form.
- index-genome: Tool to find the potential targets in a genome based on the PAM and saving them in a set of .bin files, to avoid the need to perform this recursive task on every search with the same PAM.
- search: Tool created to perform off-target search on a genome (with or without variants) or genome index (with or without variants), to help design safe and secure guides to use in a wet-lab CRISPR/Cas experiment.
- annotate-results: Tool to annotate results found during search, with functional annotations (promoter, chromatin accessibility, insulator, etc)
- generate-report: Tool to generate a graphical report with annotated and overall mismatch and bulge profile for a given guide, to help the assessment of guides and to have an easy view over guides and their behave.

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
- set of chromosomes with SNPs, copy of original chromosomes with variants SNPs
- set of chromosomes with INDELs, copy of original chromosomes with variants INDELs
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
- python3 crispritz.py search genome_library/NGG_HG_19/ pam/pamNGG.txt guides/1guides.txt result.ngg_hg19 -index -mm 4 -bDNA 1 -bRNA 1

**INPUT**
- genome directory, the directory of the folder containing genome files (.bin allowed)
- PAM, text file containing the PAM (the PAM file must be written in this form, NNNNNNNNNNNNNNNNNNNNNGG 3, the PAM must be long as future guides, in the example 20Ns, plus the selected PAM, in the example NGG, the '3' means that the PAM sequence is represented by the last three letters. REMEMBER THAT LENGTH OF INPUTTED GUIDES AND PAM MUST BE THE SAME.
- Guide/s, text file containing a list of guides (guides must be written in this form, TCACCCAGGCTGGAATACAGNNN, a sequence containing the actual guide you want to search, TCACCCAGGCTGGAATACAG, plus a number of Ns equal to the PAM length used to create the genome, in the exaple NNN. REMEMBER THAT LENGTH OF INPUTTED GUIDES AND PAM MUST BE THE SAME.
- name of result file, the name selected to nominate the output files (in the example 'result.ngg_hg19')
- index, flag set to use the search method on index-genome
- mm #, flag set to select the mismatch threshold, in the example 4
- bDNA #, flag set to select the DNA bulge threshold, in the example 1
- bRNA, flag set to select the RNA bulge threshold, in the example 1

**OUTPUT**
- list of targets and off-targets, a text file containing all the targets and off-targets found on genome with genomic coordinates
- set of profile files, set of files containing detailed profile of every inputed guide, with count of bp/mismatches and a grouped representation of targets and off-targets founded per every mismatch count and bulges count
```

**SEARCH ON A GENOME WITHOUT INDEXING (SEARCH WITH MISMATCHES ONLY)**
```
- python3 crispritz.py search chroms/ pam/pamNGG.txt guides/1guides.txt result.hg19 -mm 4 -t

**INPUT**
- directory, the directory of the folder containing genome files (.fasta allowed)
- PAM, text file containing the PAM (the PAM file must be written in this form, NNNNNNNNNNNNNNNNNNNNNGG 3, the PAM must be long as future guides, in the example 20Ns, plus the selected PAM, in the example NGG, the '3' means that the PAM sequence is represented by the last three letters.
- guide/s, text file containing a list of guides (guides must be written in this form, TCACCCAGGCTGGAATACAGNNN, a sequence containing the actual guide you want to search, TCACCCAGGCTGGAATACAG, plus a number of Ns equal to the PAM length used to create the genome, in the exaple NNN. REMEMBER THAT LENGTH OF INPUTTED GUIDES AND PAM MUST BE THE SAME.
- name of result file, the name selected to nominate the output files (in the example 'result.hg19')
- mm #, flag set to select the mismatch threshold, in the example 4
- r or -p or -t, select which type of output you want, -r output only file containing list of targets and off-targets, -p output only files containing the profile(overall mismatch profile), -t output the complete set of files, so targets list and profile files (in the example '-t' selected)

**OUTPUT**
- list of targets and off-targets, a text file containing all the targets and off-targets found on genome with genomic coordinates
- set of profile files, set of files containing detailed profile of every inputed guide, with count of bp/mismatches and a grouped representation of targets and off-targets founded per every mismatch count
```

**ANNOTATE RESULTS**
```
- python3 crispritz.py annotate-results guides/1guides.txt result.hg19.targets.txt result.hg19.annotated -exons chroms_bed/hg19_exon.bed -introns chroms_bed/hg19_intron.bed -dnase chroms_bed/hg19_dnase.bed -ctcf chroms_bed/hg19_ctcf.bed -promoter chroms_bed/hg19_promoter.bed

**INPUT**
- guide/s, text file containing a list of guides used to create the result files
- list of targets and off-targets, text file created after search containing the list of targets and off-targets
- name of output file, chosen name of the output files (in the example result.hg19.annotated)
- exons, flag to input a .BED file containing genomic intervals annotated as exons
- introns, flag to input a .BED file containing genomic intervals annotated as introns
- dnase, flag to input a .BED file containing genomic intervals annotated as dnase
- ctcf, flag to input a .BED file containing genomic intervals annotated as ctcf
- promoter, flag to input a .BED file containing genomic intervals annotated as promoters

**OUTPUT**
- annotated list of targets and off-targets, text file containing annotated rows present in the input file (if a row is not annotated, is not written in the file)
- set of count files, set of files with counting information, one per genomic annotation, those files are formatted as a matrix, with guides as rows and mismatch as columns, every cell contain the number of annotated targets for the guide with that mismatch
```

**GENERATE REPORT**
```
- python3 crispritz.py generate-report AGCCGGGACGTCAGCACTACNNN -mm 3-5 -profile results/pcsk9/reference/pcsk9.reference.profile.xls -extprofile results/pcsk9/reference/pcsk9.reference.extended_profile.xls -exons results/pcsk9/reference/pcsk9.reference.annotated.ExonsCount.txt -introns results/pcsk9/reference/pcsk9.reference.annotated.IntronsCount.txt -ctcf results/pcsk9/reference/pcsk9.reference.annotated.CTCFCount.txt -dnase results/pcsk9/reference/pcsk9.reference.annotated.DNAseCount.txt -promoters results/pcsk9/reference/pcsk9.reference.annotated.PromotersCount.txt -gecko -sumone results/pcsk9/reference/pcsk9.reference.annotated.SummaryCount.txt -sumtwo results/pcsk9/enriched/pcsk9.enriched.annotated.SummaryCount.txt

**INPUT**
- guide, guide present in the profile and extended profile file
- profile, flag to input a profile file generated by the search tool, containing information about behave of every guide previously analyzed
- extprofile, flag to input an extended profile file generated by the search tool, containing matrix logo of every guide previously analyzed
- exons, flag to input a exons.count file, generated in a previously annotation phase
- introns, flag to input a introns.count file, generated in a previously annotation phase
- dnase, flag to input a dnase.count file, generated in a previously annotation phase
- ctcf, flag to input a ctcf.count file, generated in a previously annotation phase
- promoter, flag to input a promoter.count file, generated in a previously annotation phase
- gecko, flag set to perform a comparison with the gecko pre-computed profile files, useful to understand behave of selected guide with respect to a well-known guide library [optional]
- sumone, flag to input a reference summary count, generated in a previously annotation phase [optional]
- sumtwo, flag to input a enriched summary count, generated in a previously annotation phase [optional]

**OUTPUT**
- detailed graphical report for inputted guide, an image to visualize how the guide behave in terms of targets founded in each genomic annotation and in general (total targets count)
- barplot with percentage increasing from Variant/Reference Genome, a barplot to visualize the percentage gain terms of founded targets in total, and grouped by every genomic annotation used in the annotation phase [optional, generated only if sumone and sumtwo inputted]
```

SOFTWARE REQUIREMENTS:
- PYTHON 3.*
- GCC 5.4 or later
- bcftools 1.2 or later
- intervaltree-3.0.2 or later
