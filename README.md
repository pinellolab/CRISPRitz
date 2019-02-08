# CRISPRITZ

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

**Call crispritz.py without input print the help, with a list of the callable methods**
```
- 'python3 crispritz.py'
help:
        crispritz add-variants <vcfFilesDirectory> <genomeDirectory>
        crispritz index-genome <name_genome> <genomeDirectory> <pamFile>
        crispritz search <genomeDirectory> <pamFile> <guideFile> <resultFile> {-db,-bf} (choose the type of search, db with index-genome, bf with fasta files) -mm <mm_num> [-bRNA <bRNA_num> | -bDNA <bDNA_num>] [-th <num_thread>] {-r,-p,-t} (write only off-targets results,write only profiles, write both)
        crispritz annotate-results <guidesFile> <resultsFile> <outputFile> -exons <exonsbedFile> -introns <intronsbedFile> -ctcf <ctcfbedFile> -dnase <dnasebedFile> -promoters <promotersbedFile>
        crispritz generate-report <guide> -mm <mm_num or range mm_min-mm_max> -profile <guideProfile> -extprofile <guideExtendedProfile> -exons <exonsCountFile> -introns <intronsCountFile> -ctcf <CTCFCountFile> -dnase <DNAseCountFile> -promoters <promotersCountFile> [-gecko (to use gecko pre-computed profile)]
```
**To create a variant version of the genome, you need to call add-variants and input the directory with vcf_files (also vcf.gz are accepted) and the directory of the reference genome (in .fasta files), the output will be saved in a new created folder called variants_genome, containing two copies of the inputted chromosome files, one with SNPs only and one with INDELs only**
```
- 'python3 crispritz.py add-variants vcf_files/ chroms_hg19/'
```
**To create a genome index, you need to call index-genome and input the name you choose for the directory of the result index genome (e.g, HG_19), the directory with the genome (in .fasta files) and a text file with the PAM you want to use (e.g.,NGG)**
```
- 'python3 crispritz.py index-genome HG_19 chroms_hg19/ pam/pamNGG.txt'
```

SOFTWARE REQUIREMENTS:
- PYTHON 3.*
- GCC 5.4 or later
- bcftools 1.2 or later
