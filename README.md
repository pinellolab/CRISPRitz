# CRISPRitz 
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/crispritz/README.html)

CRISPRitz is a software package containing 5 different tools dedicated to perform predictive analysis and result assessement on CRISPR/Cas experiments.

![fig_pipeline_breve-1](https://user-images.githubusercontent.com/32717860/53018798-93e99700-3453-11e9-85a9-89eba18b18bb.png)

With this aim in mind we created this package, containing the following software tools:

- add-variants: Tool to encode genomic variants using the IUPAC code. The ouput is a set of files representing the genome enriched with variants in IUPAC notation form.
- index-genome: Tool to find the candidate targets in a genome starting from a PAM. The ouput is a set of files, containing all the sequences of candidate targets extracted from the genome.
- search: Tool to perform off-target search on a genome (with or without variants) or genome index (with or without variants). The ouput is a set of files, one is a list of all targets and off-targets found, the others are profile files containing detailed information for each guide , like bp/mismatches and on/off-targets count.
- annotate-results: Tool to annotate results found during search with functional annotations (promoter, chromatin accessibility, insulator, etc). The output is a set of files, one is the list of on/off-targets with the annotation type, the others are files containing counts for each guide, the counts are the total on/off-targets found with the specific mismatch threshold and the specific annotation.
- generate-report: Tool to generate a graphical report with annotated and overall mismatch and bulge profile for a given guide. The output is a graphical representation of the input guide behaviour.
- process-data: Tool to process and compare the off-target results done on the reference and enriched genomes. Informations about samples from the 1000 genome project are also provided for each target, giving a multi-level insight on the guide behaviour at general, population and sample level.

# CRISPRitz Installation and Usage
The two fastest way to use CRISPRitz is through the installation of Docker or Conda.
Here we summarize the steps to install CRISPRitz with Docker and Conda.

## Installation (Phase 1)
**Conda installation (Linux and MacOS):**
- Open a terminal window
- Paste this command into the terminal (Linux):
    ```
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --output Miniconda3-latest-Linux-x86_64.sh
    ```
- Paste this command into the terminal (MacOS):
    ```
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh --output Miniconda3-latest-MacOSX-x86_64.sh
    ```
- If the file is correctly downloaded you now need to execute it to complete the installation, so paste this command into the terminal:
    - Linux
    ```
    bash Miniconda3-latest-Linux-x86_64.sh
    ```
    - MacOS
    ```    
    bash Miniconda3-latest-MacOSX-x86_64.sh
    ```
- Press ENTER when requested and yes when an answer is requested, in this way you allow conda to set all the directories in your HOME path for an easy use
- After the complete installation you will receive this message “Thank you for installing Miniconda3!” to certify the correct installation.
- Now you need to close the terminal window you are using and open a new one, to allow the system to start conda.
- In the new terminal window you should see something like this:
    ```
    (base) user@nameofPC:~$
    ```
    If you read the "(base)" like this, conda is loaded correctly and you can start using it.
- Now you need to set the channels to allow conda to access different repositories and set the default version of python to version 3.6, so paste these commands into the terminal you just opened:
    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install python=3.6
    ```
- Now, you can install CRISPRitz by typing the command:
    ```
    conda install crispritz
    ```
- To test your installation, type the command:
    ```
    crispritz.py
    ```
- After the execution of the command you should see a list of CRISPRitz tools.
![crispritz.py_help](https://user-images.githubusercontent.com/40895152/63214203-8452be80-c115-11e9-88e2-4613ba8c3718.png)
Now the software is installed and ready to be used.

**Docker installation:  
Note: if you are using MasOS or Windows, you just need to download the installer file
and follow the on screen instructions.  
https://docs.docker.com/docker-for-windows/install/ (Windows)  
https://docs.docker.com/docker-for-mac/install/ (MacOS)**

**Ubuntu installation guide:**
- Open a terminal window
- Paste this command to update the index repository:
    ```
    sudo apt-get update
    ```
- Paste this command to allow package installation over HTTPS:
    ```
    sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
    ```
- Paste this command to add the docker key:
    ```
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    ```
- Paste this command to set the correct version of docker for your system:
    ```
    sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    Stable"
    ```
- Paste this command to update the index repository another time, to make sure everything is ready and set to install docker:
    ```
    sudo apt-get update
    ```
- Then paste this command to finally install docker:
    ```
    sudo apt-get install docker-ce docker-ce-cli containerd.io
    ```
- Paste this last command to check if the installation is complete and functional:
    ```
    sudo docker run hello-world
    ```
- If this message is printed, everything is perfectly installed
![docker hello world](https://user-images.githubusercontent.com/40895152/63214349-769e3880-c117-11e9-8ee2-d754096b3aca.png)
- Now, we need to do some more steps to complete the settings. Paste this command to create a user group for docker user:
    ```
    sudo groupadd docker
    ```
- Paste this command to add your current user to the created group:
    ```
    sudo usermod -aG docker $USER
    ```
- Now you need to restart your machine or the virtual environment, to re-evaluate the user groups.
- One last command to test if the group is well configured. Paste this command:
    ```
    docker run hello-world
    ```
- If the previous “hello from docker” message is printed, everything is perfectly set.

## Post installation test (Phase <a name="phase2">2</a>):
**Conda:**
- Download and run this script if you have installed CRISPRitz with Conda:
    ```
    curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/auto_test_crispritz_conda.sh --output auto_test_crispritz_conda.sh
    ```
- Write this command to execute the script:
    ```
    bash auto_test_crispritz_conda.sh
    ```
- Wait until this confirmation message appears:  
“EVERY TEST PASSED!!! ENJOY CRISPRitz”

**Docker:**
- Download and run this script if you have installed CRISPRitz with Docker:
    ```
    curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/auto_test_crispritz_docker.sh --output auto_test_crispritz_docker.sh
    ```
- Write this command to execute the script:
    ```
    bash auto_test_crispritz_docker.sh
    ```
- Wait until this confirmation message appears:  
“EVERY TEST PASSED!!! ENJOY CRISPRitz”

## Usage (Phase 3):
Here is a brief guide to help use CRISPRitz, **if you already execute the post installation test
(Phase [2](#phase2)), and you obtain a positive result, you have all the necessary file in the
test_crispritz directory and you can skip this list of steps.**  
If you did not execute the test, follow these few steps to download the necessary files to try
CRISPRitz.  
Download test files (ONE TIME STEP):
- The script will download the chr22 from UCSC (hg19), the correspondent VCF file from
the 1000 Genome Project, a directory containing some test guides, a directory
containing some PAM sequences and a directory of pre-computed genomic annotations
for the hg19 genome.
- Download the script with this command:
    ```
    curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/download_test_files.sh --output download_test_files.sh
    ```
- Write this command to execute the script:
    ```
    bash download_test_files.sh
    ```
- The script will download every necessary file to test the software, we download only one
chromosome and one vcf file, to save time. All the examples can be run on an entire
genome, if you want to use the entire hg19 genome, you only need to add chromosomes
into the `hg19_ref` directory.
- Write this command to enter the test directory:
    ```
    cd test_crispritz/
    ```
- Now you are ready to execute the following example functions.

**<a name="Add-Variant">3.1</a> CRISPRitz Add-Variant Tool**  
This tool is created to insert variants in a fasta genome.  
Input:
- Directory containing a genome in fasta format, need to be separated into single
chromosome files.
- Directory containing VCF files, need to be separated into single chromosome files
(multi-sample files will be collapsed into one fake individual).

Output:
- Directory containing a duplicate of the original genome in fasta format, separated into
single chromosome files with added SNPs in IUPAC notation
- Directory containing a duplicate of the original genome in fasta format, separated into
single chromosome files with added INDELs.

Example call:
- Conda
    ```
    crispritz.py add-variants hg19_1000genomeproject_vcf/ hg19_ref/
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py add-variants hg19_1000genomeproject_vcf/ hg19_ref/
    ```

Detailed input:  
`hg19_1000genomeproject_vcf/` is the directory containing the vcf files.  
`hg19_ref/` is the directory containing the fasta files.


**<a name="Index-Genome">3.2</a> CRISPRitz Index-Genome Tool**  
This tool is created to generate an index genome (similar to the bwa-index step). This step is
time consuming (from 30 to 60 minutes) but helps to save a lot of execution time while
searching with lot of guides and with the support of bulges (DNA and RNA). If do not need to
search with bulges, skip this passage.  
Input:
- Name of the genome to create (e.g. `hg19_ref`).
- Directory containing a genome in fasta format, need to be separated into single
chromosome files.
- Text file containing the PAM (including a number of Ns equal to the guide length) and a
space separated number indicating the length of the PAM sequence (e.g. Cas9 PAM is
NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by 20 Ns and
NGG, followed by number 3, representing the length of the PAM sequence.
- Number of bulges to include in the database to perform the following search (i.e. the max
number bulges allowed for DNA and RNA when searching on the database)
- Number of threads to use for the analysis (Optional)

Output:
- Directory containing an index genome in .bin format, separated into single chromosome
files, containing all the candidate targets for a selected PAM, adding also characters to
perform bulge search.

Example call:
- Conda
    ```
    crispritz.py index-genome hg19_ref hg19_ref/ pam/pamNGG.txt -bMax 2 -th 4
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py index-genome hg19_ref hg19_ref/ pam/pamNGG.txt -bMax 2 -th 4
    ```
Detailed input:  
`hg19_ref` is the name of the output directory containing the index genome.  
`hg19_ref/` is the directory containing the fasta files.  
`pam/pamNGG.txt/` is a text file containing the PAM sequence.  
`-bMax 2` is the max number of bulges to allow for following searches when creating the
database (e.g if -bMax 2, all the following searches with the created index can be performed
with max 2 RNA bulges and 2 DNA bulges)  
`-th 4` is the number of threads to use (Optional)


**<a name="Search">3.3</a> CRISPRitz Search Tool**  
This tool is created to search on a fasta genome or an index genome.
There are two kinds of searches permitted with CRISPRitz;  
The first and simplest one uses a common fasta genome and it’s developed to perform fast,
on-the-fly searches with mismatches only.  
The second search type, uses the before generated index genome (Phase [3.2](#Index-Genome)), to perform
searches with lot of guides and also with bulges support.  

**<a name="Search_mm">3.3.1</a> Mismatches only search:**  
Input:
- Directory containing a genome in fasta format, need to be separated into single
chromosome files.
- Text file containing the PAM sequence (including a number of Ns equal to the guide
length) and a space separated number indicating the length of the PAM sequence (e.g.
Cas9 PAM is NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by
20 Ns and NGG, followed by 3, representing the length of the PAM sequence.
- Text file containing one or more guides (including a number of Ns equal to the length of
the PAM sequence) (e.g. TCACCCAGGCTGGAATACAGNNN, the last 3 Ns represents
the space occupied by the PAM in the real sequence)
- Name of the output file (e.g. `emx1.hg19`)
- Number of allowed mismatches (e.g. `-mm 4`)
- Output type (-r off-targets list only, -p profile only, -t everything) (e.g `-t`)
- Scores (-scores followed by the directory of the fasta genome, to perform the score after
the search with score calculation based on Doench 2016 and CFD, the two scoring
methods supports only NGG PAM and 23 long guides) (e.g `-scores hg19_ref/`)
- Number of threads to use for the analysis (Optional)

Output:
- Set of result files, including:
    - Targets file, containing all genomic targets for the guide set
    - Profile file, containing a matrix-like representation of guides behaviour (bp/mm, total on-/off- target, targets per mismatch threshold)
    - Extended profile file, containing the motif matrix for each guide and each mismatch threshold, useful to create visual analysis of the guides behaviour
    - Targets file with associated CFD score

Example call:
- Conda
    ```
    crispritz.py search hg19_ref/ pam/pamNGG.txt guides/EMX1.txt emx1.hg19 -mm 4 -t -scores hg19_ref/
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search hg19_ref/ pam/pamNGG.txt guides/EMX1.txt emx1.hg19 -mm 4 -t -scores hg19_ref/
    ```

Detailed input:  
`hg19_ref/` is the directory containing the fasta files.  
`pam/pamNGG.txt/` is a text file containing the PAM sequence.  
`guides/EMX1.txt` is a text file containing the EXM1 guide  
`emx1.hg19` is the output file name  
`-mm 4` to select the mismatch threshold  
`-t` to select the output type  
`-scores hg19_ref/` to activate the calculation of score (Doench 2016 and CFD)  
`-th 4` is the number of threads to use (Optional)

**<a name="Search_mm_bul">3.3.2</a> Mismatches + Bulges search:**  
Input:
- Directory containing an index genome in .bin format, separated into single chromosome
files (Phase [3.2](#Index-Genome)).
- Text file containing the PAM sequence (including a number of Ns equal to the guide
length) and a space separated number indicating the length of the PAM sequence (e.g.
Cas9 PAM is NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by
20 Ns and NGG, followed by 3, representing the length of the PAM sequence.
- Text file containing one or more guides (including a number of Ns equal to the length of
the PAM sequence) (e.g. TCACCCAGGCTGGAATACAGNNN, the last 3 Ns represents
the space occupied by the PAM in the real sequence)
- Name of output file (e.g. `emx1.hg19`)
- Tag to activate index search (`-index`)
- Number of allowed mismatches (e.g. `-mm 4`)
- Size of DNA bulges and/or RNA bulges (e.g. `-bDNA 1 -bRNA 1`)
- Output type (-r off-targets list only, -p profile only, -t everything) (e.g `-t`)
- Scores (-scores followed by the directory of the fasta genome, to perform the score after
the search with score calculation based on Doench 2016 and CFD, the two scoring
methods supports only NGG PAM and 23 long guides) (e.g `-scores hg19_ref/`)
- Number of threads to use for the analysis (Optional)

Output:
- Set of result files, including:
    - Targets file, containing all genomic targets for the guides set
    - Profile file, containing a matrix-like representation of guides behaviour (bp/mm, total on-/off- target, targets per mismatch threshold)
    - Extended profile file, containing the motif matrix for each guide and each mismatch threshold, useful to create visual analysis of the guides behaviour
    - Targets file with associated CFD score
    
Example call:
- Conda
    ```
    crispritz.py search genome_library/NGG_hg19_ref/ pam/pamNGG.txt guides/EMX1.txt emx1.hg19 -index -mm 4 -bDNA 1 -bRNA 1 -t -scores hg19_ref/
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search genome_library/NGG_hg19_ref/ pam/pamNGG.txt guides/EMX1.txt emx1.hg19 -index -mm 4 -bDNA 1 -bRNA 1 -t -scores hg19_ref/
    ```

Detailed input:  
`genome_library/NGG_hg19_ref/` is the directory containing the fasta files.  
`pam/pamNGG.txt/` is a text file containing the PAM sequence.  
`guides/EMX1.txt` is a text file containing the EXM1 guide  
`emx1.hg19` is the output file name  
`-index` tag to activate the index search  
`-bDNA 1` DNA bulges threshold  
`-bRNA 1` RNA bulges threshold  
`-mm 4` Mismatches threshold  
`-t` to select the output type  
`-scores hg19_ref/` to activate the calculation of score (Doench 2016 and CFD)  
`-th 4` is the number of threads to use (Optional)


**<a name="Annotation">3.4</a> CRISPRitz Annotation Tool:**  
This tool is created to perform genomic annotation on results obtained during the search phase.  
Input:
- Targets file, containing all genomic targets for the guides set (Phase [3.3.1](#Search_mm) / [3.3.2](#Search_mm_bul))
- Bed file containing the annotations
- Name of output file
- Samples ID file, containing the list of samples with their associated Population and Superpopulation (Optional)

Output:
- Set of files, including:
    - Targets file with annotation (identical file as the targets file in input) with an added column containing the annotations).
    - One summary file, counting all the annotations per mismatch number.
    
Example call:
- Conda
    ```
    crispritz.py annotate-results emx1.hg19.targets.txt annotations.bed emx1.hg19.annotated
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py annotate-results emx1.hg19.targets.txt annotations.bed emx1.hg19.annotated
    ```

Detailed input:  
`emx1.hg19.targets.txt` is the text file containing targets from previous search  
`annotations.bed` is the text file containing the genomic annotations  
`emx1.hg19.annotated` name of the output file  
`--change-ID samples_1000genomeproject.txt` file containing the samples and their associated Population and Superpopulation


**<a name="Generate-Report">3.5</a> CRISPRitz Generate-Report Tool**  
This tool is created to generate a visual representation of guide behaviour such as on-/off- target
activity in specific genomic regions, total number of on-/off- targets in reference and
variant genome and so on.  
Input:
- A guide present in the analyzed set (Phase [3.3.1](#Search_mm) / [3.3.2](#Search_mm_bul))
(e.g. `GAGTCCGAGCAGAAGAAGAANNN`)
- Number of mismatches to analyze (e.g. `-mm 4`)
- Annotation summary file, containing the counting of all the annotations per mismatch number (Phase [3.4](#Annotation))
- Extended profile file (Phase [3.3.1](#Search_mm) / [3.3.2](#Search_mm_bul))
- Tag to activate gecko dataset comparison (e.g. `-gecko`)
- Annotation reference summary file, containing the counting of all the annotations per mismatch number (See Post-Process Phase for more informations) 

Output:
- Pdf file containing the radar chart and motif logo for a guide, the radar chart shows how
much the guide is similar, in terms of number of targets found, to all guides in its dataset
(or the gecko dataset if selected).
- Barplot with a distribution of on-/off- targets in each annotation and a comparison
between variant and reference genome, in terms of total targets found.


Example call:
- Conda
    ```
    crispritz.py generate-report GAGTCCGAGCAGAAGAAGAANNN -mm 4 -annotation emx1.hg19.annotated.Annotation.summary.txt -extprofile emx1.hg19.extended_profile.xls -gecko
    ```
    
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py generate-report GAGTCCGAGCAGAAGAAGAANNN -mm 4 -annotation emx1.hg19.annotated.Annotation.summary.txt -extprofile emx1.hg19.extended_profile.xls -gecko
    ```

Detailed input:  
`GAGTCCGAGCAGAAGAAGAANNN` is a guide present in the result file, the one you want to analyze and print visualization files  
`-mm 4` Mismatches threshold  
`-annotation emx1.hg19.annotated.Annotation.summary.txt` is the file containing the counting of all the annotations per mismatch number
`-extprofile emx1.hg19.extended_profile.xls` is the xls file containing information detailed information about guides, used to construct the motif logo   
`-gecko` tag to activate the gecko dataset comparison, the results of your test guide, will be
compared with results from a previous computed analysis on gecko library.  


**<a name="Process-Data">3.6</a> CRISPRitz Process-Data Tool**  
This tools processes and compares the results obtained in the reference and enriched genomes, providing a more complete overview of the activity of
the input guides. The first step aggregates similar off-targets based on their position on the chromosome, then a list of samples ID is associated to each off-target haplotype. Finally, summaries about the guide activity at general, population and sample level are produced.

Input:
- Targets file of the reference genome, containing all genomic targets for the guides set (Phase [3.3.1](#Search_mm) / [3.3.2](#Search_mm_bul)) 
- Targets file of the enriched genome, containing all genomic targets for the guides set (Phase [3.3.1](#Search_mm) / [3.3.2](#Search_mm_bul)) 
- Text file containing the PAM sequence (including a number of Ns equal to the guide
length) and a space separated number indicating the length of the PAM sequence (e.g.
Cas9 PAM is NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by
20 Ns and NGG, followed by 3, representing the length of the PAM sequence.
- Text file containing one or more guides (including a number of Ns equal to the length of
the PAM sequence) (e.g. TCACCCAGGCTGGAATACAGNNN, the last 3 Ns represents
the space occupied by the PAM in the real sequence)
- Bed file containing the annotations
- Name of output file
- Directory containing `.json` files (Dictionaries), used for sample ID association. The directory is generated using the `--sample-create` option, that takes in input the directory containing the VCF files.
- Directory containing a genome in fasta format, need to be separated into single
chromosome files.
- Samples ID file, containing the list of samples with their associated Population and Superpopulation (Optional)

Output:
- Targets file, containing all post-processed genomic targets for each guide in the guides set (`targets.GUIDE.txt` files)
- Summary at Superpopulation, Population and Sample level counting all the annotations per mismatch number (`sample_annotation` files)
- Summary at guide level counting all the annotations per mismatch number (`sumref.Annotation.summary.txt` for reference genome, `Annotation.summary.txt` for variant genome)
- Count of targets, for each mismatch + bulge value, divided by Population (`PopulationDistribution.txt`), and represented by a series of barplots (`populations_distribution_GUIDE.png` files)
- Summaries of the guide activity based on Genomic Position (`summary_by_position.txt`), Samples (`summary_by_samples.txt`), Genome (`summary_by_guide`) ,and Guide (`general_target_count.txt`) 


Example Call:
- Conda
    ```
    crispritz.py process-data -reftarget emx1.hg19.ref.targets.txt -vartarget emx1.hg19.var.targets.txt pam/pamNGG.txt guides/EMX1.txt annotations.bed emx1.hg19.final -sample dictionaryDirectory/ -refgenome hg19_ref/
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py process-data -reftarget emx1.hg19.ref.targets.txt -vartarget emx1.hg19.var.targets.txt pam/pamNGG.txt guides/EMX1.txt annotations.bed emx1.hg19.final -sample dictionaryDirectory/ -refgenome hg19_ref/
    ```

Detailed Input:

`emx1.hg19.ref.targets.txt` is the text file containing targets from previous search on the reference genome
`emx1.hg19.var.targets.txt` is the text file containing targets from previous search on the enriched genome
`pam/pamNGG.txt` is a text file containing the PAM sequence.  
`guides/EMX1.txt` is a text file containing the EXM1 guide
`annotations.bed` is the text file containing the genomic annotations
`emx1.hg19.final` name of the output file
`dictionaryDirectory/` directory containing the `.json` files used for sample ID association. If no directory is available, change the option `-sample dictionaryDirectory/` to `--sample-create <vcfFileDirectory>` in order to create the `.json` files. In this case, the `vcfFileDirectory` is the directory containing the vcf files used for genome enrichment (Phase [3.1](#Add-Variant)).
`hg19_ref/` is the directory containing the fasta files


**Output examples:**
- Targets file, containing all genomic targets for the guides set
![example_targets](https://user-images.githubusercontent.com/32717860/53101471-19d81180-352a-11e9-9ee3-69de580c5e3f.PNG)
- Profile file, containing a matrix-like representation of guides behaviour (bp/mm, total on-/off- target, targets per mismatch threshold)
![profile](https://user-images.githubusercontent.com/40895152/63215013-866e4a80-c120-11e9-8855-c63c2a6e991a.png)
- Extended profile file, containing the motif matrix for each guide and each mismatch
threshold, useful to create visual analysis of the guides behaviour
![ex_profile](https://user-images.githubusercontent.com/40895152/63215020-a30a8280-c120-11e9-98ed-10a6ab47bc54.png)





- Pdf file containing the radar chart and motif logo for a guide, the radar chart shows how
much the guide is similar, in terms of number of targets found, to all guides in its dataset
(or the gecko dataset if selected).
![fig_medium_guide-1](https://user-images.githubusercontent.com/32717860/53101072-5b1bf180-3529-11e9-8c9a-cb5895f2c6c0.png)

- Barplot with a distribution of on-/off- targets in each annotation and a comparison
between variant and reference genome, in terms of total targets found.
<img width="645" alt="Capture" src="https://user-images.githubusercontent.com/32717860/62580929-e9055200-b8a7-11e9-8c6d-a83e229905b0.PNG">
