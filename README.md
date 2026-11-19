# CRISPRitz 
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/crispritz/README.html)

CRISPRitz is the search engine of [CRISPRme+](https://github.com/pinellolab/crisprme-plus). It
holds three tools, which enrich a genome with variants, index it, and search it for the targets
and off-targets of CRISPR/Cas guides:

- **add-variants** writes a copy of a genome with the SNPs of a VCF in IUPAC notation (the
  enriched genome) and, on request, the dictionaries that list each SNP's carriers and the
  fake-indel genome, each indel allele a contig of its own (section [3.1](#Add-Variant)).
- **index-genome** builds the index of a genome, with or without variants, for a PAM: every
  candidate target, stored for searches with up to a given number of DNA and RNA bulges
  (section [3.2](#Index-Genome)).
- **search** searches a genome for the targets of one or more guides: a FASTA genome with
  mismatches only, or an index with mismatches and bulges. It writes the targets and a profile
  of each guide (section [3.3](#Search)).

`crispritz.py index-info` prints the format and build parameters of an index (section
[3.2.1](#Index-Format)), and `crispritz.py version` the version and the search hooks it has
(section [3.3.4](#Search_regions)).

CRISPRitz scores, annotates and reports nothing. 2.9.0 removed `scores`, `annotate-results`,
`generate-report` and the `-scores` option of `search` (`CHANGELOG.md`, 2.9.0), and
`process-data` had gone before 2.6.0. CRISPRme+ does all of that: CFD and CRISPR-Bulge scores, functional
annotation, the samples and populations that carry each target, and its report (section
[3.4](#CRISPRme)).

# CRISPRitz Installation and Usage
The two fastest way to use CRISPRitz is through the installation of Docker or Conda.
Here we summarize the steps to install CRISPRitz with Docker and Conda.

> ### Platform support
> - **Linux (x86_64 and arm64):** fully supported natively — install via **Bioconda / Conda** (`conda install crispritz`).
> - **macOS and Windows:** use the **Docker** image (see below). The image is multi-architecture and runs **natively on Apple Silicon Macs** (`linux/arm64`, no emulation).
> - **Native macOS note:** `add-variants` (the C++ enricher) builds and runs natively on Apple Silicon, but the index engine (`index-genome` / `search`, i.e. `buildTST`) relies on GCC's libstdc++ parallel extensions (`__gnu_parallel`), which Apple `clang` does not provide. For indexing/search on macOS, **use the Docker image** (it ships a Linux toolchain). This is also why the Bioconda package is **Linux-only** (`linux-64` + `linux-aarch64`).

> ### CRISPRitz and CRISPRme+ versions
> CRISPRitz and CRISPRme+ share their major.minor version; their patch numbers are
> independent. CRISPRme+ X.Y.z requires CRISPRitz X.Y.\*: its Bioconda package pins
> `crispritz >=X.Y.0,<X.(Y+1)`, and its Docker image builds the CRISPRitz tag it was tested
> with. A breaking change in either package (an index format, the enricher's output) bumps
> the minor of both. CRISPRitz 2.10.0 goes with CRISPRme+ 2.10: it writes index format 2,
> which CRISPRme+ 2.8.1 refuses and 2.8.0 aborts on (section [3.2.1](#Index-Format)). Its
> Bioconda package therefore declares `run_constrained: crisprme >=2.10.0`, so that conda does
> not install it beside an earlier CRISPRme+ (`docs/RELEASING.md`, Version policy).

## Installation (Phase 1)
**Conda installation (Linux only):**
- Open a terminal window
- Paste this command into the terminal (Linux):
    ```
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --output Miniconda3-latest-Linux-x86_64.sh
    ```
- If the file is correctly downloaded you now need to execute it to complete the installation, so paste this command into the terminal:
    - Linux
    ```
    bash Miniconda3-latest-Linux-x86_64.sh
    ```
- Press ENTER when requested and yes when an answer is requested, in this way you allow conda to set all the directories in your HOME path for an easy use
- After the complete installation you will receive this message “Thank you for installing Miniconda3!” to certify the correct installation.
- Now you need to close the terminal window you are using and open a new one, to allow the system to start conda.
- In the new terminal window you should see something like this:
    ```
    (base) user@nameofPC:~$
    ```
    If you read the "(base)" like this, conda is loaded correctly and you can start using it.
- Now you need to set the channels to allow conda to access different repositories and set the default version of python to version 3.11, so paste these commands into the terminal you just opened:
    ```
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    conda install python=3.11
    ```
- Now, you can install CRISPRitz by typing the command:
    ```
    conda install crispritz
    ```
- To test your installation, type the command:
    ```
    crispritz.py
    ```
- After the execution of the command you should see the version and the list of CRISPRitz
  tools:
    ```
    version 2.10.0
    help:
    ...
    crispritz.py add-variants FUNCTION TO ADD VARIANTS DATA TO A FASTA GENOME
    crispritz.py index-genome FUNCTION TO CREATE GENOME INDEX TO PERFORM FAST SEARCHES WITH BULGES
    crispritz.py search FUNCTION TO PERFORM SEARCHES ON A GENOME (INDEXED OR PLAIN FASTA)
    crispritz.py index-info FUNCTION TO PRINT THE FORMAT AND BUILD PARAMETERS OF A GENOME INDEX
    ...
    ```
Now the software is installed and ready to be used.

**Docker installation (Linux, Mac OSX and Windows):  
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
Check that the installed CRISPRitz runs and is the version you expect:
```
crispritz.py version
```
It prints the version and the search hooks this version has:
```
CRISPRitz v 2.10.0
search hooks: skip-bed keep-bed leaf-dump
```
With Docker: `docker run --rm pinellolab/crispritz crispritz.py version`. Then run the examples
of Phase [3](#phase3) on the test files: each tool writes the outputs its section lists.

The scripts `test_scripts/auto_test_crispritz_conda.sh` and
`test_scripts/auto_test_crispritz_docker.sh`, which earlier versions of this README pointed to,
still run the tools 2.9.0 removed and check no exit status, so their closing message
("EVERY TEST PASSED") says nothing about the installation.

## Usage (Phase <a name="phase3">3</a>):
Here is a brief guide to help use CRISPRitz. Download the test files first (one time step):
- The script downloads chr22 of hg38 from UCSC (into `hg38_ref/`), the phased chr22 VCF of the
1000 Genomes Project on GRCh38 (biallelic SNVs, into `hg38_1000genomeproject_vcf/`), the EMX1
guide (`EMX1.sgRNA.txt`), an NGG PAM file for 20-nt guides (`20bp-NGG-SpCas9.txt`), and an hg38
annotation BED that no CRISPRitz tool reads since 2.9.0, all into `test_crispritz/`.
- Download the script with this command:
    ```
    curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/download_test_files.sh --output download_test_files.sh
    ```
- Write this command to execute the script:
    ```
    bash download_test_files.sh
    ```
- The script downloads one chromosome and one VCF file, to save time. All the examples can be
run on an entire genome: add the other chromosomes, one FASTA file each, into `hg38_ref/`, and
their VCFs into `hg38_1000genomeproject_vcf/`.
- Write this command to enter the test directory:
    ```
    cd test_crispritz/
    ```
- Now you are ready to execute the following example functions.

**<a name="Add-Variant">3.1</a> CRISPRitz Add-Variant Tool**  
This tool is created to insert variants in a fasta genome.  
Input:
- Directory containing VCF files, need to be separated into single chromosome files
  (multi-sample files will be collapsed into one fake individual), each gzipped and named
  with its chromosome as a part between dots: the last part that holds `chr`
  (`ALL.chr22.<...>.vcf.gz` is chr22's).
- Directory containing a genome in fasta format, one `<chromosome>.fa` file per chromosome
(`chr22.fa`), the chromosome as the VCF's name gives it.
- `true`, to write the dictionaries and the fake-indel genome as well (Optional; CRISPRme+
always passes it).

Output, in `variants_genome/` of the working directory:
- `SNPs_genome/<genome directory>_enriched/`: a copy of the genome, one
`<chromosome>.enriched.fa` file per chromosome, with the SNPs in IUPAC notation (a chromosome
without a VCF is copied as it is).
- With `true`, for each chromosome with a VCF: the SNP dictionary
`SNPs_genome/my_dict_<chromosome>.json`, which lists each SNP's carriers, the indel log
`SNPs_genome/log<chromosome>.txt`, and the fake-indel genome
`fake_<VCF directory>_<chromosome>/fake<chromosome>.fa`, each indel allele a contig of its own.

Each indel allele (an indel or a multi-base record such as an MNV, of a record whose FILTER is
`PASS` or `.` and that some sample carries) is written as a contig of its own, its fake-indel
window: since 2.10.0, 30 reference bases, the record's REF (from POS on) with the ALT in its
place, then 30 reference bases (fewer where the chromosome ends), one `N` between two
contigs; up to 2.9.1, 25 bases before and 27 after. A target of L bases (guide and PAM) with
b DNA bulges holds a record's change on its window only where L + b <= 31 (L + b <= 27, 26
at an MNV, up to 2.9.1). The indel log (`log<chr>.txt`) records each window: its `CHR` field,
`<chr>_<start>-<end>_<id>`, holds the window's 0-based start and end (the end exclusive and
not cut at the chromosome's end), so a reader takes its left flank as `POS - 1 - start` and
reads the windows of either release; `crispritz.py`'s `INDEL_LEFT` and `INDEL_RIGHT` state
the window an enrichment not made yet will have. An indel less than 31 bases from the
chromosome start (POS 30 or less) is skipped, with a `WARNING` on stdout.

Example call:
- Conda
    ```
    crispritz.py add-variants hg38_1000genomeproject_vcf/ hg38_ref/ true
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py add-variants hg38_1000genomeproject_vcf/ hg38_ref/ true
    ```

Detailed input:  
`hg38_1000genomeproject_vcf/` is the directory containing the vcf files.  
`hg38_ref/` is the directory containing the fasta files.  
`true` writes the SNP dictionaries, the indel logs and the fake-indel genome too.

Multi-allelic sites: a site may be one VCF record listing several ALTs or one record per ALT
(as `bcftools norm -m -any` writes it). Since 2.10.0, either way, the code of a position is
the IUPAC letter of its reference base (the reference is read upper-cased, as before, so a
soft-masked base is a base, and a chromosome with a VCF is written in upper case) and of
every ALT that is one base (`A`, `C`, `G` or `T`, in either case: see below) of its records
whose FILTER is `PASS` or `.` and whose REF is one base, in any record order, and its SNP
dictionary entry lists every ALT of one character of those records (a `*` or an `N` too,
which add nothing to the code); up to 2.9.1 a split site kept only its last record's ALT. A
reference IUPAC code adds no base of its own, and a single base leaves the reference letter
in place; a reference `N` stays `N`, whatever the VCF lists there, as does a byte that is no
base or IUPAC code (such as the `X` RepeatMasker `-x` masks with), which `index-genome` and
`search` read as `N`; and a position whose ALTs, with its reference base when that is A, C,
G or T, cover all four bases is coded `Q`, a fifth symbol `index-genome` and `search` read
as any base: its IUPAC letter, `N`, is a gap to them, which matches nothing (an index file
holding `Q` is flagged in its header, so that no engine that cannot read it searches it; a
reference FASTA must not hold `Q`, in either case, and `add-variants` refuses one). Apply an
allele-frequency filter after splitting, so it acts per ALT. Each record's `INFO/AF` is read
by its key (up to 2.9.1 its place came from the file's first record). Each dictionary entry
lists its own record's carriers, their genotypes in split form: `1` is the entry's ALT and
`0` the reference or another ALT, a missing allele and the separators as they were. So a
split record's genotypes are kept verbatim, and a record listing several ALTs is written ALT
by ALT: a `1|2` carrier of `C>A,T` reads `1|0` under A and `0|1` under T, where 2.9.1 wrote
`1|2` under both. One allele on two records is one entry with the carriers of both. Split
records phased one by one can put two ALTs of a position on one haplotype slot (`0|1` under
each), and whatever builds haplotypes from the dictionary has to resolve that.
Details: [`docs/MULTIALLELIC_ENRICHMENT.md`](docs/MULTIALLELIC_ENRICHMENT.md).

REF and ALT in lower case: VCF 4.x reads their bases whatever their case, and since 2.10.0
so does `add-variants`, which reads a record's REF and every ALT in upper case. An ALT holding
`<`, a symbolic allele, is kept as written: VCF 4.3 and 4.5 let a tool change the case of an
allele "except for IDs, which are case sensitive". `*` and `.` hold no letter. A breakend ALT
(holding `[` or `]`) is upper-cased whole, the name of its mate's contig with it; that changes
only the text of its log row and window, which `index-genome` reads upper-cased. So a SNP ALT
written in lower case (`C>g`) is coded. An SNV record that differs from another only by case
is the same allele (one dictionary entry with the carriers of both, as for two records written
alike), and two indel records that do keep a window and a log row each, now under one name, as
two written alike do. An indel's fake-indel window, its indel log row and every dictionary
entry name the record in upper case; rsID, AF and genotypes are read as before. Up to 2.9.1
the alleles were read as written: a lower-case SNP ALT was not coded while the dictionary
listed its carriers, and a window held its ALT, and a log row its REF and ALT, in the VCF's
case. Lower case comes in through normalisation: `bcftools norm -f` writes a soft-masked
reference's lower case into the indels it left-aligns. HPRC normalised as CRISPRme builds it
holds 108,178 such records over its 24 chromosomes, 1,546 of them on chr22, every one an indel
or MNV; HF's raw vcfwave VCFs hold none.


**<a name="Index-Genome">3.2</a> CRISPRitz Index-Genome Tool**  
This tool is created to generate an index genome (similar to the bwa-index step). Indexing
takes time and memory, but saves a lot of execution time while searching with lot of guides
and with the support of bulges (DNA and RNA): measured on hg38 (docs/INDEX_FORMAT.md section
10), chr22's NGG index with `-bMax 2` took 6 to 16 s and 1.5 GB, and chr1's PAM-less index of
an 8-nt flank 18 minutes and 103 GB. If do not need to search with bulges, skip this passage.  
Input:
- Name of the genome to create (e.g. `hg38_ref`).
- Directory containing a genome in fasta format, one `<chromosome>.fa` file per chromosome
(`index-genome` indexes no other file).
- Text file containing the PAM (including a number of Ns equal to the guide length) and a
space separated number indicating the length of the PAM sequence (e.g. Cas9 PAM is
NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by 20 Ns and
NGG, followed by number 3, representing the length of the PAM sequence.
- Number of bulges to include in the database to perform the following search (i.e. the max
number bulges allowed for DNA and RNA when searching on the database)
- Number of threads to use for the analysis (Optional; default 4)

Output:
- Directory `genome_library/<PAM>_<bulges>_<name>/` in the working directory
(`genome_library/NGG_2_hg38_ref/` here): the index genome in .bin format, one or more files per
chromosome, containing all the candidate targets for a selected PAM, adding also characters to
perform bulge search.

Example call:
- Conda
    ```
    crispritz.py index-genome hg38_ref hg38_ref/ 20bp-NGG-SpCas9.txt -bMax 2 -th 4
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py index-genome hg38_ref hg38_ref/ 20bp-NGG-SpCas9.txt -bMax 2 -th 4
    ```
Detailed input:  
`hg38_ref` is the name of the index genome, the last part of its directory's name.  
`hg38_ref/` is the directory containing the fasta files.  
`20bp-NGG-SpCas9.txt` is a text file containing the PAM sequence.  
`-bMax 2` is the max number of bulges to allow for following searches when creating the
database: with `-bMax 2`, every search of the index can allow up to 2 RNA bulges and 2 DNA
bulges. It is at most 31. `crispritz.py search` lowers a larger `-bDNA` or `-bRNA` to the
number the index directory's name states (the `2` of `NGG_2_hg38_ref`), with a `WARNING`.
Since 2.10.0, `searchTST` refuses a search of more bulges than an index's header states
(format 2) or than its own `max_bulges` argument, naming both numbers, so
`crispritz.py search` stops too on a format-2 index whose directory was renamed to state more
bulges than its files were built with. Before, such a call crashed or exited 0, often without
alignments that an index built for those bulges returns (`CHANGELOG.md`, 2.10.0, Fixed).  
`-th 4` is the number of threads to use (Optional)

**<a name="Index-Format">3.2.1</a> Index format 2 (CRISPRitz 2.10.0 and later).**  
An index built by CRISPRitz 2.10.0 or later is in *index format 2*, specified in
[docs/INDEX_FORMAT.md](docs/INDEX_FORMAT.md). It holds every target whose guide region (of
the length the index was built for) and PAM fit in the sequence, with no `N` in the guide
region, including targets next to a contig end or a run of `N`: the extra bases a DNA bulge,
or a longer search guide, needs beyond the guide (`-bMax` of them) no longer have to exist,
they are simply not available to an alignment that would need them. Up to 2.9.1 a target
whose guide ended within `-bMax` nt of a contig end or an `N` run, on its side away from the
PAM, was dropped, even a perfect match; a PAM on a contig's last bases was never indexed,
nor a 5′ PAM on the `+` strand at a contig's first base, and a reverse-strand PAM at a
contig's first base was reported on the `+` strand. Each `.bin` now starts with a header
naming its format and build parameters, which `crispritz.py index-info` prints without
loading the index (it exits 1, naming the file and the reason, on a `.bin` whose first and
last bytes or size show that any search would refuse it):

```
crispritz.py index-info genome_library/NGG_2_hg38_ref/
```

Compatibility: CRISPRitz 2.10.0 searches indexes built by earlier versions (format 1) exactly
as before, except at a contig's first base: a reverse-strand PAM there is reported on the
correct strand, and a spurious match the old builder could store there is dropped. Such
indexes keep missing the boundary targets until rebuilt. An index built by 2.10.0 cannot be
searched by CRISPRitz 2.9.1 or earlier: 2.9.1 refuses it with an error, earlier versions
abort with one. An index of an enriched genome may hold `Q`, the four-allele symbol; the
files that do are flagged in their header (`index-info` reports `four_allele`).


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
- Name of the output file (e.g. `emx1.hg38`)
- Number of allowed mismatches (e.g. `-mm 4`)
- Output type (-r off-targets list only, -p profile only, -t everything) (e.g `-t`)
- Number of threads to use for the analysis (Optional; by default every available one)

Output:
- Set of result files, including:
    - Targets file, containing all genomic targets for the guide set
    - Profile file, containing a matrix-like representation of guides behaviour (bp/mm, total on-/off- target, targets per mismatch threshold)
    - Extended profile file, containing the motif matrix for each guide and each mismatch threshold, useful to create visual analysis of the guides behaviour

CRISPRitz scores no target: CRISPRme+ computes the CFD and CRISPR-Bulge scores (section
[3.4](#CRISPRme)).

Example call:
- Conda
    ```
    crispritz.py search hg38_ref/ 20bp-NGG-SpCas9.txt EMX1.sgRNA.txt emx1.hg38 -mm 4 -t -th 4
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search hg38_ref/ 20bp-NGG-SpCas9.txt EMX1.sgRNA.txt emx1.hg38 -mm 4 -t -th 4
    ```

Detailed input:  
`hg38_ref/` is the directory containing the fasta files.  
`20bp-NGG-SpCas9.txt` is a text file containing the PAM sequence.  
`EMX1.sgRNA.txt` is a text file containing the EMX1 guide  
`emx1.hg38` is the output file name  
`-mm 4` to select the mismatch threshold  
`-t` to select the output type  
`-th 4` is the number of threads to use (Optional)

**<a name="Search_mm_bul">3.3.2</a> Mismatches + Bulges search:**  
Input:
- Directory containing an index genome in .bin format, as `index-genome` writes it
(`genome_library/<PAM>_<bulges>_<name>/`, Phase [3.2](#Index-Genome)): `crispritz.py search`
reads the index's number of bulges from that name.
- Text file containing the PAM sequence (including a number of Ns equal to the guide
length) and a space separated number indicating the length of the PAM sequence (e.g.
Cas9 PAM is NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by
20 Ns and NGG, followed by 3, representing the length of the PAM sequence.
- Text file containing one or more guides (including a number of Ns equal to the length of
the PAM sequence) (e.g. TCACCCAGGCTGGAATACAGNNN, the last 3 Ns represents
the space occupied by the PAM in the real sequence)
- Name of output file (e.g. `emx1.hg38.bulges`)
- Number of allowed mismatches (e.g. `-mm 4`)
- Size of DNA bulges and/or RNA bulges (e.g. `-bDNA 1 -bRNA 1`). Either option selects the
index search, `-bDNA 0 -bRNA 0` one without bulges; without them `crispritz.py search` runs
the mismatches only search of 3.3.1, which reads no index and stops with an `ERROR` on an index
directory. (`-index`, which earlier versions of this README passed, is ignored.)
- Output type (-r off-targets list only, -p profile only, -t everything) (e.g `-t`)
- Number of threads to use for the analysis (Optional; default 1)

Output:
- Set of result files, including:
    - Targets file, containing all genomic targets for the guides set
    - Profile file, containing a matrix-like representation of guides behaviour (bp/mm, total on-/off- target, targets per mismatch threshold)
    - Extended profile file, containing the motif matrix for each guide and each mismatch threshold, useful to create visual analysis of the guides behaviour
    
Example call:
- Conda
    ```
    crispritz.py search genome_library/NGG_2_hg38_ref/ 20bp-NGG-SpCas9.txt EMX1.sgRNA.txt emx1.hg38.bulges -mm 4 -bDNA 1 -bRNA 1 -t -th 4
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search genome_library/NGG_2_hg38_ref/ 20bp-NGG-SpCas9.txt EMX1.sgRNA.txt emx1.hg38.bulges -mm 4 -bDNA 1 -bRNA 1 -t -th 4
    ```

Detailed input:  
`genome_library/NGG_2_hg38_ref/` is the directory containing the index genome.  
`20bp-NGG-SpCas9.txt` is a text file containing the PAM sequence.  
`EMX1.sgRNA.txt` is a text file containing the EMX1 guide  
`emx1.hg38.bulges` is the output file name  
`-bDNA 1` DNA bulges threshold  
`-bRNA 1` RNA bulges threshold  
`-mm 4` Mismatches threshold  
`-t` to select the output type  
`-th 4` is the number of threads to use (Optional)


**<a name="Search_prune_pam">3.3.3</a> One generous-flank PAM-less index, any PAM via `--prune-pam`:**  

Building (and storing/distributing) a separate index for every PAM (NGG, NRG, NGA,
NNGRRT, TTTV, …) is wasteful: the costly part of an index is the protospacer tree,
not the PAM. CRISPRitz can instead build **one PAM-less index with a generous all-N
flank of width K**, then serve **any 3′ PAM of length L ≤ K** on the fly by N-padding the
search PAM up to K and filtering candidates at leaf-emit time. The targets are identical
(the same rows, coordinates and sequences) to those of a dedicated per-PAM index built for
the same guide length and `-bMax` (at contig ends too if the index is in format 2: see the
properties in §[3.3.3.2](#Search_prune_pam_search)); a 5′ PAM needs an index of exactly its
width — see §[3.3.3.2](#Search_prune_pam_search).

**<a name="Search_prune_pam_build">3.3.3.1</a> Build the generous-flank index.**  
The stored flank width **K is exactly the PAM length declared in the PAM file**
(`index-genome`/`search` read the trailing number). To serve any PAM up to K, build
with an **all-N PAM file** whose total length is `spacer_length + K` and whose declared
number is `K`. Example for a 20-nt spacer and a generous **K = 8**
(recommended default — covers NGG, NGA, NNGRRT, NNNRRT and almost every characterized
3′ PAM, with headroom):

```
# pam_generous8.txt  (20 N spacer + 8 N flank, 3' PAM)
NNNNNNNNNNNNNNNNNNNNNNNNNNNN 8
```

```bash
crispritz.py index-genome hg38_generous8 hg38_dir/ pam_generous8.txt -bMax 2
# -> genome_library/NNNNNNNN_2_hg38_generous8/
```

The index is built for the guide length of its PAM file (20 here). A search PAM file may
declare another guide length (e.g. 21 N's before `NNGRRT` for SaCas9): the search runs,
and says so on stdout, but its rows are those of an index built for 20-nt guides, which
with bulges can differ from those of an index built for the search's guide length (and a
shorter guide misses targets at contig ends and `N` runs) — see
[docs/INDEX_FORMAT.md](docs/INDEX_FORMAT.md) section 5.1. Build the generous index with
the guide length you search with.

A 5′ (Cas12a-style) PAM is served only by a PAM-less index of exactly its width, declared
with a negative number (e.g. `NNNNNNNNNNNNNNNNNNNNNNNN -4` for TTTV): a shorter 5′ PAM off a
wider flank is refused, because the reported sequence would carry the pad N's in its PAM.

**<a name="Search_prune_pam_search">3.3.3.2</a> Search any PAM ≤ K with `--prune-pam`.**  
Pass the real, short PAM file (its declared length is **L**, the true PAM length) and add
`--prune-pam`. The engine N-pads a 3′ search PAM up to the index flank width K at its
**distal** end (a real PAM is always adjacent to the protospacer, so its letters occupy the
proximal positions and the pad-N's cover the outer stored positions the short PAM does not
constrain), then skips any candidate whose concrete flank is incompatible with the search
PAM **before it is written or scored**:

```bash
# NGG off the K=8 generous index (no dedicated NGG index needed)
crispritz.py search genome_library/NNNNNNNN_2_hg38_generous8/ pamNGG.txt guides.txt out_ngg \
    -mm 4 -bDNA 1 -bRNA 1 -th 4 -t --prune-pam

# NGA, NNGRRT … all served off the SAME index, just swap the PAM file (each with the
# index's 20 N's, e.g. NNNNNNNNNNNNNNNNNNNNNNGRRT 6). Keep -bDNA/-bRNA:
# without them crispritz.py dispatches to the brute-force search, which does not read an index.
crispritz.py search genome_library/NNNNNNNN_2_hg38_generous8/ pamNGA.txt    guides.txt out_nga    -mm 4 -bDNA 1 -bRNA 1 -th 4 -t --prune-pam
crispritz.py search genome_library/NNNNNNNN_2_hg38_generous8/ pamNNGRRT.txt guides.txt out_sauri  -mm 4 -bDNA 1 -bRNA 1 -th 4 -t --prune-pam
```

Properties (enforced by `test/run_pam_filter_test.sh` and `test/run_boundary_test.sh`):
- **Equivalent to a dedicated index of the same guide length and `-bMax`.** The targets
  returned for e.g. NGG off the K = 8 index are the *same set of rows* as a dedicated
  3-flank NGG index's built for the same guide length and `-bMax`, coordinates and
  sequences included: the DNA/crRNA columns carry the L-base PAM a dedicated index reports.
  This holds at contig ends too for an index in format 2 ([3.2.1](#Index-Format)), whose
  PAM-less windows may run off the contig on the flank side, and for a search guide of
  any length. It does not hold against a dedicated index built for another guide length:
  21-nt NNGRRT guides off a 20-nt K = 8 index can give, with bulges, other rows than a
  21-nt NNGRRT index (docs/INDEX_FORMAT.md section 5.1). An index built by CRISPRitz
  2.9.1 or earlier (format 1) required the whole K-base flank to fit, so off it a PAM with
  1 to K - L bases to a contig's end (+ strand) or 0 to K - L - 1 bases from its start
  (- strand) is missing; K = L is unaffected.
  The profile outputs (`-p`/`-t`) are not comparable between the two: they count the
  engine's duplicate rows, whose number depends on the index tree's shape, and they label
  the guide with the search PAM padded to K.
- **A shorter 5′ PAM off a wider flank is refused**; an equal-width 5′ index is supported.
- **`L > K` is a hard error** (the index never stored those extra flank bases — rebuild
  with a wider K).
- **Default (no `--prune-pam`) is unchanged / byte-identical** to upstream.

**<a name="Search_prune_pam_size">3.3.3.3</a> Does a wider flank blow up the index? No — it is linear.**  
The flank is stored 4-bit-packed (2 bases/byte) per leaf, so a wider flank adds a *constant,
tiny* per-leaf cost, not a multiplicative one. Measured on a 5 Mb hg38 chr22 slice
(reference-only, mismatch index, ~10 M leaves):

| flank K | on-disk `.bin` | bytes / leaf |
|--------:|---------------:|-------------:|
| 3       | 216.0 MB       | 21.60        |
| 8       | 236.0 MB       | 23.60        |

Going 3 → 8 adds exactly **+2 bytes/leaf** (+9.3%) — the 5 extra bases pack into 2 extra
bytes. Extrapolated to a full hg38 reference index (~6.2 × 10⁹ leaves): ≈134 GB (K = 3)
→ ≈146 GB (K = 8), i.e. **+≈12 GB (+9.3%)**. A variant/enriched index has more leaves but
the same +2 bytes/leaf, so the percentage is unchanged. **A generous K = 8 flank is a
~10% size premium that eliminates the need to store one index per PAM — not an explosion.**

**<a name="Search_regions">3.3.4</a> Skipping or keeping regions: `--skip-bed`, `--keep-bed`, `--leaf-dump` (2.10.0).**  
An indexed search (with `-bDNA`/`-bRNA`, 0 allowed) can be given a BED of intervals
(tab-separated contig, 0-based start, end; contigs of the index only, which
`crispritz.py index-info` lists as `contigs`):

```bash
# every site whose core meets an interval is skipped; the rest is searched as usual
crispritz.py search genome_library/NGG_2_hg38_ref/ pamNGG.txt guides.txt out_skip \
    -mm 4 -bDNA 1 -bRNA 1 -th 4 -r --skip-bed dense.bed --leaf-dump out_skip.dump.tsv
# the inverse: only the sites whose core meets an interval are searched
crispritz.py search genome_library/NGG_2_hg38_ref/ pamNGG.txt guides.txt out_keep \
    -mm 4 -bDNA 1 -bRNA 1 -th 4 -r --keep-bed dense.bed
```

A site's **core** is its bulge-free alignment's search PAM and the search-guide bases next
to it (`Gs + L` bases: on the K=8 index searched for NGG, 23, without the 5 distal flank
bases; for a 21-nt guide on a 20-nt index, 24, with the pad base next to the guide). Read
from an output line, it is `[Position + d - r, + n)` for a PAM-right line (a 3′ PAM on
`+`, a 5′ PAM on `-`) and `[Position, + n)` for a PAM-left one, `d` and `r` the gaps in
its crRNA and DNA, `n` the crRNA's length less `d`. With `--skip-bed` the search returns
exactly the lines it returns without the BED whose core meets no interval; with
`--keep-bed`, exactly the others: the two partition it. The skipped subtrees of the index
are never entered, so a skip of the variant-dense windows costs a fraction of the search
it replaces (on chr22's IGL locus with HGDP + 1000 Genomes, one guide at 4 mismatches and
2+2 bulges on the INDELS index of 2.10.0's 30/30 fake-indel windows: 3.4 s and 0.06 GB with
the BED of its dense windows; without it, over 51 GB of lines and 166 GB of memory, where
2.9.1's 25/27 windows took 24 minutes and 157 GB). `--leaf-dump` (with `--skip-bed`) writes
one line per skipped site: contig, window start, the letters the index stores for its window
on the + strand
(`-` where it stores a placeholder: a bulge-pad base from the first one off the contig or
on an `N` on, a flank base off the contig) and `L` or `R`, the window's side. A search
given a BED prints a line of four tab-separated fields, `SKIPPED LEAVES`, its count,
`KEPT LEAVES` and its count (and, with a dump, a line `LEAF DUMP` and its count, tab-separated
too), `crispritz.py version` lists `search hooks: skip-bed keep-bed leaf-dump`, and a malformed
BED, a contig the index has no tree for, a skip with a keep BED, or a dump without a skip BED
is refused before the search starts. Without these flags the search is unchanged, byte for
byte. The haplotype scan of dense regions that CRISPRme 2.10.1 is to bring uses them (CRISPRme
2.10.0 does not); the specification and its measurements are in
[docs/INDEX_FORMAT.md](docs/INDEX_FORMAT.md) section 13.


**<a name="CRISPRme">3.4</a> Scores, annotation, carriers and reports: CRISPRme+**  
CRISPRitz enumerates targets; it scores, annotates and reports nothing. The tools that did
(`scores`, `annotate-results`, `generate-report` and the `-scores` option of `search`) were
removed in 2.9.0, and `process-data`, which compared the searches of a reference and an
enriched genome and listed the samples of each target, before 2.6.0. Their work is done by
[CRISPRme+](https://github.com/pinellolab/crisprme-plus), which runs CRISPRitz and reads what
it writes (the enriched genome, the dictionaries, the indel logs, the indexes and the targets):
the CFD and CRISPR-Bulge scores, functional annotation, the samples, haplotypes and populations
that carry each target, and a report of each guide.


**Output examples:**
- Targets file, containing all genomic targets for the guides set: one line per alignment,
with its bulge type (`X` for none, `DNA` or `RNA`), the guide (crRNA, `-` at a DNA bulge) and
the target (DNA, its mismatches in lower case), its position and strand, and its mismatches,
bulges and both together (from the golden output of `test/run_golden_test.sh`):
    ```
    #Bulge type	crRNA	DNA	Chromosome	Position	Cluster Position	Direction	Mismatches	Bulge Size	Total
    X	GACTGCATTGCATGACTGCANNN	GACTGaATTGCAaGACTGCAGGG	chr1	685	685	+	2	0	2
    DNA	G-ACTGCATTGCATGACTGCANNN	aGACTGaATTGCAaGACTGCAGGG	chr1	684	685	+	3	1	4
    ```
- Profile file, containing a matrix-like representation of guides behaviour (bp/mm, total on-/off- target, targets per mismatch threshold)
![profile](https://user-images.githubusercontent.com/40895152/63215013-866e4a80-c120-11e9-8855-c63c2a6e991a.png)
- Extended profile file, containing the motif matrix for each guide and each mismatch
threshold, useful to create visual analysis of the guides behaviour
![ex_profile](https://user-images.githubusercontent.com/40895152/63215020-a30a8280-c120-11e9-98ed-10a6ab47bc54.png)
