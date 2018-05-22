# CRISPRITZ

HOW TO:

In this package we include two algorithms to perform analysis and CRISPR off-target searching on genome.
The two algorithms are created to perform different types of analysis:
1- TST-Search, this algorithm is created to perform long-run analysis, so if you want to test a huge library of guides against your genome.
To use use this algorithm, you first need to create a genome_library.

Create the library:

python3 crispritz.py index-gen HG_19_indexed ../chroms_hg19/ ../input/PAM.txt 

with this command you create your library, 
in this case the genome name will be: HG_19_indexed
the reference genome will be: ../chroms_hg19/
and the PAM used to create the index will be: ../input/PAM.txt (a single line file containing your pam, for example 'NNNNNNNNNNNNNNNNNNNNNGG 3', 3 is the PAM length, so NGG)

After that, you can start using the genome library.

Perform Search:

python3 crispritz.py index-search HG_19_indexed NGG input/10guide4mm.txt -mm 4 -bDNA 1 -bRNA 1

with this command you can search your indexed genome,
the genome searche will be: HG_19_indexed
the pam will be: NGG
the guides will be: input/10guide4mm.txt (for example, 'GCTCAGTTACGGGTAAAACTNNN')
you can set the threshold for mm,in this case: -mm 4
the number of DNA bulges: -bDNA 1
and the number of RNA bulges: -bRNA 1

Those are the two commands required to create the genome library and to start the search.

2- Brute Force, this algorithms is created to perform very fasta search with no pre-processed files, you only need a guides list, a pam, and chromosomes files in .fa format.
This is the only command you need:

python3 crispritz.py search ../chroms_hg19/ ../input/PAM.txt ../guides/10guide4mm.txt 4 -mm 4 -r

the genome used in this search will be: ../chroms_hg19/
the pam will be: ../input/PAM.txt
the guides: ../guides/10guide4mm.txt
the number of used threads(in case of multi-thread machine) will be: 4
the missmatch threshold will be: -mm 4
and you can also choose if you want to print off-targets results: -r

Those are the simple commands required to start using those algorithms.

SOFTWARE REQUIREMENTS:
PYTHON 3
GCC 5.4
