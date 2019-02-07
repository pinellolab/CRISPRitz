# CRISPRITZ

CRISPRitz is a software package containing 5 different tools dedicate to perform analysis on CRISPR/Cas off-target analysis.
The aim of the software is to help the use performing tediuos and long tasks, such as, the search for off-targets on a genome, or the assessment of guides, in an easy and fast way.

With this aim in mind we create this package, containing the following tools:

add-variants: Created to permit the user to encode genomic variants using the IUPAC code, outputin a genome with added variant in IUPAC notation form.
index-genome: Created to permit the user to find the potential targets in a genome based on the PAM and saving them in a set of .bin files, avoiding the necessity to perform this recursive task on every search with the same PAM
search: The main function of the package, created to help the user perform off-target search on a genome(with variants) or genome index to help create safe and secure guides to use in a wet-lab CRISPR/Cas experiment.
annotate-results: Created to permit the user to annotate regions founded during search with functional annotations (promoter, chromatin accessibility, insulator, etc)
generate-report: Created to permit the user to generate a graphical report with annotated and overall mismatch and bulge profile for a given guide, to help the assessment of guides and to have an easy view over guides and their behave.

![Alt text](relative/path/to/img.jpg?raw=true "Title")

SOFTWARE REQUIREMENTS:
PYTHON 3
GCC 5.4
