#!/bin/bash

make -f Makefile_conda
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/opt/crispritz
chmod -R 700 .
cp crispritz.py $PREFIX/bin
cp buildTST $PREFIX/opt/crispritz
cp searchTST $PREFIX/opt/crispritz
cp searchBruteForce $PREFIX/opt/crispritz
cp -R sourceCode/Python_Scripts $PREFIX/opt/crispritz
# The enricher is the 4th binary Makefile_conda builds, and until 2.9.1 it was the
# only one never installed -- so a packaged install silently fell back to
# enricher.py via crispritz.py's _enricher_command(), while a source install
# (CRISPRme's install_from_source.sh) used the compiled one. Install it into BOTH
# places _enricher_command() checks, in its order: PATH first, then beside the
# reference implementation. zlib is already in the recipe's build/host/run deps and
# enricher.cpp carries no OpenMP pragmas, so this ships on every platform the other
# three already do, by the same compiler invocation.
#
# Only safe as of 2.9.1: before it, the port wrote a plain base where enricher.py
# writes an IUPAC code at any position carrying both a SNP and an indel, dropping
# the SNP from the enriched genome. verify_enricher.sh case 8 covers that shape and
# runs in CI as the enricher-ports job.
cp enricher $PREFIX/bin
chmod +x $PREFIX/bin/enricher
cp enricher $PREFIX/opt/crispritz/Python_Scripts/Enrichment/enricher
