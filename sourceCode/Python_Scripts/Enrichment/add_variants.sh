#!/usr/bin/env bash
# $1 is vcf file (Eg. variants_hg38/ALL.chr1.vcf.gz)
# $2 is chrname with file extension(Eg. chr1.fa)
# $3 is chrname (Eg. chr1)
#
# Invokes the compiled C++ `enricher`, a port of enricher.py -- the reference
# implementation kept alongside this script. The C++ binary is built by
# `make -f Makefile_conda crispritz`; falls back to the copy sitting next to this
# script if it is not on PATH.
#
# "Byte-identical" is a VERIFIED claim, not an assumption: verify_enricher.sh
# diffs every output of both ports over a synthetic VCF covering 8 variant
# shapes, and runs in CI as the `enricher-ports` job. It WAS an assumption until
# 2.9.1, and it was wrong -- the C++ port dropped a SNP wherever an indel shared
# its position, which the fixture could not see because its SNP+indel case put
# both alleles in ONE multiallelic record while real VCFs are split. See the
# 2.9.1 CHANGELOG. Treat any change to either port as requiring that script to
# pass.
#
# NOTE: as of 2.9.1 conda/build.sh installs this binary into both places
# crispritz.py's `_enricher_command` checks (PATH, then beside enricher.py), so
# the compiled port is what a packaged install runs; before that it compiled the
# binary but never installed it, and packaged installs ran enricher.py. Which
# implementation production runs is therefore set by the recipe -- change it
# deliberately, with verify_enricher.sh passing.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if command -v enricher >/dev/null 2>&1; then
    ENRICHER="enricher"
elif [ -x "$SCRIPT_DIR/enricher" ]; then
    ENRICHER="$SCRIPT_DIR/enricher"
else
    echo "ERROR: 'enricher' binary not found (build it via Makefile_conda)." >&2
    exit 1
fi

# enricher <vcf.gz> <chrom.fa> <argv3-prefix> <yes|no> <argv5-path-for-fakedir-name>
# Argument order matches enricher.py exactly.
"$ENRICHER" "$1" "$2" "$3" "yes" "$1"
