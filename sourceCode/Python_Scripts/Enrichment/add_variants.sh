#!/usr/bin/env bash
# $1 is vcf file (Eg. variants_hg38/ALL.chr1.vcf.gz)
# $2 is chrname with file extension(Eg. chr1.fa)
# $3 is chrname (Eg. chr1)
#
# Invokes the compiled C++ `enricher` (a byte-identical port of enricher.py,
# the reference implementation kept alongside this script). The C++ binary is
# built by `make -f Makefile_conda crispritz` and installed on PATH; fall back
# to the copy sitting next to this script if it is not on PATH.

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
