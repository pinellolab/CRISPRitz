#!/usr/bin/env bash
#
# Regenerates test/fixtures/format1/: small format-1 indexes built by the CRISPRitz
# 2.9.1 buildTST, and what the 2.9.1 searchTST reports on them. Not run by CI (CI has
# no 2.9.1 binaries); run_format1_test.sh reads what this wrote. Kept so the fixture's
# provenance is reproducible:
#
#   OLD_BUILDTST=/path/to/2.9.1/buildTST OLD_SEARCHTST=/path/to/2.9.1/searchTST \
#     bash test/format1/make_fixture.sh
#
# The 2.9.1 binaries used for the committed fixture: buildTST md5
# 24f13ac4f581b990de7ec306c22ff292, searchTST md5 ef5f6e741a0827bbf7efbb09e6895580
# (compiled from tag v2.9.1 with g++ 14.2 -O3 on arm64; identical to the binaries in
# the crisprme:rc281 image).
#
# The contigs, each on its own (index 0 is per contig):
#   f1  interior sites on both strands (controls)
#   f2  a 3' PAM on the minus strand at base 0: index 0, a real PAM-left window, which
#       2.9.1 reports on the + strand at base max_bulges (IDX-14)
#   f3  a plus-strand window starting at 0 whose bases read as a minus-strand site that
#       matches the guide, with no PAM there: the 2.9.1 builder stores it as index 0 from
#       the wrong bases, and 2.9.1 reports it as a perfect match on the + strand
#   f4  a 5' PAM on the plus strand at base 0 (index 0 built as the other shape) and
#       f6 a 5' minus-strand window at 0 (index 0, read correctly by 2.9.1: control)
#   f5  a 3' minus-strand PAM exactly K - L = 5 from the start: index 0 of the K=8
#       PAM-less index, which 2.9.1 reports on the + strand under --prune-pam (IDX-4)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIX="${SCRIPT_DIR}/../fixtures/format1"
: "${OLD_BUILDTST:?set OLD_BUILDTST to a CRISPRitz 2.9.1 buildTST}"
: "${OLD_SEARCHTST:?set OLD_SEARCHTST to a CRISPRitz 2.9.1 searchTST}"
export OMP_NUM_THREADS=1

mkdir -p "${FIX}"
python3 - "${FIX}" <<'PY'
import os, random, sys
fix = sys.argv[1]
rng = random.Random(291)
G = "GACTTAGCCTAAGGTCTACG"
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
contigs = {
    "f1": rnd(60) + G + "TGG" + rnd(40) + rc(G + "AGG") + rnd(60),
    "f2": rc(G + "TGG") + rnd(80),
    "f3": "AAA" + rc(G) + "GG" + rnd(80),
    "f4": "TTTA" + G + rnd(80),
    "f5": rnd(5) + rc(G + "TGG") + rnd(80),
    "f6": rnd(2) + rc(G) + rc("TTTA") + rnd(80),
}
os.makedirs(os.path.join(fix, "genome"), exist_ok=True)
for name, seq in contigs.items():
    with open(os.path.join(fix, "genome", name + ".fa"), "w") as fh:
        fh.write(">%s\n%s\n" % (name, seq))
PY

N20="NNNNNNNNNNNNNNNNNNNN"
printf '%sNGG 3\n' "${N20}"        > "${FIX}/pam_NGG.txt"
printf '%sNNNNNNNN 8\n' "${N20}"   > "${FIX}/pam_N8.txt"
printf 'TTTV%s -4\n' "${N20}"      > "${FIX}/pam_TTTV.txt"
printf 'GACTTAGCCTAAGGTCTACGNNN\n' > "${FIX}/guide3.txt"
printf 'GACTTAGCCTAAGGTCTACGNNNNNNNN\n' > "${FIX}/guide8.txt"
printf 'NNNNGACTTAGCCTAAGGTCTACG\n' > "${FIX}/guide5.txt"

build() {  # <index dir name> <pam file>
  mkdir -p "${FIX}/$1"
  for fa in "${FIX}"/genome/*.fa; do ( cd "${FIX}/$1" && "${OLD_BUILDTST}" "${fa}" "$2" 1 2 > /dev/null ); done
}
build idx_NGG_2  "${FIX}/pam_NGG.txt"
build idx_N8_2   "${FIX}/pam_N8.txt"
build idx_TTTV_2 "${FIX}/pam_TTTV.txt"

mkdir -p "${FIX}/expected_291"
# the searches run_format1_test.sh repeats: <tag> <index> <guide> <pam> <mm> <bDNA> <bRNA> <prune> <index-0 contigs>
while read -r tag idx guide pam mm bd br prune _index0; do
  "${OLD_SEARCHTST}" "${FIX}/${idx}/" "${FIX}/${guide}" "${mm}" "${bd}" "${br}" "${FIX}/${pam}" \
    "${FIX}/expected_291/${tag}" r 1 2 -1 "${prune}" > /dev/null
done < "${SCRIPT_DIR}/searches.txt"
ls -la "${FIX}"/*
