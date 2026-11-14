#!/usr/bin/env bash
#
# buildTST on a contig too short to hold one target window.
#
# searchPAMonGenome's scan bounds are `genome_sequence.length() - k` on an unsigned
# length, so a sequence shorter than the window wraps them to ~2^32. For a PAM with
# a real constraint the scan then read past the sequence and segfaulted; for an
# all-N PAM -- whose scan reads nothing, since a PAM-less index must not reject a
# site on PAM content -- it ran ~2^31 iterations pushing two indices each, and was
# killed for memory after taking tens of GB. Such contigs are ordinary: the enricher
# writes a header-only fake<chr>.fa for every chromosome with no carried indel.
#
# A contig shorter than one guide plus one PAM base holds no target on either strand
# (index format 2 no longer needs the bulge pad to fit, docs/INDEX_FORMAT.md), so
# buildTST must finish at once, exit 0 and write no .bin -- while a contig exactly one
# target long, with no room for the pad, must now be indexed (until 2.9.1 the whole
# pad + guide + PAM window had to fit, which lost such targets: IDX-13), and a normal
# chromosome beside a header-only one must still be indexed, through crispritz.py
# index-genome as well.
#
# Each run is under a timeout and, on Linux, a 4 GB address-space cap, so a
# regression fails this test instead of exhausting the host's memory.
#
# Configuration: CXX (default g++), CXXFLAGS, BUILDTST (prebuilt binary),
# KEEP_WORKDIR.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC_DIR="${REPO_ROOT}/sourceCode/CRISPR-Cas-Tree"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"
MAX_BULGES=3

WORK="$(mktemp -d)"
cleanup() { [ -n "${KEEP_WORKDIR:-}" ] || rm -rf "${WORK}"; }
trap cleanup EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

BT="${BUILDTST:-${WORK}/buildTST}"
if [ -z "${BUILDTST:-}" ]; then
  echo "--- compiling buildTST ---"
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC_DIR}/mainParallel.cpp" -o "${BT}"
fi

if command -v timeout > /dev/null 2>&1; then TIMEOUT="timeout 60"
elif command -v gtimeout > /dev/null 2>&1; then TIMEOUT="gtimeout 60"
else TIMEOUT="perl -e alarm(shift);exec(@ARGV) 60"; fi
CAP=""
[ "$(uname -s)" = Linux ] && CAP="ulimit -v 4000000;"

# PAM files: 3' all-N at the two published widths, a 5' all-N, a dedicated 3' PAM
printf 'NNNNNNNNNNNNNNNNNNNNNNN 3\n'        > "${WORK}/nnn.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNNNNNNNN 8\n'   > "${WORK}/n8.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNNNN -4\n'      > "${WORK}/n4_5p.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNGG 3\n'        > "${WORK}/ngg.txt"

# contigs: header-only, 10 bp, and one guide long (20 bp): one base short of the
# smallest target, a guide plus one PAM base (an all-N PAM; a dedicated PAM needs
# guide + PAM)
mkdir -p "${WORK}/fa"
printf '>chrE\n'                            > "${WORK}/fa/chrE.fa"
printf '>chrS\nACGTACGTAC\n'                > "${WORK}/fa/chrS.fa"
printf '>chrW\nACGTTGCATCGATCGGATCC\n'      > "${WORK}/fa/chrW.fa"

bt() {  # <out dir> <fasta> <pam file>; sets RC and NBIN
  local dir=$1 fa=$2 pam=$3
  rm -rf "${dir}"; mkdir -p "${dir}"
  set +e
  ( cd "${dir}" && bash -c "${CAP} exec ${TIMEOUT} \"${BT}\" \"${fa}\" \"${pam}\" 1 ${MAX_BULGES}" > log.txt 2>&1 )
  RC=$?
  set -e
  NBIN=$(find "${dir}" -maxdepth 1 -name '*.bin' | wc -l | tr -d ' ')
}

for fa in chrE chrS chrW; do
  for pam in nnn n8 n4_5p ngg; do
    t0=$(date +%s)
    bt "${WORK}/out_${fa}_${pam}" "${WORK}/fa/${fa}.fa" "${WORK}/${pam}.txt"
    dt=$(( $(date +%s) - t0 ))
    case "${RC}" in
      0) ;;
      124|142) fail "${fa} / ${pam}: buildTST did not finish within 60 s (runaway scan)" ;;
      *) fail "${fa} / ${pam}: buildTST exited ${RC} (139 = segfault; 137 or bad_alloc = memory): $(tail -2 "${WORK}/out_${fa}_${pam}/log.txt")" ;;
    esac
    [ "${NBIN}" = 0 ] || fail "${fa} / ${pam}: a contig with no room for a target produced ${NBIN} .bin file(s)"
    echo "  ok  ${fa} / ${pam}: exit 0 in ${dt}s, no .bin"
  done
done

# a contig exactly one target long has no room for the pad and is indexed all the
# same: 23 bp holding guide + AGG for NGG, and for the all-N PAMs (any PAM fits) the
# same 23 bp and a 21-bp guide + one base
python3 - "${WORK}/fa/chrT.fa" "${WORK}/fa/chrU.fa" <<'PY'
import sys
open(sys.argv[1], "w").write(">chrT\nACGTTGCATCGATCGGATCCAGG\n")   # 23 bp: guide + AGG
open(sys.argv[2], "w").write(">chrU\nACGTTGCATCGATCGGATCCA\n")     # 21 bp: guide + 1
PY
for spec in "chrT ngg" "chrT nnn" "chrT n8" "chrT n4_5p" "chrU nnn" "chrU n8"; do
  set -- ${spec}
  bt "${WORK}/one_$1_$2" "${WORK}/fa/$1.fa" "${WORK}/$2.txt"
  [ "${RC}" = 0 ] || fail "$1 / $2: buildTST exited ${RC}"
  [ "${NBIN}" -ge 1 ] || fail "$1 / $2: a contig exactly one target long produced no .bin (the pad must not have to fit)"
  echo "  ok  $1 / $2: one-target contig indexed (${NBIN} .bin)"
done

# a real target still indexes: the same guard must not drop a window that fits
python3 - "${WORK}/fa_ok.fa" <<'PY'
import random, sys
random.seed(3)
g = [random.choice("ACGT") for _ in range(400)]
g[200:223] = list("ACTGAAATCTGTAAGCAGGCAGG")
open(sys.argv[1], "w").write(">chrOK\n" + "".join(g) + "\n")
PY
for pam in nnn n8 ngg; do
  bt "${WORK}/ok_${pam}" "${WORK}/fa_ok.fa" "${WORK}/${pam}.txt"
  [ "${RC}" = 0 ] || fail "chrOK / ${pam}: buildTST exited ${RC}"
  [ "${NBIN}" -ge 1 ] || fail "chrOK / ${pam}: a normal contig produced no .bin"
  echo "  ok  chrOK / ${pam}: indexed (${NBIN} .bin)"
done

# through the real CLI: index-genome over a normal chromosome plus a header-only one
ROOT="${WORK}/root"
mkdir -p "${ROOT}/bin" "${ROOT}/opt/crispritz" "${WORK}/genome"
cp "${REPO_ROOT}/crispritz.py" "${ROOT}/bin/crispritz.py"
cp "${BT}" "${ROOT}/opt/crispritz/buildTST"
cp "${WORK}/fa_ok.fa" "${WORK}/genome/chrOK.fa"
cp "${WORK}/fa/chrE.fa" "${WORK}/genome/chrE.fa"
mkdir -p "${WORK}/cli"
set +e
( cd "${WORK}/cli" && bash -c "${CAP} exec ${TIMEOUT} python3 \"${ROOT}/bin/crispritz.py\" index-genome tiny \"${WORK}/genome/\" \"${WORK}/nnn.txt\" -bMax ${MAX_BULGES} -th 1" > log.txt 2>&1 )
RC=$?
set -e
[ "${RC}" = 0 ] || fail "index-genome with a header-only chromosome exited ${RC}: $(tail -3 "${WORK}/cli/log.txt")"
IDX="${WORK}/cli/genome_library/NNN_${MAX_BULGES}_tiny"
n_ok=$(find "${IDX}" -maxdepth 1 -name '*chrOK*.bin' | wc -l | tr -d ' ')
n_e=$(find "${IDX}" -maxdepth 1 -name '*chrE*.bin' | wc -l | tr -d ' ')
[ "${n_ok}" -ge 1 ] || fail "index-genome did not index the normal chromosome beside the header-only one"
[ "${n_e}" = 0 ] || fail "index-genome wrote a .bin for the header-only chromosome"
echo "  ok  index-genome: exit 0, chrOK indexed (${n_ok} .bin), header-only chrE skipped"

echo "short-contig test: PASSED"
