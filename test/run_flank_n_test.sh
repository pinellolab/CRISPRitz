#!/usr/bin/env bash
#
# Test for a hard genomic 'N' in the stored PAM/flank region of a PAM-less index.
#
# A PAM-less index stores a wide flank (e.g. K=8) so that ONE index can serve any
# shorter PAM. A hard genomic 'N' converts to the bitset 0000, which intersects
# nothing, so an 'N' anywhere in that wide flank used to exclude the target even
# though a dedicated index built for a SHORTER PAM over the same guide region
# keeps it. The wider the stored flank the more targets were lost, which inverts
# the point of a PAM-less index.
#
# A target is now kept when the 'N' falls in the PAM/flank region. It is stored as
# the 0000 nibble, so the PAM actually searched decides it: within the real PAM
# width a stored 'N' still prunes (0000 intersects nothing), while the N-pad that
# widens the search PAM to the stored width is not tested at all. A hard 'N' in
# the GUIDE region is still discarded, because the tree has no encoding for 'N'
# in a guide.
#
# Properties asserted, on synthetic contigs where the 'N' offset is exact:
#
#   1. ORACLE      : the dedicated NRR index finds exactly the planted targets it
#                    should, on BOTH strands -- so the equivalence below cannot
#                    pass by both sides being equally wrong.
#   2. GUIDE REGION: a planted target with an 'N' inside the protospacer is in
#                    NEITHER index. Which end of the extracted window holds the
#                    PAM depends on the PAM side AND the strand, so getting it
#                    wrong lets a guide-region 'N' reach the tree -- which has no
#                    'N' encoding -- and silently corrupts the index.
#   3. EQUIVALENCE : PAM-less K=8 + prune NRR == dedicated NRR over a
#                    mismatch/DNA-bulge/RNA-bulge matrix.
#   4. ROUND TRIP  : searching the PAM-less index with its OWN all-N PAM emits the
#                    flank 'N' back at the planted offset, proving the 0000 nibble
#                    survives the writer and the loader (the nibble shift must
#                    depend on how many nucleotides are pending, not on their
#                    value, or a 0x0 loses its slot and shifts the rest).
#
# The 5'/front-PAM construction branches are covered by the Cas12a legs of
# run_pam_filter_test.sh and run_golden_test.sh.
#
# Equivalence compares SORTED rows: the two indices have different tree shapes, so
# row order can differ while content is identical (see run_pam_filter_test.sh).
#
# Configuration via environment variables:
#   CXX          : compiler (default: g++).
#   CXXFLAGS     : extra build flags.
#   BUILDTST     : path to a prebuilt buildTST (skips building it).
#   SEARCHTST    : path to a prebuilt searchTST (skips building it).
#   KEEP_WORKDIR : if set, keep the temp work dir.
#
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
[ -n "${KEEP_WORKDIR:-}" ] && echo "work dir: ${WORK}"

fail() { echo "FAIL: $*" >&2; exit 1; }

BT="${BUILDTST:-${WORK}/buildTST}"
if [ -z "${BUILDTST:-}" ]; then
  echo "--- compiling buildTST ---"
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC_DIR}/mainParallel.cpp" -o "${BT}"
fi
ST="${SEARCHTST:-${WORK}/searchTST}"
if [ -z "${SEARCHTST:-}" ]; then
  echo "--- compiling searchTST ---"
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC_DIR}/searchOnTST.cpp" \
    "${SRC_DIR}/detailedOutput.cpp" "${SRC_DIR}/convert.cpp" \
    -I "${SRC_DIR}/include" -o "${ST}"
fi

PROTO="ACTGAAATCTGTAAGCAGGC"     # 20nt protospacer; the guide matches it at 0 mm
PROTO_N="ACTGAAATCTNTAAGCAGGC"   # same, with a hard N inside the protospacer
PAD="$(printf 'T%.0s' $(seq 1 200))"

# contig  flank      proto  expected
#   flank offsets are counted outward from the protospacer, so for a 3' PAM
#   offsets 0,1,2 are what NRR constrains and 3..7 are the PAM-less index's extra
#   stored bases, which an NRR search must ignore. A leading 'm' plants the target
#   on the minus strand.
SPEC="\
pA AGGTTTTT P found
pB AGGTNTTT P found
pF AGGNTTTT P found
pG AGGTTTTN P found
pE NGGTTTTT P absent
pC ANGTTTTT P absent
pH AGNTTTTT P absent
pI ATTTTTTT P absent
pD AGGTTTTT N absent
mA AGGTTTTT P found
mB AGGTNTTT P found
mG AGGTTTTN P found
mE NGGTTTTT P absent
mC ANGTTTTT P absent
mD AGGTTTTT N absent"

revcomp() {
  printf '%s' "$1" | awk '{
    n = split($0, a, ""); s = "";
    for (i = n; i >= 1; i--) {
      c = a[i];
      if      (c == "A") s = s "T";
      else if (c == "T") s = s "A";
      else if (c == "C") s = s "G";
      else if (c == "G") s = s "C";
      else               s = s "N";
    }
    print s }'
}

GEN="${WORK}/genome"
mkdir -p "${GEN}"
: > "${WORK}/expected"
while read -r name flank protosel expect; do
  [ -z "${name}" ] && continue
  proto="${PROTO}"
  [ "${protosel}" = "N" ] && proto="${PROTO_N}"
  core="${proto}${flank}"                                  # 3' PAM: flank follows
  case "${name}" in m*) core="$(revcomp "${core}")" ;; esac
  printf '>%s\n%s%s%s\n' "${name}" "${PAD}" "${core}" "${PAD}" > "${GEN}/${name}.fa"
  if [ "${expect}" = "found" ]; then printf '%s\n' "${name}" >> "${WORK}/expected"; fi
done <<< "${SPEC}"
sort -o "${WORK}/expected" "${WORK}/expected"
echo "--- planted $(find "${GEN}" -name '*.fa' | wc -l | tr -d ' ') contigs, $(wc -l < "${WORK}/expected" | tr -d ' ') expected to be found ---"

printf 'NNNNNNNNNNNNNNNNNNNNNNNNNNNN 8\n' > "${WORK}/pamless.txt"   # K=8 PAM-less, 3'
printf 'NNNNNNNNNNNNNNNNNNNNNRR 3\n'      > "${WORK}/nrr.txt"       # dedicated NRR
printf '%sNNN\n'      "${PROTO}"          > "${WORK}/guide.txt"      # padded to NRR
printf '%sNNNNNNNN\n' "${PROTO}"          > "${WORK}/guide_wide.txt" # padded to K=8

build_index() { # <index_dir> <pam_file>
  local dir="$1" pam="$2" f n
  mkdir -p "${dir}"
  for f in "${GEN}"/*.fa; do
    ( cd "${dir}" && "${BT}" "${f}" "${pam}" 1 "${MAX_BULGES}" >/dev/null )
  done
  n="$(find "${dir}" -maxdepth 1 -name '*.bin' | wc -l | tr -d ' ')"
  if [ "${n}" -lt 1 ]; then fail "buildTST produced no .bin for ${pam}"; fi
}

# run_search <index_dir> <guide> <pam_file> <mm> <bDNA> <bRNA> <prune 0|1> <result>
run_search() {
  "${ST}" "$1/" "$2" "$4" "$5" "$6" "$3" "$8" t 1 "${MAX_BULGES}" -1 "$7" >/dev/null
}
contigs() { awk -F'\t' 'NR>1 {print $4}' "$1" | sort -u; }
body()    { tail -n +2 "$1" | sort; }

IDX_PL="${WORK}/idx_pamless"
IDX_DD="${WORK}/idx_nrr"
build_index "${IDX_PL}" "${WORK}/pamless.txt"
build_index "${IDX_DD}" "${WORK}/nrr.txt"

echo ""
echo "########## 1. ORACLE: the dedicated NRR index finds exactly the planted set ##########"
run_search "${IDX_DD}" "${WORK}/guide.txt" "${WORK}/nrr.txt" 0 0 0 0 "${WORK}/oracle"
got="$(contigs "${WORK}/oracle.targets.txt" | tr '\n' ' ')"
exp="$(tr '\n' ' ' < "${WORK}/expected")"
if [ "${got}" != "${exp}" ]; then
  echo "  found   : ${got}"
  echo "  expected: ${exp}"
  fail "the dedicated index does not find the planted target set"
fi
echo "  dedicated finds exactly [${got}]"

echo ""
echo "########## 2. GUIDE REGION: an N in the protospacer is in NEITHER index ##########"
run_search "${IDX_PL}" "${WORK}/guide.txt" "${WORK}/nrr.txt" 0 0 0 1 "${WORK}/pl0"
for f in "${WORK}/oracle.targets.txt" "${WORK}/pl0.targets.txt"; do
  if contigs "${f}" | grep -qxE 'pD|mD'; then
    fail "a target with an N in the protospacer reached the index (${f})"
  fi
done
echo "  protospacer-N contigs (pD, mD) absent from both indices"

echo ""
echo "########## 3. EQUIVALENCE: PAM-less K=8 + prune NRR == dedicated NRR ##########"
for cfg in "0 0 0" "2 0 0" "0 1 0" "0 0 1" "1 1 1"; do
  set -- ${cfg}
  tag="mm$1.b$2$3"
  run_search "${IDX_DD}" "${WORK}/guide.txt" "${WORK}/nrr.txt" "$1" "$2" "$3" 0 "${WORK}/D.${tag}"
  run_search "${IDX_PL}" "${WORK}/guide.txt" "${WORK}/nrr.txt" "$1" "$2" "$3" 1 "${WORK}/P.${tag}"
  nd="$(body "${WORK}/D.${tag}.targets.txt" | grep -c . || true)"
  np="$(body "${WORK}/P.${tag}.targets.txt" | grep -c . || true)"
  if [ "${nd}" -lt 1 ]; then fail "${tag}: the dedicated index returned no rows (vacuous)"; fi
  if ! diff -q <(body "${WORK}/D.${tag}.targets.txt") \
               <(body "${WORK}/P.${tag}.targets.txt") >/dev/null; then
    echo "  dedicated=${nd} pam-less=${np}"
    diff <(body "${WORK}/D.${tag}.targets.txt") \
         <(body "${WORK}/P.${tag}.targets.txt") | head -6
    fail "${tag}: PAM-less + prune does not equal dedicated"
  fi
  echo "  mm=$1 bDNA=$2 bRNA=$3: ${nd} rows, identical"
done

echo ""
echo "########## 4. ROUND TRIP: the flank N survives the writer and the loader ##########"
run_search "${IDX_PL}" "${WORK}/guide_wide.txt" "${WORK}/pamless.txt" 0 0 0 0 "${WORK}/wide"
dna="$(awk -F'\t' 'NR>1 && $4=="pB" {print $3; exit}' "${WORK}/wide.targets.txt")"
if [ -z "${dna}" ]; then fail "contig pB was not emitted by the own-PAM search"; fi
# pB plants its N at flank offset 4, i.e. position 20+4 of the emitted 28-mer
if [ "${dna:24:1}" != "N" ]; then
  fail "flank N not at the planted offset: got '${dna}' (expected N at index 24)"
fi
echo "  pB emitted as ${dna} -- N at the planted offset"

echo ""
echo "======================================"
echo "flank-N test: PASSED"
echo "======================================"
