#!/usr/bin/env bash
#
# Format-1 compatibility: indexes built by CRISPRitz 2.9.1 keep searching as with
# 2.9.1, but for the leaf at index 0 (IDX-14), which 2.9.1 read as the wrong shape.
#
# test/fixtures/format1/ holds small indexes built by the 2.9.1 buildTST and what the
# 2.9.1 searchTST reported on them (test/format1/make_fixture.sh says how). For each
# search in test/format1/searches.txt this searches the same index with this tree's
# searchTST and asserts:
#
#   1. every row of the contigs whose index-0 leaf does not produce a row (the last
#      column of searches.txt lists those that do) is byte-identical to 2.9.1's, in
#      the same order;
#   2. the rows of the listed contigs are exactly those the independent oracle
#      (test/boundary/oracle.py) lists for them: f2's minus-strand PAM at base 0 and
#      f5's minus-strand PAM K - L from the start under --prune-pam come out on the
#      minus strand at their real position instead of on the + strand, and f3's bogus
#      perfect match (a window the 2.9.1 builder stored from the wrong bases) is gone;
#   3. the whole output equals test/fixtures/format1/expected_292/, byte for byte.
#
# A format-1 index lacks the boundary sites of index format 2 (docs/INDEX_FORMAT.md);
# only a rebuild adds those (f2's site off the K=8 index, f4's 5' PAM at base 0), so
# the oracle comparison is limited to the listed contigs.
#
# Configuration: CXX, CXXFLAGS, SEARCHTST (prebuilt binary), PYTHON, KEEP_WORKDIR.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC_DIR="${REPO_ROOT}/sourceCode/CRISPR-Cas-Tree"
FIX="${SCRIPT_DIR}/fixtures/format1"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"
PYTHON="${PYTHON:-python3}"
export OMP_NUM_THREADS=1

WORK="$(mktemp -d)"
cleanup() { [ -n "${KEEP_WORKDIR:-}" ] || rm -rf "${WORK}"; }
trap cleanup EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

ST="${SEARCHTST:-${WORK}/searchTST}"
if [ -z "${SEARCHTST:-}" ]; then
  echo "--- compiling searchTST ---"
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC_DIR}/searchOnTST.cpp" \
    "${SRC_DIR}/detailedOutput.cpp" "${SRC_DIR}/convert.cpp" -I "${SRC_DIR}/include" -o "${ST}"
fi

# the fixture's .bin files are format 1: no CZI magic at the start
for f in "${FIX}"/idx_*/*.bin; do
  [ "$(head -c 4 "${f}" | od -An -tx1 | tr -d ' \n')" != "435a49fe" ] || fail "${f} is not a format-1 index"
done

n=0
while read -r tag idx guide pam mm bd br prune index0; do
  [ -z "${tag}" ] && continue
  INDEX0="$(echo "${index0}" | tr ',' '|')"   # contigs whose index-0 leaf yields a row
  "${ST}" "${FIX}/${idx}/" "${FIX}/${guide}" "${mm}" "${bd}" "${br}" "${FIX}/${pam}" \
    "${WORK}/${tag}" r 1 2 -1 "${prune}" > "${WORK}/${tag}.log" 2> "${WORK}/${tag}.err" \
    || fail "${tag}: searchTST exited non-zero: $(tail -2 "${WORK}/${tag}.err")"
  [ ! -s "${WORK}/${tag}.err" ] || fail "${tag}: searchTST wrote to stderr: $(head -2 "${WORK}/${tag}.err")"
  grep -q "INDEX FORMAT: 1 " "${WORK}/${tag}.log" || fail "${tag}: the index was not read as format 1"
  new="${WORK}/${tag}.targets.txt"
  old="${FIX}/expected_291/${tag}.targets.txt"

  # 1. the other contigs: byte-identical to 2.9.1, order included
  if ! diff <(awk -F'\t' -v re="^(${INDEX0})\$" '$4 !~ re' "${old}") \
            <(awk -F'\t' -v re="^(${INDEX0})\$" '$4 !~ re' "${new}") > "${WORK}/${tag}.d1"; then
    head "${WORK}/${tag}.d1"; fail "${tag}: rows away from index 0 differ from 2.9.1"
  fi
  kept="$(awk -F'\t' -v re="^(${INDEX0})\$" 'NR > 1 && $4 !~ re' "${new}" | wc -l | tr -d ' ')"

  # 2. the index-0 contigs: the oracle's rows for them
  oflag=""; [ "${prune}" = 1 ] && oflag="--prune"
  ipam="pam_NGG.txt"; case "${idx}" in idx_N8_2) ipam="pam_N8.txt";; idx_TTTV_2) ipam="pam_TTTV.txt";; esac
  "${PYTHON}" "${SCRIPT_DIR}/boundary/oracle.py" --genome "${FIX}/genome" --guides "${FIX}/${guide}" \
    --pam "${FIX}/${pam}" --index-pam "${FIX}/${ipam}" --mb 2 --mm "${mm}" --bdna "${bd}" --brna "${br}" \
    ${oflag} | awk -F'\t' -v re="^(${INDEX0})\$" '$4 ~ re' > "${WORK}/${tag}.oracle0"
  if ! diff <(awk -F'\t' -v re="^(${INDEX0})\$" 'NR > 1 && $4 ~ re' "${new}" | LC_ALL=C sort) \
            "${WORK}/${tag}.oracle0" > "${WORK}/${tag}.d2"; then
    head "${WORK}/${tag}.d2"; fail "${tag}: the index-0 contigs' rows differ from the oracle"
  fi
  fixed="$(wc -l < "${WORK}/${tag}.oracle0" | tr -d ' ')"
  was="$(awk -F'\t' -v re="^(${INDEX0})\$" 'NR > 1 && $4 ~ re' "${old}" | wc -l | tr -d ' ')"

  # 3. the whole output, against the committed expectation
  diff -q "${FIX}/expected_292/${tag}.targets.txt" "${new}" > /dev/null \
    || fail "${tag}: output differs from test/fixtures/format1/expected_292/${tag}.targets.txt"

  echo "  ok  ${tag}: ${kept} rows as 2.9.1; index-0 contigs: ${was} rows in 2.9.1 -> ${fixed} (== oracle)"
  n=$((n + 1))
done < "${SCRIPT_DIR}/format1/searches.txt"
[ "${n}" -gt 0 ] || fail "no searches run"

echo "format-1 compatibility test: PASSED (${n} searches)"
