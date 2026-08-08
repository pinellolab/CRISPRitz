#!/usr/bin/env bash
#
# Minimal end-to-end test for the CRISPRitz TST index/search binaries.
#
# It builds a TST index (buildTST) over a tiny committed synthetic genome and
# a PAM, then runs a search (searchTST) with mismatches + a DNA bulge, and
# asserts the search exits 0 and produces a non-empty targets file with the
# expected 10-column header.
#
# This script is deliberately dependency-light: it drives the compiled C++
# binaries directly (no conda / no crispritz.py), so it runs in a plain
# gcc image or under AddressSanitizer.
#
# Configuration via environment variables:
#   PAM        : which PAM fixture to use: "NGG" (default) or "WTN".
#   CXX        : compiler to use for building the binaries (default: g++).
#   CXXFLAGS   : extra flags appended to the build (e.g. ASan flags).
#   BUILDTST   : path to a prebuilt buildTST binary (skips building it).
#   SEARCHTST  : path to a prebuilt searchTST binary (skips building it).
#   MM         : number of mismatches (default 3).
#   BDNA       : number of DNA bulges (default 1).
#   KEEP_WORKDIR : if set, do not delete the temporary work dir.
#
set -euo pipefail

# --- locate paths -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC_DIR="${REPO_ROOT}/sourceCode/CRISPR-Cas-Tree"
FIX_DIR="${SCRIPT_DIR}/fixtures"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"
MM="${MM:-3}"
BDNA="${BDNA:-1}"
PAM="${PAM:-NGG}"

case "${PAM}" in
  NGG) PAM_FILE="${FIX_DIR}/pamNGG.txt" ;;
  WTN) PAM_FILE="${FIX_DIR}/pamWTN.txt" ;;
  *)   echo "ERROR: unknown PAM '${PAM}' (expected NGG or WTN)"; exit 2 ;;
esac

GUIDE_FILE="${FIX_DIR}/guide.txt"
GENOME_FA="${FIX_DIR}/genome/chr1.fa"

for f in "${PAM_FILE}" "${GUIDE_FILE}" "${GENOME_FA}"; do
  if [[ ! -s "${f}" ]]; then
    echo "ERROR: missing/empty fixture: ${f}"; exit 2
  fi
done

echo "=== CRISPRitz minimal test ==="
echo "PAM        : ${PAM} (${PAM_FILE})"
echo "guide      : ${GUIDE_FILE}"
echo "genome     : ${GENOME_FA}"
echo "mismatches : ${MM}    DNA bulges: ${BDNA}"
echo "compiler   : ${CXX} ${CXXFLAGS}"

# --- work dir ---------------------------------------------------------------
WORKDIR="$(mktemp -d)"
cleanup() {
  if [[ -z "${KEEP_WORKDIR:-}" ]]; then rm -rf "${WORKDIR}"; fi
}
trap cleanup EXIT
echo "workdir    : ${WORKDIR}"

# --- build binaries (unless prebuilt) ---------------------------------------
# The build invocation mirrors Makefile_conda. Note: mainParallel.cpp
# #include's pamGenerator.cpp, which is self-contained (aho-Corasick.cpp is
# intentionally not needed), and the search binary does NOT require boost.
BUILDTST_BIN="${BUILDTST:-${WORKDIR}/buildTST}"
SEARCHTST_BIN="${SEARCHTST:-${WORKDIR}/searchTST}"

if [[ -z "${BUILDTST:-}" ]]; then
  echo "--- compiling buildTST ---"
  # shellcheck disable=SC2086
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} \
    "${SRC_DIR}/mainParallel.cpp" \
    -o "${BUILDTST_BIN}"
fi

if [[ -z "${SEARCHTST:-}" ]]; then
  echo "--- compiling searchTST ---"
  # shellcheck disable=SC2086
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} \
    "${SRC_DIR}/searchOnTST.cpp" \
    "${SRC_DIR}/detailedOutput.cpp" \
    "${SRC_DIR}/convert.cpp" \
    -I "${SRC_DIR}/include" \
    -o "${SEARCHTST_BIN}"
fi

# --- build the index --------------------------------------------------------
# buildTST args: <genome.fa> <pam.txt> <threads> <max_bulges>
# It writes .bin files into the current working directory, so run it inside a
# dedicated index dir.
INDEX_DIR="${WORKDIR}/index"
mkdir -p "${INDEX_DIR}"
MAX_BULGES=2
echo "--- building index (buildTST) ---"
(
  cd "${INDEX_DIR}"
  "${BUILDTST_BIN}" "${GENOME_FA}" "${PAM_FILE}" 1 "${MAX_BULGES}"
)

BIN_COUNT="$(find "${INDEX_DIR}" -maxdepth 1 -name '*.bin' | wc -l | tr -d ' ')"
if [[ "${BIN_COUNT}" -lt 1 ]]; then
  echo "ERROR: buildTST produced no .bin index files"; exit 1
fi
echo "index .bin files: ${BIN_COUNT}"

# --- run the search ---------------------------------------------------------
# searchTST args:
#   <indexDir> <guideFile> <mm> <bDNA> <bRNA> <pamFile> <resultName> <r|p|t> <threads> <maxBulges>
RESULT="${WORKDIR}/result"
echo "--- running search (searchTST) ---"
"${SEARCHTST_BIN}" \
  "${INDEX_DIR}/" \
  "${GUIDE_FILE}" \
  "${MM}" \
  "${BDNA}" \
  0 \
  "${PAM_FILE}" \
  "${RESULT}" \
  r \
  1 \
  "${MAX_BULGES}"

# searchTST returning here means exit 0 (set -e would have aborted otherwise).
echo "searchTST completed with exit 0"

# --- assertions -------------------------------------------------------------
TARGETS="${RESULT}.targets.txt"

if [[ ! -f "${TARGETS}" ]]; then
  echo "ASSERT FAIL: targets file not created: ${TARGETS}"; exit 1
fi
if [[ ! -s "${TARGETS}" ]]; then
  echo "ASSERT FAIL: targets file is empty: ${TARGETS}"; exit 1
fi

EXPECTED_HEADER="#Bulge type	crRNA	DNA	Chromosome	Position	Cluster Position	Direction	Mismatches	Bulge Size	Total"
ACTUAL_HEADER="$(head -n 1 "${TARGETS}")"
if [[ "${ACTUAL_HEADER}" != "${EXPECTED_HEADER}" ]]; then
  echo "ASSERT FAIL: unexpected header."
  echo "  expected: ${EXPECTED_HEADER}"
  echo "  actual  : ${ACTUAL_HEADER}"
  exit 1
fi

NCOL="$(head -n 1 "${TARGETS}" | awk -F'\t' '{print NF}')"
if [[ "${NCOL}" -ne 10 ]]; then
  echo "ASSERT FAIL: expected 10 columns in header, got ${NCOL}"; exit 1
fi

DATA_ROWS="$(tail -n +2 "${TARGETS}" | grep -c . || true)"
echo "targets file: ${TARGETS}"
echo "data rows (excluding header): ${DATA_ROWS}"
echo "--- targets preview ---"
head -n 6 "${TARGETS}" || true

# For the NGG case we planted real on/near-target sites, so we expect >=1 hit.
# For degenerate/other PAMs we only require the run to complete and produce a
# well-formed file (the #105 regression is about NOT crashing).
if [[ "${PAM}" == "NGG" && "${DATA_ROWS}" -lt 1 ]]; then
  echo "ASSERT FAIL: expected at least one NGG target row, got ${DATA_ROWS}"; exit 1
fi

echo "=== PASS: search completed, well-formed non-empty targets file ==="
