#!/usr/bin/env bash
#
# Golden-output regression test for the CRISPRitz TST search engine.
#
# For each configured Cas system {cas9, cas12a} it:
#   1. builds a TST index (buildTST) over the committed synthetic genome and
#      the system's PAM,
#   2. runs a search (searchTST) with mismatches + a DNA bulge (single-thread,
#      OMP_NUM_THREADS=1 for determinism),
#   3. diffs the produced ".targets.txt" byte-for-byte against the committed
#      golden reference in test/fixtures/golden/,
#   4. exits non-zero on ANY difference.
#
# This locks the search output byte-for-byte so a regression in the search
# engine (e.g. CRISPRme #105, or any change to mismatch/bulge/PAM handling)
# is caught on every push. It is the search-engine analog of CRISPRme's two
# precomputed brute-force benchmarks (Cas9 NGG + Cas12a TTTV).
#
# Dependency-light: drives the compiled C++ binaries directly (no conda /
# crispritz.py), so it runs in a plain gcc image, matching CI.
#
# Configuration via environment variables:
#   CXX        : compiler (default: g++).
#   CXXFLAGS   : extra build flags.
#   BUILDTST   : path to a prebuilt buildTST (skips building it).
#   SEARCHTST  : path to a prebuilt searchTST (skips building it).
#   MM         : mismatches (default 3).
#   BDNA       : DNA bulges (default 1).
#   MAX_BULGES : index/search max bulges (default 2).
#   SYSTEMS    : space-separated subset of {cas9 cas12a} (default: both).
#   KEEP_WORKDIR : if set, keep the temp work dir.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC_DIR="${REPO_ROOT}/sourceCode/CRISPR-Cas-Tree"
FIX_DIR="${SCRIPT_DIR}/fixtures"
GOLDEN_DIR="${FIX_DIR}/golden"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"
MM="${MM:-3}"
BDNA="${BDNA:-1}"
MAX_BULGES="${MAX_BULGES:-2}"
SYSTEMS="${SYSTEMS:-cas9 cas12a}"

GENOME_FA="${FIX_DIR}/genome/chr1.fa"

# Determinism: single-threaded search so target ordering is stable.
export OMP_NUM_THREADS=1

echo "=== CRISPRitz golden search regression test ==="
echo "genome     : ${GENOME_FA}"
echo "systems    : ${SYSTEMS}"
echo "mismatches : ${MM}    DNA bulges: ${BDNA}    max bulges: ${MAX_BULGES}"
echo "compiler   : ${CXX} ${CXXFLAGS}"
echo "threads    : OMP_NUM_THREADS=${OMP_NUM_THREADS}"

# --- work dir ---------------------------------------------------------------
WORKDIR="$(mktemp -d)"
cleanup() { if [[ -z "${KEEP_WORKDIR:-}" ]]; then rm -rf "${WORKDIR}"; fi; }
trap cleanup EXIT
echo "workdir    : ${WORKDIR}"

# --- build binaries (unless prebuilt) ---------------------------------------
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

# --- per-system fixture map -------------------------------------------------
# system -> "pam_file guide_file golden_file"
pam_for()    { case "$1" in cas9) echo "${FIX_DIR}/pam_cas9.txt";;    cas12a) echo "${FIX_DIR}/pam_cas12a.txt";;    esac; }
guide_for()  { case "$1" in cas9) echo "${FIX_DIR}/guide_cas9.txt";;  cas12a) echo "${FIX_DIR}/guide_cas12a.txt";;  esac; }
golden_for() { case "$1" in cas9) echo "${GOLDEN_DIR}/cas9.targets.txt";; cas12a) echo "${GOLDEN_DIR}/cas12a.targets.txt";; esac; }

run_one() {
  local sys="$1"
  local pam guide golden
  pam="$(pam_for "${sys}")"
  guide="$(guide_for "${sys}")"
  golden="$(golden_for "${sys}")"

  echo ""
  echo "################ ${sys} ################"
  for f in "${pam}" "${guide}" "${golden}" "${GENOME_FA}"; do
    if [[ ! -s "${f}" ]]; then echo "ERROR: missing/empty fixture: ${f}"; exit 2; fi
  done
  echo "PAM   : ${pam}"
  echo "guide : ${guide}"
  echo "golden: ${golden}"

  # --- build index ---
  local index_dir="${WORKDIR}/${sys}_index"
  mkdir -p "${index_dir}"
  echo "--- building index (buildTST) ---"
  ( cd "${index_dir}" && "${BUILDTST_BIN}" "${GENOME_FA}" "${pam}" 1 "${MAX_BULGES}" >/dev/null )
  local bin_count
  bin_count="$(find "${index_dir}" -maxdepth 1 -name '*.bin' | wc -l | tr -d ' ')"
  if [[ "${bin_count}" -lt 1 ]]; then echo "ERROR: buildTST produced no .bin index files for ${sys}"; exit 1; fi
  echo "index .bin files: ${bin_count}"

  # --- run search ---
  local result="${WORKDIR}/${sys}_result"
  echo "--- running search (searchTST) ---"
  "${SEARCHTST_BIN}" \
    "${index_dir}/" \
    "${guide}" \
    "${MM}" \
    "${BDNA}" \
    0 \
    "${pam}" \
    "${result}" \
    r \
    1 \
    "${MAX_BULGES}" >/dev/null

  local targets="${result}.targets.txt"
  if [[ ! -s "${targets}" ]]; then echo "ASSERT FAIL: empty/missing targets file for ${sys}: ${targets}"; exit 1; fi

  # --- diff against golden ---
  echo "--- diff vs golden ---"
  if diff -u "${golden}" "${targets}"; then
    local rows
    rows="$(tail -n +2 "${targets}" | grep -c . || true)"
    echo "OK: ${sys} search output matches golden byte-for-byte (${rows} target rows)."
  else
    echo ""
    echo "ASSERT FAIL: ${sys} search output DIFFERS from golden reference."
    echo "  golden : ${golden}"
    echo "  actual : ${targets}"
    echo "  (Set KEEP_WORKDIR=1 to inspect; if this change is intentional,"
    echo "   regenerate the golden and commit it.)"
    exit 1
  fi
}

for sys in ${SYSTEMS}; do
  case "${sys}" in
    cas9|cas12a) run_one "${sys}" ;;
    *) echo "ERROR: unknown system '${sys}' (expected cas9 or cas12a)"; exit 2 ;;
  esac
done

echo ""
echo "=== PASS: all golden search outputs match byte-for-byte ==="
