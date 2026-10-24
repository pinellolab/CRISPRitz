#!/usr/bin/env bash
#
# Test for the on-the-fly PAM filter (crispritz.py search --prune-pam / the
# CRISPRITZ_PAM_PRUNE env var / the 12th positional arg to searchTST).
#
# The PAM filter lets a PAM-less index (e.g. NNN, 3' / NNNN, 5') serve a more
# restrictive search PAM (e.g. NGG / NRG / NRR / NGA / NAG / NNNN / TTTV) by
# skipping candidates whose concrete flanking PAM does not match the search
# PAM, BEFORE they are emitted. The result must be identical to searching a
# dedicated index built for that restrictive PAM.
#
# This script drives the compiled C++ binaries directly (no conda / crispritz.py),
# matching the other test/ scripts and CI, and self-checks every property:
#
#   1. DEFAULT no-op          : NNN index + NGG PAM with the filter OFF emits the
#                               SAME rows as an NGG-PAM search on that NNN index
#                               with no filter -> byte-identical to upstream.
#   2. EQUIVALENCE (3' PAM)    : NNN index + NGG PAM + filter == dedicated NGG
#                               index, as a SET (same targets/coords/scores) --
#                               the committed cas9 golden.
#   3. EQUIVALENCE (5' PAM)    : NNNN index + TTTV PAM + filter == dedicated TTTV
#                               (Cas12a) index, as a SET -- the cas12a golden.
#
# NOTE on ordering: the filtered PAM-less search produces the SAME target SET as
# a dedicated-PAM index, but the row ORDER within the file can differ, because
# the two indices have different tree shapes / cluster chains (the dedicated
# index never stored the excluded-PAM leaves). The content (targets, coords,
# mismatches, bulge sizes, direction) is identical. The equivalence checks
# therefore compare SORTED rows. Downstream CRISPRme+ clusters/sorts targets, so
# file order is not semantically meaningful. The DEFAULT no-op check (same index,
# filter off vs the upstream invocation) IS asserted byte-for-byte.
#   4. SUBSET + PAM-correctness: for NGG/NRG/NRR/NGA/NAG/NNNN the filtered rows
#                               are a SUBSET of the unfiltered NNN rows AND every
#                               kept row's emitted DNA PAM satisfies the search PAM
#                               AND every dropped row's PAM does NOT (exact set).
#   5. BULGE path              : the same SUBSET + PAM-correctness holds with
#                               mismatches + DNA + RNA bulges (saveIndices2).
#   6. LENGTH GUARD            : a 3-nt NNN index searched with a 4-nt PAM and the
#                               filter on exits non-zero with a clear error.
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
FIX_DIR="${SCRIPT_DIR}/fixtures"
GOLDEN_DIR="${FIX_DIR}/golden"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"

GENOME_FA="${FIX_DIR}/genome/chr1.fa"
GUIDE_CAS9="${FIX_DIR}/guide_cas9.txt"      # GACTGCATTGCATGACTGCANNN   (3' PAM)
GUIDE_CAS12A="${FIX_DIR}/guide_cas12a.txt"  # NNNNCTTGTCAAGGCTATTGGTCA  (5' PAM)

# Determinism: single-threaded search so target ordering is stable.
export OMP_NUM_THREADS=1

echo "=== CRISPRitz on-the-fly PAM-filter test ==="
echo "genome   : ${GENOME_FA}"
echo "compiler : ${CXX} ${CXXFLAGS}"

for f in "${GENOME_FA}" "${GUIDE_CAS9}" "${GUIDE_CAS12A}" \
         "${GOLDEN_DIR}/cas9.targets.txt" "${GOLDEN_DIR}/cas12a.targets.txt"; do
  if [[ ! -s "${f}" ]]; then echo "ERROR: missing/empty fixture: ${f}"; exit 2; fi
done

WORKDIR="$(mktemp -d)"
cleanup() { if [[ -z "${KEEP_WORKDIR:-}" ]]; then rm -rf "${WORKDIR}"; fi; }
trap cleanup EXIT
echo "workdir  : ${WORKDIR}"

# --- build binaries (unless prebuilt) ---------------------------------------
BUILDTST_BIN="${BUILDTST:-${WORKDIR}/buildTST}"
SEARCHTST_BIN="${SEARCHTST:-${WORKDIR}/searchTST}"

if [[ -z "${BUILDTST:-}" ]]; then
  echo "--- compiling buildTST ---"
  # shellcheck disable=SC2086
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC_DIR}/mainParallel.cpp" -o "${BUILDTST_BIN}"
fi
if [[ -z "${SEARCHTST:-}" ]]; then
  echo "--- compiling searchTST ---"
  # shellcheck disable=SC2086
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} \
    "${SRC_DIR}/searchOnTST.cpp" "${SRC_DIR}/detailedOutput.cpp" "${SRC_DIR}/convert.cpp" \
    -I "${SRC_DIR}/include" -o "${SEARCHTST_BIN}"
fi

MAX_BULGES=2

# --- PAM fixtures (written here so the filter PAMs are self-contained) ------
# 3' PAMs (pam_at_start = false): 20 N for the guide + a 3-nt PAM, length "3".
PAMDIR="${WORKDIR}/pams"; mkdir -p "${PAMDIR}"
printf 'NNNNNNNNNNNNNNNNNNNNNNN 3\n' > "${PAMDIR}/NNN.txt"   # PAM-less (3' index base)
printf 'NNNNNNNNNNNNNNNNNNNNNGG 3\n' > "${PAMDIR}/NGG.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNRG 3\n' > "${PAMDIR}/NRG.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNRR 3\n' > "${PAMDIR}/NRR.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNGA 3\n' > "${PAMDIR}/NGA.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNAG 3\n' > "${PAMDIR}/NAG.txt"
printf 'NNNNNNNNNNNNNNNNNNNNNNNN 4\n' > "${PAMDIR}/NNNN.txt" # 4-nt 3' PAM (length-guard)
# 5' PAMs (pam_at_start = true, encoded by the NEGATIVE PAM length): 4-nt PAM.
printf 'NNNNNNNNNNNNNNNNNNNNNNNN -4\n' > "${PAMDIR}/NNNN5.txt" # PAM-less 5' index base
printf 'TTTVNNNNNNNNNNNNNNNNNNNN -4\n' > "${PAMDIR}/TTTV.txt"

# --- helpers ----------------------------------------------------------------
build_index() { # <pam_file> <index_dir>
  local pam="$1" dir="$2"
  mkdir -p "${dir}"
  ( cd "${dir}" && "${BUILDTST_BIN}" "${GENOME_FA}" "${pam}" 1 "${MAX_BULGES}" >/dev/null )
  local n; n="$(find "${dir}" -maxdepth 1 -name '*.bin' | wc -l | tr -d ' ')"
  if [[ "${n}" -lt 1 ]]; then echo "ERROR: buildTST produced no .bin for ${pam}"; exit 1; fi
}

# run_search <index_dir> <guide> <pam_file> <mm> <bDNA> <bRNA> <prune 0|1> <result>
run_search() {
  local idx="$1" guide="$2" pam="$3" mm="$4" bdna="$5" brna="$6" prune="$7" res="$8"
  "${SEARCHTST_BIN}" "${idx}/" "${guide}" "${mm}" "${bdna}" "${brna}" \
    "${pam}" "${res}" t 1 "${MAX_BULGES}" -1 "${prune}" >/dev/null
}

rows() { tail -n +2 "$1" | grep -c . || true; }

# Compare two targets files as SETS (sorted data rows, header excluded).
same_set() { # <a> <b>
  diff <(tail -n +2 "$1" | sort) <(tail -n +2 "$2" | sort) >/dev/null
}

# The PAM region a search uses, by orientation:
#  - 3' PAM (default): the LAST 3 (or 4) chars of the DNA column.
#  - 5' PAM (Cas12a) : the FIRST 4 chars of the DNA column.
# dna_pam <targets_file> <pamlen> <at_start 0|1>  -> prints one PAM per data row
dna_pam() {
  local f="$1" plen="$2" atstart="$3"
  if [[ "${atstart}" == "1" ]]; then
    tail -n +2 "${f}" | awk -F'\t' -v p="${plen}" '{s=$3; gsub(/-/,"",s); print toupper(substr(s,1,p))}'
  else
    tail -n +2 "${f}" | awk -F'\t' -v p="${plen}" '{s=$3; gsub(/-/,"",s); print toupper(substr(s,length(s)-p+1,p))}'
  fi
}

# IUPAC membership test: does concrete base B satisfy IUPAC code C?
pam_ok() { # <code> <base>
  local c="$1" b="$2"
  case "${c}" in
    A) [[ "${b}" == A ]] ;;
    C) [[ "${b}" == C ]] ;;
    G) [[ "${b}" == G ]] ;;
    T) [[ "${b}" == T ]] ;;
    N) true ;;
    R) [[ "${b}" == A || "${b}" == G ]] ;;
    Y) [[ "${b}" == C || "${b}" == T ]] ;;
    S) [[ "${b}" == G || "${b}" == C ]] ;;
    W) [[ "${b}" == A || "${b}" == T ]] ;;
    K) [[ "${b}" == G || "${b}" == T ]] ;;
    M) [[ "${b}" == A || "${b}" == C ]] ;;
    B) [[ "${b}" != A ]] ;;
    D) [[ "${b}" != C ]] ;;
    H) [[ "${b}" != G ]] ;;
    V) [[ "${b}" != T ]] ;;
    *) false ;;
  esac
}

# Assert: every PAM string in file $1 satisfies the IUPAC PAM pattern $2.
assert_all_pam_match() { # <pam_list_file> <pattern>
  local list="$1" pat="$2" line i ok
  while IFS= read -r line; do
    [[ -z "${line}" ]] && continue
    for ((i=0;i<${#pat};i++)); do
      if ! pam_ok "${pat:$i:1}" "${line:$i:1}"; then
        echo "ASSERT FAIL: kept PAM '${line}' violates pattern '${pat}' at pos $i"; return 1
      fi
    done
  done < "${list}"
  return 0
}

# Assert: every PAM string in file $1 does NOT satisfy the pattern (dropped set).
assert_none_pam_match() { # <pam_list_file> <pattern>
  local list="$1" pat="$2" line i bad
  while IFS= read -r line; do
    [[ -z "${line}" ]] && continue
    bad=0
    for ((i=0;i<${#pat};i++)); do
      if ! pam_ok "${pat:$i:1}" "${line:$i:1}"; then bad=1; break; fi
    done
    if [[ "${bad}" -eq 0 ]]; then
      echo "ASSERT FAIL: dropped PAM '${line}' actually MATCHES pattern '${pat}' (should have been kept)"; return 1
    fi
  done < "${list}"
  return 0
}

PASS=0; run_assert() { if "$@"; then PASS=$((PASS+1)); else echo "  -> FAILED: $*"; exit 1; fi; }

#############################################################################
echo ""
echo "########## build PAM-less indices ##########"
NNN_IDX="${WORKDIR}/idx_NNN";    build_index "${PAMDIR}/NNN.txt"   "${NNN_IDX}"
NNNN5_IDX="${WORKDIR}/idx_NNNN5"; build_index "${PAMDIR}/NNNN5.txt" "${NNNN5_IDX}"
# dedicated indices for the equivalence golden cross-check
NGG_IDX="${WORKDIR}/idx_NGG";    build_index "${FIX_DIR}/pam_cas9.txt"   "${NGG_IDX}"
TTTV_IDX="${WORKDIR}/idx_TTTV";  build_index "${FIX_DIR}/pam_cas12a.txt" "${TTTV_IDX}"
echo "indices built: NNN(3'), NNNN(5'), NGG(3' dedicated), TTTV(5' dedicated)"

#############################################################################
echo ""
echo "########## 1. DEFAULT no-op (filter OFF == upstream) ##########"
# NGG search on the NNN index with the filter OFF must be byte-identical whether
# we pass prune=0 explicitly or omit the 12th arg entirely (true upstream path).
run_search "${NNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/NGG.txt" 3 1 0 0 "${WORKDIR}/d_off"
"${SEARCHTST_BIN}" "${NNN_IDX}/" "${GUIDE_CAS9}" 3 1 0 "${PAMDIR}/NGG.txt" \
  "${WORKDIR}/d_upstream" t 1 "${MAX_BULGES}" -1 >/dev/null   # 11-arg upstream call
if diff -q "${WORKDIR}/d_off.targets.txt" "${WORKDIR}/d_upstream.targets.txt" >/dev/null; then
  echo "OK: prune=0 byte-identical to the 11-arg upstream invocation ($(rows "${WORKDIR}/d_off.targets.txt") rows)"; PASS=$((PASS+1))
else
  echo "ASSERT FAIL: filter OFF differs from upstream"; exit 1
fi

#############################################################################
echo ""
echo "########## 2. EQUIVALENCE 3' PAM: NNN+filter == dedicated NGG (golden) ##########"
run_search "${NNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/NGG.txt" 3 1 0 1 "${WORKDIR}/eq_ngg"
# Independently reconfirm the dedicated NGG index reproduces the committed golden.
run_search "${NGG_IDX}" "${GUIDE_CAS9}" "${FIX_DIR}/pam_cas9.txt" 3 1 0 0 "${WORKDIR}/ded_ngg"
# the dedicated NGG index still reproduces the committed golden byte-for-byte
if diff -u "${GOLDEN_DIR}/cas9.targets.txt" "${WORKDIR}/ded_ngg.targets.txt" >/dev/null; then
  echo "OK: dedicated NGG index reproduces the committed cas9 golden (byte-for-byte)"; PASS=$((PASS+1))
else echo "ASSERT FAIL: dedicated NGG index does not match cas9 golden"; exit 1; fi
# the NNN index + filter reproduces the SAME target SET (see NOTE on ordering)
if same_set "${GOLDEN_DIR}/cas9.targets.txt" "${WORKDIR}/eq_ngg.targets.txt"; then
  echo "OK: NNN index + NGG + --prune-pam == dedicated NGG golden as a SET ($(rows "${WORKDIR}/eq_ngg.targets.txt") rows, identical targets/coords/scores)"; PASS=$((PASS+1))
else
  echo "ASSERT FAIL: NNN+filter != dedicated NGG golden (set difference):"
  diff <(tail -n +2 "${GOLDEN_DIR}/cas9.targets.txt" | sort) <(tail -n +2 "${WORKDIR}/eq_ngg.targets.txt" | sort)
  exit 1
fi

#############################################################################
echo ""
echo "########## 3. EQUIVALENCE 5' PAM (Cas12a TTTV): NNNN+filter == dedicated TTTV (golden) ##########"
run_search "${NNNN5_IDX}" "${GUIDE_CAS12A}" "${PAMDIR}/TTTV.txt" 3 1 0 1 "${WORKDIR}/eq_tttv"
run_search "${TTTV_IDX}" "${GUIDE_CAS12A}" "${FIX_DIR}/pam_cas12a.txt" 3 1 0 0 "${WORKDIR}/ded_tttv"
if diff -u "${GOLDEN_DIR}/cas12a.targets.txt" "${WORKDIR}/ded_tttv.targets.txt" >/dev/null; then
  echo "OK: dedicated TTTV index reproduces the committed cas12a golden (byte-for-byte)"; PASS=$((PASS+1))
else echo "ASSERT FAIL: dedicated TTTV index does not match cas12a golden"; exit 1; fi
if same_set "${GOLDEN_DIR}/cas12a.targets.txt" "${WORKDIR}/eq_tttv.targets.txt"; then
  echo "OK: NNNN(5') index + TTTV + --prune-pam == dedicated TTTV golden as a SET ($(rows "${WORKDIR}/eq_tttv.targets.txt") rows). pam_at_start (5' PAM) path CONFIRMED."; PASS=$((PASS+1))
else
  echo "ASSERT FAIL: NNNN+filter != dedicated TTTV golden (pam_at_start bug; set difference):"
  diff <(tail -n +2 "${GOLDEN_DIR}/cas12a.targets.txt" | sort) <(tail -n +2 "${WORKDIR}/eq_tttv.targets.txt" | sort)
  exit 1
fi

#############################################################################
echo ""
echo "########## 4. SUBSET + PAM-correctness (mismatch path, 3' PAMs) ##########"
# Baseline: unfiltered NNN search (emit everything). Each restrictive PAM's
# filtered output must be a subset AND exactly partition the baseline by PAM.
run_search "${NNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/NNN.txt" 3 0 0 0 "${WORKDIR}/base_mm"
BASE_N="$(rows "${WORKDIR}/base_mm.targets.txt")"
echo "baseline NNN (mm only): ${BASE_N} rows"
for P in NGG NRG NRR NGA NAG NNNN; do
  PATLEN=3; [[ "${P}" == NNNN ]] && PATLEN=4
  # NNNN is a 4-nt pattern -> needs a 4-nt index; skip against the 3-nt NNN here
  # (covered by the length guard in section 6). Use its 3-nt prefix meaning: run
  # only the 3-nt PAMs against the 3-nt NNN index.
  [[ "${P}" == NNNN ]] && continue
  run_search "${NNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/${P}.txt" 3 0 0 1 "${WORKDIR}/mm_${P}"
  # subset: every filtered row present in baseline
  if ! comm -23 <(tail -n +2 "${WORKDIR}/mm_${P}.targets.txt" | sort) \
                 <(tail -n +2 "${WORKDIR}/base_mm.targets.txt" | sort) | grep -q .; then
    : # empty diff = subset holds
  else
    echo "ASSERT FAIL: ${P} filtered rows are not a subset of the NNN baseline"; exit 1
  fi
  # partition: kept PAMs all match pattern; dropped PAMs all fail it
  PAT="${P}"
  dna_pam "${WORKDIR}/mm_${P}.targets.txt" "${PATLEN}" 0 > "${WORKDIR}/kept_${P}.pams"
  comm -23 <(tail -n +2 "${WORKDIR}/base_mm.targets.txt" | sort) \
           <(tail -n +2 "${WORKDIR}/mm_${P}.targets.txt" | sort) > "${WORKDIR}/dropped_${P}.rows"
  # convert dropped rows to PAMs
  awk -F'\t' -v p="${PATLEN}" '{s=$2; gsub(/-/,"",s); print toupper(substr(s,length(s)-p+1,p))}' "${WORKDIR}/dropped_${P}.rows" > "${WORKDIR}/dropped_${P}.pams"
  run_assert assert_all_pam_match  "${WORKDIR}/kept_${P}.pams"    "${PAT}"
  run_assert assert_none_pam_match "${WORKDIR}/dropped_${P}.pams" "${PAT}"
  echo "OK: ${P} subset + exact PAM partition ($(rows "${WORKDIR}/mm_${P}.targets.txt") kept / ${BASE_N} baseline)"
done

# NNNN (4-nt) exercised against a 4-nt PAM-less index to prove 4-nt filtering.
NNNN_IDX="${WORKDIR}/idx_NNNN4"; build_index "${PAMDIR}/NNNN.txt" "${NNNN_IDX}"
run_search "${NNNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/NNNN.txt" 3 0 0 1 "${WORKDIR}/mm_NNNN"
run_search "${NNNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/NNNN.txt" 3 0 0 0 "${WORKDIR}/mm_NNNN_off"
if diff -q "${WORKDIR}/mm_NNNN.targets.txt" "${WORKDIR}/mm_NNNN_off.targets.txt" >/dev/null; then
  echo "OK: NNNN (all-N 4-nt) filter keeps everything (degenerate PAM == no-op), $(rows "${WORKDIR}/mm_NNNN.targets.txt") rows"; PASS=$((PASS+1))
else echo "ASSERT FAIL: NNNN all-N filter dropped rows (should keep all)"; exit 1; fi

#############################################################################
echo ""
echo "########## 5. BULGE path (mm + DNA + RNA bulges, saveIndices2) ##########"
run_search "${NNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/NNN.txt" 3 1 1 0 "${WORKDIR}/base_bulge"
BASE_B="$(rows "${WORKDIR}/base_bulge.targets.txt")"
run_search "${NNN_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/NGG.txt" 3 1 1 1 "${WORKDIR}/bulge_NGG"
if comm -23 <(tail -n +2 "${WORKDIR}/bulge_NGG.targets.txt" | sort) \
             <(tail -n +2 "${WORKDIR}/base_bulge.targets.txt" | sort) | grep -q .; then
  echo "ASSERT FAIL: NGG bulge filtered rows are not a subset of the NNN bulge baseline"; exit 1
fi
dna_pam "${WORKDIR}/bulge_NGG.targets.txt" 3 0 > "${WORKDIR}/kept_bulge_NGG.pams"
run_assert assert_all_pam_match "${WORKDIR}/kept_bulge_NGG.pams" "NGG"
echo "OK: bulge path subset + PAM-correct ($(rows "${WORKDIR}/bulge_NGG.targets.txt") kept / ${BASE_B} baseline)"

#############################################################################
echo ""
echo "########## 6. LENGTH GUARD (search PAM len != index PAM len) ##########"
# 3-nt NNN index searched with a 4-nt PAM + filter ON must fail fast, non-zero.
set +e
"${SEARCHTST_BIN}" "${NNN_IDX}/" "${GUIDE_CAS9}" 3 0 0 "${PAMDIR}/NNNN.txt" \
  "${WORKDIR}/guard" t 1 "${MAX_BULGES}" -1 1 > "${WORKDIR}/guard.log" 2>&1
RC=$?
set -e
if [[ "${RC}" -ne 0 ]] && grep -qi "length to match the index PAM length" "${WORKDIR}/guard.log"; then
  echo "OK: length-mismatch exits non-zero (rc=${RC}) with a clear error"; PASS=$((PASS+1))
else
  echo "ASSERT FAIL: expected non-zero exit + clear error on PAM length mismatch (rc=${RC})"
  cat "${WORKDIR}/guard.log"; exit 1
fi

echo ""
echo "=== PASS: on-the-fly PAM filter test (${PASS} assertions) ==="
