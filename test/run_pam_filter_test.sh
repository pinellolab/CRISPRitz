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

# A wider-flank index legitimately emits a WIDER DNA/crRNA PAM region than a
# dedicated L-flank index (e.g. an 8-flank NGG search reports "...GGGACGGC" where
# a 3-flank NGG index reports "...GGG"). The off-target LOCUS is identical. To
# compare an index with stored flank width K against a dedicated L-flank index,
# trim the trailing (K-L) flank characters from the crRNA (col 2) and DNA (col 3)
# columns of the K-flank file, then compare the SET. The flank region never
# contains bulge dashes, so a plain trailing-char trim is exact.
#   same_locus_set_3prime <wide_file> <dedicated_file> <K> <L>
same_locus_set_3prime() {
  local wide="$1" ded="$2" K="$3" L="$4" trim=$(( $3 - $4 ))
  diff <(tail -n +2 "${wide}" | awk -F'\t' -v OFS='\t' -v t="${trim}" \
           '{$2=substr($2,1,length($2)-t); $3=substr($3,1,length($3)-t); print}' | sort) \
       <(tail -n +2 "${ded}"  | sort) >/dev/null
}
# 5' variant: the wider flank is at the FRONT, so (a) trim the leading (K-L) chars
# from crRNA (col 2) and DNA (col 3), and (b) a 5' match's reported START coordinate
# includes the flank, so a K-flank match starts (K-L) bp LEFT of the equivalent
# L-flank match -- add (K-L) back to the Position (col 5) and Cluster Position
# (col 6) to recover the dedicated-index coordinate. (3' PAMs don't shift because
# the match start is the protospacer end-independent left edge of the protospacer.)
same_locus_set_5prime() {
  local wide="$1" ded="$2" K="$3" L="$4" trim=$(( $3 - $4 ))
  diff <(tail -n +2 "${wide}" | awk -F'\t' -v OFS='\t' -v t="${trim}" \
           '{$2=substr($2,t+1); $3=substr($3,t+1); $5=$5+t; $6=$6+t; print}' | sort) \
       <(tail -n +2 "${ded}"  | sort) >/dev/null
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
echo "########## 6. LENGTH GUARD (search PAM len > index flank width -> hard error) ##########"
# 3-nt NNN index searched with a 4-nt PAM + filter ON must fail fast, non-zero:
# the index never stored a 4th flank base, so a 4-nt PAM cannot be served (L > K).
set +e
"${SEARCHTST_BIN}" "${NNN_IDX}/" "${GUIDE_CAS9}" 3 0 0 "${PAMDIR}/NNNN.txt" \
  "${WORKDIR}/guard" t 1 "${MAX_BULGES}" -1 1 > "${WORKDIR}/guard.log" 2>&1
RC=$?
set -e
if [[ "${RC}" -ne 0 ]] && grep -qi "LONGER than the index" "${WORKDIR}/guard.log"; then
  echo "OK: L>K (search PAM longer than index flank) exits non-zero (rc=${RC}) with a clear error"; PASS=$((PASS+1))
else
  echo "ASSERT FAIL: expected non-zero exit + clear error on L>K PAM length (rc=${RC})"
  cat "${WORKDIR}/guard.log"; exit 1
fi

#############################################################################
# ===== CONFIGURABLE GENEROUS FLANK: one wide PAM-less index serves ANY =====
# ===== PAM up to the stored flank width K, shorter PAMs N-padded.      =====
#############################################################################
# Build ONE generous K=8 3' PAM-less index (20 guide-N + 8 flank-N, length "8")
# and serve 3-, 4- and 6-nt PAMs off it via N-padding. The win: no per-PAM index.
printf 'NNNNNNNNNNNNNNNNNNNNNNNNNNNN 8\n'  > "${PAMDIR}/N8.txt"    # PAM-less K=8 build/base
printf 'NNNNNNNNNNNNNNNNNNNNNGGNNNNN 8\n'  > "${PAMDIR}/NGG8.txt"  # NGG as an explicit 8-wide PAM (control)

echo ""
echo "########## 7. GENEROUS FLANK build (K=8 PAM-less 3' index) ##########"
N8_IDX="${WORKDIR}/idx_N8"; build_index "${PAMDIR}/N8.txt" "${N8_IDX}"
# confirm the stored flank width really is 8 (the .bin prefix)
PREFIX8="$(find "${N8_IDX}" -maxdepth 1 -name '*.bin' -printf '%f\n' | head -1 | cut -d_ -f1)"
if [[ "${#PREFIX8}" -eq 8 ]]; then
  echo "OK: K=8 index built; stored flank prefix='${PREFIX8}' (8 nt)"; PASS=$((PASS+1))
else echo "ASSERT FAIL: expected an 8-nt flank prefix, got '${PREFIX8}' (${#PREFIX8} nt)"; exit 1; fi

echo ""
echo "########## 8. N-PADDING: short PAMs (<=8) served off the K=8 index ##########"
# baseline: the full unfiltered K=8 PAM-less search (emit everything), mm only.
run_search "${N8_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/N8.txt" 3 0 0 0 "${WORKDIR}/base8"
BASE8="$(rows "${WORKDIR}/base8.targets.txt")"
echo "baseline K=8 PAM-less (mm only): ${BASE8} rows"
# proximal flank of a 3' DNA column: the 8-wide flank is the LAST 8 chars; a
# length-L 3' PAM constrains its PROXIMAL L positions, i.e. the first L of that
# 8-char window = chars [len-8+1 .. len-8+L] (1-based). (NOT the last L -- those
# are distal and were N-padded.)
prox8_pam() { # <targets_file> <L>
  local f="$1" L="$2"
  tail -n +2 "${f}" | awk -F'\t' -v L="${L}" '{s=$3; gsub(/-/,"",s); print toupper(substr(s,length(s)-8+1,L))}'
}
# Each of these short 3' PAMs (length < 8) must N-pad up to 8 and still partition
# the baseline exactly by their PROXIMAL-L relevant flank letters.
for spec in "NGG:3" "NRG:3" "NRR:3" "NNNN:4" "NGGNRN:6"; do
  P="${spec%%:*}"; L="${spec##*:}"
  # write an L-length 3' PAM file: 20 guide-N + the L PAM letters, number L. The
  # binary N-pads up to the index flank width (8) on the distal/end side itself.
  printf 'NNNNNNNNNNNNNNNNNNNN%s %d\n' "${P}" "${L}" > "${PAMDIR}/gf_${P}.txt"
  run_search "${N8_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/gf_${P}.txt" 3 0 0 1 "${WORKDIR}/gf_${P}"
  # subset of the K=8 baseline?
  if comm -23 <(tail -n +2 "${WORKDIR}/gf_${P}.targets.txt" | sort) \
               <(tail -n +2 "${WORKDIR}/base8.targets.txt" | sort) | grep -q .; then
    echo "ASSERT FAIL: ${P} (N-padded) rows are not a subset of the K=8 baseline"; exit 1
  fi
  # exact partition by the PROXIMAL L flank positions
  prox8_pam "${WORKDIR}/gf_${P}.targets.txt" "${L}" > "${WORKDIR}/gfkept_${P}.pams"
  comm -23 <(tail -n +2 "${WORKDIR}/base8.targets.txt" | sort) \
           <(tail -n +2 "${WORKDIR}/gf_${P}.targets.txt" | sort) > "${WORKDIR}/gfdrop_${P}.rows"
  awk -F'\t' -v L="${L}" '{s=$2; gsub(/-/,"",s); print toupper(substr(s,length(s)-8+1,L))}' "${WORKDIR}/gfdrop_${P}.rows" > "${WORKDIR}/gfdrop_${P}.pams"
  run_assert assert_all_pam_match  "${WORKDIR}/gfkept_${P}.pams" "${P}"
  run_assert assert_none_pam_match "${WORKDIR}/gfdrop_${P}.pams" "${P}"
  echo "OK: ${P} (len ${L}) N-padded to K=8 -> subset + exact PAM partition ($(rows "${WORKDIR}/gf_${P}.targets.txt") kept / ${BASE8} baseline)"
done

echo ""
echo "########## 9. EQUIVALENCE: NGG on the K=8 flank == NGG on a dedicated K=3 index ##########"
# The definitive proof the N-padding orientation is correct: searching NGG off the
# GENEROUS 8-flank PAM-less index yields the SAME target SET as the dedicated
# 3-flank NGG index (the committed cas9 golden). mm + bulges, to exercise both
# emit paths (saveIndices + saveIndices2).
run_search "${N8_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/gf_NGG.txt" 3 1 0 1 "${WORKDIR}/gf_ngg_vs"
# The 8-flank index emits an 8-char flank; trim the distal 5 (= K-L = 8-3) chars
# off the crRNA/DNA columns, then the off-target LOCI must be the SAME SET as the
# dedicated 3-flank NGG index (the committed cas9 golden).
if same_locus_set_3prime "${WORKDIR}/gf_ngg_vs.targets.txt" "${GOLDEN_DIR}/cas9.targets.txt" 8 3; then
  echo "OK: 8-flank NGG (N-padded) == dedicated 3-flank NGG golden, SAME LOCI ($(rows "${WORKDIR}/gf_ngg_vs.targets.txt") rows; distal 5 flank bases trimmed). N-padding orientation (3'/end) CONFIRMED."; PASS=$((PASS+1))
else
  echo "ASSERT FAIL: 8-flank NGG != dedicated 3-flank NGG golden loci (N-pad orientation bug); set difference:"
  diff <(tail -n +2 "${WORKDIR}/gf_ngg_vs.targets.txt" | awk -F'\t' -v OFS='\t' '{$2=substr($2,1,length($2)-5); $3=substr($3,1,length($3)-5); print}' | sort) \
       <(tail -n +2 "${GOLDEN_DIR}/cas9.targets.txt" | sort)
  exit 1
fi

echo ""
echo "########## 10. GENEROUS FLANK 5' PAM: K=8 Cas12a index, TTTV N-padded ##########"
# 5' generous flank: 8-N PAM + 24-N guide region (negative length encodes 5').
printf 'NNNNNNNNNNNNNNNNNNNNNNNNNNNN -8\n' > "${PAMDIR}/N8_5.txt"
N8_5_IDX="${WORKDIR}/idx_N8_5"; build_index "${PAMDIR}/N8_5.txt" "${N8_5_IDX}"
# TTTV as a 4-nt 5' PAM (declared length 4), N-padded (on the FRONT) to the
# 8-flank index by the binary. The PAM FILE still declares the real 4-nt PAM:
# total length = 20-N protospacer + TTTV = 24, number -4 (so pamlen-pamlimit=20).
printf 'TTTVNNNNNNNNNNNNNNNNNNNN -4\n' > "${PAMDIR}/gf_TTTV.txt"
run_search "${N8_5_IDX}" "${GUIDE_CAS12A}" "${PAMDIR}/gf_TTTV.txt" 3 1 0 1 "${WORKDIR}/gf_tttv"
# 5' flank is at the FRONT; trim the distal 4 (= K-L = 8-4) leading chars, then
# the loci must match the dedicated 4-flank TTTV golden.
if same_locus_set_5prime "${WORKDIR}/gf_tttv.targets.txt" "${GOLDEN_DIR}/cas12a.targets.txt" 8 4; then
  echo "OK: 8-flank TTTV (5', front-N-padded) == dedicated 4-flank TTTV golden, SAME LOCI ($(rows "${WORKDIR}/gf_tttv.targets.txt") rows; distal 4 flank bases trimmed). 5' pam_at_start N-pad orientation CONFIRMED."; PASS=$((PASS+1))
else
  echo "ASSERT FAIL: 8-flank TTTV != dedicated TTTV golden loci (5' N-pad orientation bug); set difference:"
  diff <(tail -n +2 "${WORKDIR}/gf_tttv.targets.txt" | awk -F'\t' -v OFS='\t' '{$2=substr($2,5); $3=substr($3,5); $5=$5+4; $6=$6+4; print}' | sort) \
       <(tail -n +2 "${GOLDEN_DIR}/cas12a.targets.txt" | sort)
  exit 1
fi

echo ""
echo "########## 11. DEFAULT-OFF on the K=8 index is byte-identical (no regression) ##########"
# Filter OFF on the generous index == the 11-arg upstream invocation, byte-for-byte.
run_search "${N8_IDX}" "${GUIDE_CAS9}" "${PAMDIR}/N8.txt" 3 1 0 0 "${WORKDIR}/gf_off"
"${SEARCHTST_BIN}" "${N8_IDX}/" "${GUIDE_CAS9}" 3 1 0 "${PAMDIR}/N8.txt" \
  "${WORKDIR}/gf_upstream" t 1 "${MAX_BULGES}" -1 >/dev/null
if diff -q "${WORKDIR}/gf_off.targets.txt" "${WORKDIR}/gf_upstream.targets.txt" >/dev/null; then
  echo "OK: K=8 index, filter OFF byte-identical to the 11-arg upstream call ($(rows "${WORKDIR}/gf_off.targets.txt") rows)"; PASS=$((PASS+1))
else echo "ASSERT FAIL: K=8 default-off differs from upstream"; exit 1; fi

echo ""
echo "=== PASS: on-the-fly PAM filter test (${PASS} assertions) ==="
