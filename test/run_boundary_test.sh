#!/usr/bin/env bash
#
# Boundary test: sites at contig ends and N runs (index format 2, docs/INDEX_FORMAT.md).
#
# Until 2.9.1 an index lost, or misreported, sites near a contig end or an N run:
#   IDX-15  a PAM on a contig's last bases was in no index and missed by the
#           brute-force scan;
#   IDX-13  a site whose guide lay within max_bulges nt of a contig end or an N run
#           was dropped, even a perfect match (the whole pad + guide + PAM window had
#           to fit and be free of N);
#   IDX-4   with --prune-pam off a PAM-less index whose flank K is wider than the PAM
#           L, PAMs within K - L nt of a contig end were lost;
#   IDX-14  a leaf at index 0 was read as the wrong strand.
#
# make_genome.py plants perfect, 1-4-mismatch and 1-bulge sites at every distance
# d = 0..4 from a contig start, a contig end and an N run, on both strands, for a 3'
# PAM and a 5' PAM, with the PAM itself 0..9 bases from a contig end or an N run.
# Every index is built and searched through crispritz.py (index-genome, then search
# with -bDNA/-bRNA, with and without --prune-pam), and every search must equal, as a
# set of rows, oracle.py: an independent brute-force aligner that reads the FASTA
# directly, enumerates the sites of docs/INDEX_FORMAT.md section 3 and aligns with
# searchTST's alignment rules. Every 0-bulge search must also equal crispritz.py
# search without bulges, i.e. searchBruteForce, which uses no index.
#
# Configuration via environment variables:
#   CXX, CXXFLAGS     : compiler and extra flags (e.g. ASan)
#   BUILDTST, SEARCHTST, SEARCHBF : prebuilt binaries (skip compiling)
#   QUICK=1           : a smaller matrix (CI); default is the full matrix
#   KEEP_WORKDIR=1    : keep the work directory
#   KEEP_GOING=1      : report every differing search and fail at the end (to see
#                       what an older engine gets wrong); default stops at the first
#   PYTHON            : python3 to run crispritz.py and the oracle (default python3)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC="${REPO_ROOT}/sourceCode"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"
PYTHON="${PYTHON:-python3}"
QUICK="${QUICK:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

WORK="$(mktemp -d)"
cleanup() { [ -n "${KEEP_WORKDIR:-}" ] || rm -rf "${WORK}"; }
trap cleanup EXIT
[ -n "${KEEP_WORKDIR:-}" ] && echo "work dir: ${WORK}"
fail() { echo "FAIL: $*" >&2; exit 1; }
n_bad=0
differs() {  # a search that differs from its reference: stop, or count it
  if [ -n "${KEEP_GOING:-}" ]; then echo "  DIFFERS: $*"; n_bad=$((n_bad + 1)); else fail "$*"; fi
}

# --- binaries, laid out as an installed package: <root>/bin/crispritz.py and
#     <root>/opt/crispritz/<binary> (crispritz.py finds them at bin/../opt/crispritz)
ROOT="${WORK}/root"
mkdir -p "${ROOT}/bin" "${ROOT}/opt/crispritz"
cp "${REPO_ROOT}/crispritz.py" "${ROOT}/bin/crispritz.py"
OPT="${ROOT}/opt/crispritz"
if [ -n "${BUILDTST:-}" ]; then cp "${BUILDTST}" "${OPT}/buildTST"; else
  echo "--- compiling buildTST ---"
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC}/CRISPR-Cas-Tree/mainParallel.cpp" -o "${OPT}/buildTST"
fi
if [ -n "${SEARCHTST:-}" ]; then cp "${SEARCHTST}" "${OPT}/searchTST"; else
  echo "--- compiling searchTST ---"
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC}/CRISPR-Cas-Tree/searchOnTST.cpp" \
    "${SRC}/CRISPR-Cas-Tree/detailedOutput.cpp" "${SRC}/CRISPR-Cas-Tree/convert.cpp" \
    -I "${SRC}/CRISPR-Cas-Tree/include" -o "${OPT}/searchTST"
fi
if [ -n "${SEARCHBF:-}" ]; then cp "${SEARCHBF}" "${OPT}/searchBruteForce"; else
  echo "--- compiling searchBruteForce ---"
  ( cd "${SRC}/CRISPRofiler" && "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} main.cpp profiling.cpp \
      guide_searching.cpp pam_searching.cpp pre_computation.cpp reading.cpp analysis.cpp \
      -o "${OPT}/searchBruteForce" )
fi
CZ="${ROOT}/bin/crispritz.py"

# --- genomes, one for each PAM side
for side in 3 5; do
  "${PYTHON}" "${SCRIPT_DIR}/boundary/make_genome.py" "${WORK}/g${side}" "${side}" > /dev/null
done
echo "--- genomes: $(ls "${WORK}/g3/genome" | wc -l | tr -d ' ') contigs (3' PAM), $(ls "${WORK}/g5/genome" | wc -l | tr -d ' ') contigs (5' PAM) ---"

# --- PAM files: <name> <declared PAM line>
PAMS="${WORK}/pams"; mkdir -p "${PAMS}"
N20="NNNNNNNNNNNNNNNNNNNN"
mkpam() { printf '%s %s\n' "$2" "$3" > "${PAMS}/$1.txt"; }
mkpam NGG    "${N20}NGG" 3
mkpam NRG    "${N20}NRG" 3
mkpam NG     "${N20}NG" 2
mkpam NNGRRT "${N20}NNGRRT" 6
mkpam NNN    "${N20}NNN" 3
mkpam N8     "${N20}NNNNNNNN" 8
mkpam TTTV   "TTTV${N20}" -4
mkpam NNNN5  "NNNN${N20}" -4

pamless() { [ -z "$(sed 's/ .*//' "${PAMS}/$1.txt" | tr -d N)" ]; }  # an all-N PAM

GUIDE="GACTTAGCCTAAGGTCTACG"
GUIDE2="GACTTAGCCTAAGGTCTAGC"   # two mismatches next to a 3' PAM: a second guide
guides_for() {  # <pam name> -> guide file padded to that PAM's width
  local name=$1 line w side
  line="$(cat "${PAMS}/${name}.txt")"
  w="${line##* }"
  side=3; [ "${w#-}" != "${w}" ] && side=5 && w="${w#-}"
  local pad; pad="$(printf 'N%.0s' $(seq 1 "${w}"))"
  local f="${WORK}/guides_${name}.txt"
  if [ "${side}" = 3 ]; then printf '%s%s\n%s%s\n' "${GUIDE}" "${pad}" "${GUIDE2}" "${pad}" > "${f}"
  else printf '%s%s\n%s%s\n' "${pad}" "${GUIDE}" "${pad}" "${GUIDE2}" > "${f}"; fi
  echo "${f}"
}

# --- the matrix
# index: <side> <index PAM> <max bulges> | <search PAM>:<prune 0|1> ...
if [ -n "${QUICK}" ]; then
  INDEXES="
3 NGG 0 | NGG:0
3 NGG 2 | NGG:0
3 N8 2 | NGG:1 NG:1 NNGRRT:1 N8:0
5 TTTV 2 | TTTV:0
5 NNNN5 2 | TTTV:1 NNNN5:0"
  SEARCHES="0,0,0 4,0,0 1,1,0 1,0,1 1,2,2"
else
  INDEXES="
3 NGG 0 | NGG:0
3 NGG 1 | NGG:0
3 NGG 2 | NGG:0
3 NGG 3 | NGG:0
3 NRG 2 | NRG:0
3 NNN 2 | NGG:1 NRG:1 NG:1 NNN:0
3 N8 2 | NGG:1 NRG:1 NG:1 NNGRRT:1 N8:0
3 N8 1 | NGG:1 N8:0
5 TTTV 0 | TTTV:0
5 TTTV 2 | TTTV:0
5 NNNN5 2 | TTTV:1 NNNN5:0"
  SEARCHES="0,0,0 4,0,0 1,1,0 1,0,1 2,1,1 0,2,0 0,0,2 1,2,2 0,3,0 1,3,1"
fi

n_cmp=0; n_bf=0; n_rows=0
LIB="${WORK}/lib"; mkdir -p "${LIB}"
body() { tail -n +2 "$1" | LC_ALL=C sort; }

while IFS='|' read -r ispec sspecs; do
  [ -z "${ispec// /}" ] && continue
  set -- ${ispec}
  side=$1; ipam=$2; mb=$3
  genome="${WORK}/g${side}/genome"
  name="b${side}${ipam}${mb}"
  ( cd "${LIB}" && "${PYTHON}" "${CZ}" index-genome "${name}" "${genome}/" "${PAMS}/${ipam}.txt" \
      -bMax "${mb}" -th 1 > "${WORK}/${name}.index.log" 2>&1 ) \
    || fail "index-genome ${ipam} -bMax ${mb}: $(tail -3 "${WORK}/${name}.index.log")"
  idx="$(ls -d "${LIB}/genome_library/"*"_${mb}_${name}")"
  echo "--- index ${ipam}, -bMax ${mb}: $(ls "${idx}" | wc -l | tr -d ' ') .bin files"
  for sspec in ${sspecs}; do
    spam="${sspec%%:*}"; prune="${sspec##*:}"
    gfile="$(guides_for "${spam}")"
    flag=""; [ "${prune}" = 1 ] && flag="--prune-pam"
    for s in ${SEARCHES}; do
      IFS=, read -r mm bd br <<< "${s}"
      [ "${bd}" -le "${mb}" ] && [ "${br}" -le "${mb}" ] || continue
      tag="${name}.${spam}.p${prune}.mm${mm}.b${bd}${br}"
      ( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${idx}/" "${PAMS}/${spam}.txt" "${gfile}" "${tag}" \
          -mm "${mm}" -bDNA "${bd}" -bRNA "${br}" -th 1 -r ${flag} > "${tag}.log" 2> "${tag}.err" ) \
        || fail "${tag}: search exited non-zero: $(tail -3 "${WORK}/${tag}.err")"
      [ ! -s "${WORK}/${tag}.err" ] || fail "${tag}: search wrote to stderr: $(head -3 "${WORK}/${tag}.err")"
      oflag=""; [ "${prune}" = 1 ] && oflag="--prune"
      "${PYTHON}" "${SCRIPT_DIR}/boundary/oracle.py" --genome "${genome}" --guides "${gfile}" \
        --pam "${PAMS}/${spam}.txt" --index-pam "${PAMS}/${ipam}.txt" --mb "${mb}" \
        --mm "${mm}" --bdna "${bd}" --brna "${br}" ${oflag} > "${WORK}/${tag}.oracle"
      if ! diff <(body "${WORK}/${tag}.targets.txt") "${WORK}/${tag}.oracle" > "${WORK}/${tag}.diff"; then
        [ -n "${KEEP_GOING:-}" ] || { echo "  engine (<) vs oracle (>), first lines:"; head -12 "${WORK}/${tag}.diff"; }
        differs "${tag}: indexed search != oracle ($(grep -c '^<' "${WORK}/${tag}.diff" || true) extra, $(grep -c '^>' "${WORK}/${tag}.diff" || true) missing)"
      fi
      rows="$(wc -l < "${WORK}/${tag}.oracle" | tr -d ' ')"
      [ "${rows}" -gt 0 ] || fail "${tag}: no rows (vacuous comparison)"
      n_cmp=$((n_cmp + 1)); n_rows=$((n_rows + rows))
      # 0 bulges: the brute-force scan must agree too (it needs real PAM bases, so a
      # PAM-less index searched without --prune-pam, whose flank is unconstrained, is
      # left out)
      if [ "${bd}" = 0 ] && [ "${br}" = 0 ] && { [ "${prune}" = 1 ] || ! pamless "${ipam}"; }; then
        ( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${genome}/" "${PAMS}/${spam}.txt" "${gfile}" "${tag}.bf" \
            -mm "${mm}" -th 1 -r > "${tag}.bf.log" 2> "${tag}.bf.err" ) || fail "${tag}: brute-force search failed"
        [ ! -s "${WORK}/${tag}.bf.err" ] || fail "${tag}: brute-force search wrote to stderr"
        if ! diff <(body "${WORK}/${tag}.targets.txt") <(body "${WORK}/${tag}.bf.targets.txt") > "${WORK}/${tag}.bf.diff"; then
          [ -n "${KEEP_GOING:-}" ] || { echo "  indexed (<) vs brute force (>):"; head -12 "${WORK}/${tag}.bf.diff"; }
          differs "${tag}: indexed 0-bulge search != searchBruteForce ($(grep -c '^<' "${WORK}/${tag}.bf.diff" || true) only indexed, $(grep -c '^>' "${WORK}/${tag}.bf.diff" || true) only brute force)"
        fi
        n_bf=$((n_bf + 1))
      fi
    done
  done
  echo "    ok: every search of ${ipam} -bMax ${mb} == oracle (${n_cmp} comparisons so far)"
done <<< "${INDEXES}"

[ "${n_bad}" = 0 ] || fail "${n_bad} searches differ from their reference"
echo ""
echo "boundary test: PASSED (${n_cmp} indexed searches == oracle, ${n_rows} rows; ${n_bf} of them also == searchBruteForce)"
