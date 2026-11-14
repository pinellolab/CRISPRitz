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
# search without bulges, i.e. searchBruteForce, which uses no index. Every row must be
# its contig's own bases at its position (check_rows.py).
#
# Some searches use a guide longer or shorter than the index's (section 5.1: e.g.
# CRISPRme's 21-nt SaCas9 NNGRRT guides and 23-nt Cas12a TTTV guides off a 20-nt
# index). A longer guide must equal the brute force at 0 bulges; a shorter one finds
# a subset of it, missing only the sites whose index-length guide region does not
# fit, which oracle.py's model of a dedicated index for the shorter guide must equal.
# And --prune-pam off a PAM-less index must give the same rows as a dedicated index
# built for the same guide length and -bMax (PAIRS).
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
mkpam NNGRRT21 "N${N20}NNGRRT" 6   # CRISPRme's 21bp-NNGRRT-SaCas9.txt: 21-nt guides
mkpam NGG21  "N${N20}NGG" 3
mkpam NGG19  "${N20#N}NGG" 3
mkpam NNN    "${N20}NNN" 3
mkpam N8     "${N20}NNNNNNNN" 8
mkpam TTTV   "TTTV${N20}" -4
mkpam TTTV23 "TTTV${N20}NNN" -4     # CRISPRme's 23bp-TTTV-Cas12a.txt: 23-nt guides
mkpam TTTV19 "TTTV${N20#N}" -4
mkpam NNNN5  "NNNN${N20}" -4

pamless() { [ -z "$(sed 's/ .*//' "${PAMS}/$1.txt" | tr -d N)" ]; }  # an all-N PAM

GUIDE="GACTTAGCCTAAGGTCTACG"
GUIDE2="GACTTAGCCTAAGGTCTAGC"   # two mismatches next to a 3' PAM: a second guide
EXT="TCAG"                      # make_genome.py's E cases: the bases beyond the far end
glen() {  # <pam name> -> the guide length of that PAM file
  local line motif w
  line="$(cat "${PAMS}/$1.txt")"; motif="${line%% *}"; w="${line##* }"; w="${w#-}"
  echo $(( ${#motif} - w ))
}
fit() {  # <guide> <side> <length>: the 20-nt guide cut or extended at its far end
  local g=$1 side=$2 n=$3 x
  if [ "${n}" -le 20 ]; then
    if [ "${side}" = 3 ]; then echo "${g:$((20 - n))}"; else echo "${g:0:${n}}"; fi
  else
    x="${EXT:0:$((n - 20))}"
    if [ "${side}" = 3 ]; then  # EXT reads outward: reversed in front of a 3' guide
      local r="" i
      for (( i = ${#x} - 1; i >= 0; i-- )); do r="${r}${x:i:1}"; done
      echo "${r}${g}"
    else echo "${g}${x}"; fi
  fi
}
guides_for() {  # <pam name> -> guide file of that PAM's guide length, padded to its width
  local name=$1 line w side n
  line="$(cat "${PAMS}/${name}.txt")"
  w="${line##* }"
  side=3; [ "${w#-}" != "${w}" ] && side=5 && w="${w#-}"
  n="$(glen "${name}")"
  local pad; pad="$(printf 'N%.0s' $(seq 1 "${w}"))"
  local f="${WORK}/guides_${name}.txt" g1 g2
  g1="$(fit "${GUIDE}" "${side}" "${n}")"; g2="$(fit "${GUIDE2}" "${side}" "${n}")"
  if [ "${side}" = 3 ]; then printf '%s%s\n%s%s\n' "${g1}" "${pad}" "${g2}" "${pad}" > "${f}"
  else printf '%s%s\n%s%s\n' "${pad}" "${g1}" "${pad}" "${g2}" > "${f}"; fi
  echo "${f}"
}

# --- the matrix
# index: <side> <index PAM> <max bulges> | <search PAM>:<prune 0|1> ...
# (NNGRRT21, NGG21, TTTV23: guides longer than the index's 20 nt; NGG19, TTTV19:
# shorter)
# PAIRS: <index>.<search PAM>.p1=<index>.<search PAM>.p0 -- a --prune-pam search and a
# dedicated index of the same guide length and -bMax, which must give the same rows
# for every search both ran.
if [ -n "${QUICK}" ]; then
  INDEXES="
3 NGG 0 | NGG:0
3 NGG 2 | NGG:0 NGG21:0 NGG19:0
3 NNGRRT 2 | NNGRRT21:0
3 N8 2 | NGG:1 NG:1 NNGRRT:1 NNGRRT21:1 N8:0
5 TTTV 2 | TTTV:0
5 TTTV 3 | TTTV23:0
5 NNNN5 2 | TTTV:1 NNNN5:0
5 NNNN5 3 | TTTV23:1"
  SEARCHES="0,0,0 4,0,0 1,1,0 1,0,1 1,2,2"
  PAIRS="b3N82.NGG.p1=b3NGG2.NGG.p0 b3N82.NNGRRT21.p1=b3NNGRRT2.NNGRRT21.p0
b5NNNN52.TTTV.p1=b5TTTV2.TTTV.p0 b5NNNN53.TTTV23.p1=b5TTTV3.TTTV23.p0"
else
  INDEXES="
3 NGG 0 | NGG:0
3 NGG 1 | NGG:0
3 NGG 2 | NGG:0 NGG21:0 NGG19:0
3 NGG 3 | NGG:0
3 NRG 2 | NRG:0
3 NNGRRT 2 | NNGRRT:0 NNGRRT21:0
3 NNGRRT21 2 | NNGRRT21:0
3 NNN 2 | NGG:1 NRG:1 NG:1 NNN:0
3 N8 2 | NGG:1 NRG:1 NG:1 NNGRRT:1 NNGRRT21:1 NGG21:1 NGG19:1 N8:0
3 N8 1 | NGG:1 N8:0
5 TTTV 0 | TTTV:0
5 TTTV 2 | TTTV:0 TTTV19:0
5 TTTV 3 | TTTV:0 TTTV23:0
5 NNNN5 2 | TTTV:1 TTTV19:1 NNNN5:0
5 NNNN5 3 | TTTV23:1"
  SEARCHES="0,0,0 4,0,0 1,1,0 1,0,1 2,1,1 0,2,0 0,0,2 1,2,2 0,3,0 1,3,1"
  PAIRS="b3N82.NGG.p1=b3NGG2.NGG.p0 b3N82.NRG.p1=b3NRG2.NRG.p0 b3NNN2.NGG.p1=b3NGG2.NGG.p0
b3NNN2.NRG.p1=b3NRG2.NRG.p0 b3N82.NNGRRT.p1=b3NNGRRT2.NNGRRT.p0
b3N82.NNGRRT21.p1=b3NNGRRT2.NNGRRT21.p0 b3N82.NGG21.p1=b3NGG2.NGG21.p0
b3N82.NGG19.p1=b3NGG2.NGG19.p0 b5NNNN52.TTTV.p1=b5TTTV2.TTTV.p0
b5NNNN52.TTTV19.p1=b5TTTV2.TTTV19.p0 b5NNNN53.TTTV23.p1=b5TTTV3.TTTV23.p0"
fi

n_cmp=0; n_bf=0; n_rows=0; n_sub=0; n_pair=0; n_chk=0
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
  # index-info accepts every .bin of the build (exit 0, nothing on stderr: searchTST
  # would not refuse one), and a format-2 header holds the build's parameters
  ( cd "${WORK}" && "${PYTHON}" "${CZ}" index-info "${idx}" > "${name}.info.json" 2> "${name}.info.err" < /dev/null ) \
    || fail "index-info ${name}: $(head -3 "${WORK}/${name}.info.err")"
  [ ! -s "${WORK}/${name}.info.err" ] || fail "index-info ${name} wrote to stderr"
  "${PYTHON}" - "${WORK}/${name}.info.json" "${PAMS}/${ipam}.txt" "${mb}" <<'PY' || fail "index-info ${name}: a header is not the build's"
import json, sys
info = json.load(open(sys.argv[1]))
motif, limit = open(sys.argv[2]).read().split()
limit, mb = int(limit), int(sys.argv[3])
pam = motif[:-limit] if limit < 0 else motif[-limit:]
want = {"format": 2, "pam_at_start": limit < 0, "pamless": set(pam) == {"N"},
        "guide_len": len(motif) - abs(limit), "pam_len": abs(limit), "max_bulges": mb,
        "pam": pam, "complete": True}
assert info["files"], info
for name, f in info["files"].items():
    assert "error" not in f and f["num_leaves"] > 0, (name, f)
    if f["format"] == 2:
        assert {k: f[k] for k in want} == want, (name, f, want)
    else:  # format 1: an engine of 2.9.1 or earlier (the KEEP_GOING measurement)
        assert f["format"] == 1 and f["guide_len"] == want["guide_len"], (name, f)
PY
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
      # every row is its contig's own bases at its position (independent of the oracle)
      if ! "${PYTHON}" "${SCRIPT_DIR}/boundary/check_rows.py" "${genome}" "${WORK}/${tag}.targets.txt" \
           > "${WORK}/${tag}.chk"; then
        [ -n "${KEEP_GOING:-}" ] || head -5 "${WORK}/${tag}.chk"
        differs "${tag}: $(tail -1 "${WORK}/${tag}.chk")"
      fi
      n_chk=$((n_chk + $(tail -1 "${WORK}/${tag}.chk" | awk '{print $1}')))
      # 0 bulges: the brute-force scan must agree too (it needs real PAM bases, so a
      # PAM-less index searched without --prune-pam, whose flank is unconstrained, is
      # left out). A search guide shorter than the index's finds a subset of the
      # brute force: the sites whose index-length guide region fits; the brute force
      # then equals the oracle's dedicated index for the shorter guide.
      if [ "${bd}" = 0 ] && [ "${br}" = 0 ] && { [ "${prune}" = 1 ] || ! pamless "${ipam}"; }; then
        ( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${genome}/" "${PAMS}/${spam}.txt" "${gfile}" "${tag}.bf" \
            -mm "${mm}" -th 1 -r > "${tag}.bf.log" 2> "${tag}.bf.err" ) || fail "${tag}: brute-force search failed"
        [ ! -s "${WORK}/${tag}.bf.err" ] || fail "${tag}: brute-force search wrote to stderr"
        if [ "$(glen "${spam}")" -ge "$(glen "${ipam}")" ]; then
          if ! diff <(body "${WORK}/${tag}.targets.txt") <(body "${WORK}/${tag}.bf.targets.txt") > "${WORK}/${tag}.bf.diff"; then
            [ -n "${KEEP_GOING:-}" ] || { echo "  indexed (<) vs brute force (>):"; head -12 "${WORK}/${tag}.bf.diff"; }
            differs "${tag}: indexed 0-bulge search != searchBruteForce ($(grep -c '^<' "${WORK}/${tag}.bf.diff" || true) only indexed, $(grep -c '^>' "${WORK}/${tag}.bf.diff" || true) only brute force)"
          fi
        else
          only_idx="$(LC_ALL=C comm -23 <(body "${WORK}/${tag}.targets.txt") <(body "${WORK}/${tag}.bf.targets.txt") | wc -l | tr -d ' ')"
          only_bf="$(LC_ALL=C comm -13 <(body "${WORK}/${tag}.targets.txt") <(body "${WORK}/${tag}.bf.targets.txt") | wc -l | tr -d ' ')"
          [ "${only_idx}" = 0 ] || differs "${tag}: ${only_idx} rows of the indexed search are not in searchBruteForce's"
          [ "${only_bf}" -gt 0 ] || differs "${tag}: a guide shorter than the index's found every brute-force row (no cut site lost?)"
          "${PYTHON}" "${SCRIPT_DIR}/boundary/oracle.py" --genome "${genome}" --guides "${gfile}" \
            --pam "${PAMS}/${spam}.txt" --index-pam "${PAMS}/${spam}.txt" --mb 0 \
            --mm "${mm}" --bdna 0 --brna 0 > "${WORK}/${tag}.oracle_ded"
          diff -q <(body "${WORK}/${tag}.bf.targets.txt") "${WORK}/${tag}.oracle_ded" > /dev/null \
            || differs "${tag}: searchBruteForce != the oracle's dedicated index for the shorter guide"
          n_sub=$((n_sub + 1))
        fi
        n_bf=$((n_bf + 1))
      fi
    done
  done
  echo "    ok: every search of ${ipam} -bMax ${mb} == oracle (${n_cmp} comparisons so far)"
done <<< "${INDEXES}"

# --prune-pam off a PAM-less index == a dedicated index of the same guide length and -bMax
for pair in ${PAIRS}; do
  a="${pair%%=*}"; b="${pair##*=}"
  for s in ${SEARCHES}; do
    IFS=, read -r mm bd br <<< "${s}"
    fa="${WORK}/${a}.mm${mm}.b${bd}${br}.targets.txt"; fb="${WORK}/${b}.mm${mm}.b${bd}${br}.targets.txt"
    [ -f "${fa}" ] && [ -f "${fb}" ] || continue
    diff -q <(body "${fa}") <(body "${fb}") > /dev/null \
      || differs "${a} != ${b} at -mm ${mm} -bDNA ${bd} -bRNA ${br}: --prune-pam and the dedicated index differ"
    n_pair=$((n_pair + 1))
  done
done
[ "${n_pair}" -gt 0 ] || fail "no --prune-pam/dedicated pair was compared"

[ "${n_bad}" = 0 ] || fail "${n_bad} searches differ from their reference"
echo ""
echo "boundary test: PASSED (${n_cmp} indexed searches == oracle, ${n_rows} rows, ${n_chk} rows checked against their contig; ${n_bf} 0-bulge searches checked against searchBruteForce, ${n_sub} of them with a shorter guide; ${n_pair} --prune-pam searches == the dedicated index)"
