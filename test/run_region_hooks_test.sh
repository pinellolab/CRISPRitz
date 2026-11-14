#!/usr/bin/env bash
#
# Region hooks test: crispritz.py search --skip-bed, --keep-bed and --leaf-dump
# (docs/INDEX_FORMAT.md section 13).
#
# A skip BED must drop exactly the lines whose core -- the search guide and PAM,
# bulge-free, read from the line itself -- meets an interval, and a keep BED keep
# exactly those; the leaf dump must hold the stored letters of every leaf the skip
# drops, as the FASTA has them, with '-' where the index stores a placeholder. Checked
# by test/region_hooks/check_hooks.py, which computes everything from the output lines
# and the FASTA, on synthetic genomes (test/region_hooks/make_genome.py) whose BEDs put
# single bases on, and just beside, each planted site's core: on its pad, the edges of
# its guide region and PAM, and the distal flank of an 8-wide window. The searches:
#   * a dedicated NGG index: 20-nt guides at 0 and 2+2 bulges, a 21-nt guide (the pad
#     route: the core takes a pad base) and a 19-nt one;
#   * a PAM-less K=8 index under --prune-pam: NGG and NRG (the 5 distal flank bases are
#     no part of the core), 21-nt NNGRRT (both), and an 8-N PAM without --prune-pam;
#   * a 5' PAM: a dedicated TTTV index at 0 and 2+2 bulges and with 23-nt guides, and a
#     PAM-less 5' K=4 index under --prune-pam;
#   * fake-indel genomes (windows separated by one N, LEFT/RIGHT 25/27 and 30/30) with
#     their dense windows as the BED: NGG, and NNGRRT with 21-, 22- and 23-nt guides,
#     where a core reaches across the separator into the window before (section 13.5);
#     and the 21-nt guide again with the BED whose intervals hold the separator after
#     each dense window (the haplotype-scan design's form), whose core it then meets;
#   * indexes built by CRISPRitz 2.9.1 (format 1: test/fixtures/format1), searched by
#     searchTST with its 13th-15th arguments.
# For each: the search without a BED (U), with the skip BED and the leaf dump (S), with
# the keep BED (K), with the skip BED and CRISPRITZ_SKIP_PRUNE=0 (S0): S and K partition
# U by THE core, in U's order; S0 = S; SKIPPED/KEPT LEAVES add up to the index's leaves;
# the dump is complete and letter-exact. Where THE core and the prototype's leaf rule
# disagree (the K=8 index and the pad route), lines of both kinds must occur. With 4
# threads the skip gives the same lines and the same dump. Then the rule of section
# 13.5 for a scan of the dense fake-indel windows alone (test/region_hooks/hr27_bound.py,
# the dense windows' intervals with and without the separators beside them): on the pad
# route site by site (PAM-right and PAM-left sites 0-5 bases from the window's edge, 25/27
# and 30/30 windows, 21- to 26-nt guides), and where a PAM runs across a separator (a
# PAM-less index searched without --prune-pam; windows of 28, 29 and 30). Then the
# refusals (a malformed BED, an unknown contig, a skip with a keep BED, a dump without a
# skip BED, a BED on a brute-force search), an empty BED, a BED with comments and CRLF, the
# environment variables (crispritz.py does not pass them on), a searchTST without the
# hooks (crispritz.py stops), and the capability (crispritz.py version, SEARCH_HOOKS).
#
# Configuration via environment variables:
#   CXX, CXXFLAGS      : compiler and extra flags (e.g. ASan)
#   BUILDTST, SEARCHTST: prebuilt binaries (skip compiling)
#   SEARCHTST_REF      : a searchTST without the hooks (e.g. 4579026's): every search
#                        without a BED must then give its output byte for byte
#   QUICK=1            : one mismatch fewer in the searches with 3 or more bulges (CI's
#                        ASan job); default is the full matrix
#   KEEP_WORKDIR=1     : keep the work directory
#   PYTHON             : python3 to run crispritz.py and the checks (default python3)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC="${REPO_ROOT}/sourceCode"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"
PYTHON="${PYTHON:-python3}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
unset CRISPRITZ_SKIP_BED CRISPRITZ_KEEP_BED CRISPRITZ_LEAF_DUMP CRISPRITZ_SKIP_PRUNE

WORK="$(mktemp -d)"
cleanup() { [ -n "${KEEP_WORKDIR:-}" ] || rm -rf "${WORK}"; }
trap cleanup EXIT
[ -n "${KEEP_WORKDIR:-}" ] && echo "work dir: ${WORK}"
fail() { echo "FAIL: $*" >&2; exit 1; }

# --- binaries, laid out as an installed package: <root>/bin/crispritz.py and
#     <root>/opt/crispritz/<binary> (crispritz.py finds them at bin/../opt/crispritz)
ROOT="${WORK}/root"
OPT="${ROOT}/opt/crispritz"
mkdir -p "${ROOT}/bin" "${OPT}"
cp "${REPO_ROOT}/crispritz.py" "${ROOT}/bin/crispritz.py"
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
CZ="${ROOT}/bin/crispritz.py"
ST="${OPT}/searchTST"

# --- genomes, guides, BEDs, PAM files
D="${WORK}/data"
"${PYTHON}" "${SCRIPT_DIR}/region_hooks/make_genome.py" "${D}" > /dev/null
P="${WORK}/pam"; mkdir -p "${P}"
n() { printf 'N%.0s' $(seq 1 "$1"); }
printf '%sGG 3\n'   "$(n 21)" > "${P}/ngg20.txt"
printf '%sGG 3\n'   "$(n 22)" > "${P}/ngg21.txt"
printf '%sGG 3\n'   "$(n 20)" > "${P}/ngg19.txt"
printf '%sRG 3\n'   "$(n 21)" > "${P}/nrg20.txt"
printf '%sGRRT 6\n' "$(n 22)" > "${P}/nngrrt20.txt"
printf '%sGRRT 6\n' "$(n 23)" > "${P}/nngrrt21.txt"
printf '%sGRRT 6\n' "$(n 24)" > "${P}/nngrrt22.txt"
printf '%sGRRT 6\n' "$(n 25)" > "${P}/nngrrt23.txt"
printf '%s 8\n'     "$(n 28)" > "${P}/n8.txt"
printf 'TTTV%s -4\n' "$(n 20)" > "${P}/tttv20.txt"
printf 'TTTV%s -4\n' "$(n 23)" > "${P}/tttv23.txt"
printf '%s -4\n'    "$(n 24)" > "${P}/nnnn5.txt"

# cutbed <index-info json> <BED> <out>: the BED cut to the contigs the index holds
# (index-info's contigs: a contig with no site has no tree, and a BED naming it is
# refused, below), as a caller does
cutbed() {
  "${PYTHON}" - "$1" "$2" "$3" <<'PY'
import json, sys
contigs = set(json.load(open(sys.argv[1]))["contigs"])
with open(sys.argv[2]) as fi, open(sys.argv[3], "w") as fo:
    for line in fi:
        if line.startswith("#") or line.split("\t")[0] in contigs:
            fo.write(line)
PY
}

# --- indexes: name genome pam -bMax
LIB="${WORK}/lib"; mkdir -p "${LIB}"
declare -A IDX LEAVES
while read -r name genome ipam mb; do
  ( cd "${LIB}" && "${PYTHON}" "${CZ}" index-genome "${name}" "${D}/${genome}/" "${P}/${ipam}.txt" \
      -bMax "${mb}" -th 1 > "${WORK}/${name}.index.log" 2>&1 ) \
    || fail "index-genome ${name}: $(tail -3 "${WORK}/${name}.index.log")"
  IDX[${name}]="$(ls -d "${LIB}/genome_library/"*"_${mb}_${name}")"
  ( cd "${WORK}" && "${PYTHON}" "${CZ}" index-info "${IDX[${name}]}" > "${name}.info.json" ) \
    || fail "index-info ${name}"
  LEAVES[${name}]="$("${PYTHON}" -c "import json,sys; print(sum(f['num_leaves'] for f in json.load(open(sys.argv[1]))['files'].values()))" "${WORK}/${name}.info.json")"
  cutbed "${WORK}/${name}.info.json" "${D}/${genome}.bed" "${WORK}/${name}.bed"  # the genome's BED
  echo "--- index ${name}: ${ipam} -bMax ${mb} on ${genome}, $(ls "${IDX[${name}]}" | wc -l | tr -d ' ') .bin files, ${LEAVES[${name}]} leaves"
done <<'EOF'
i3ngg g3 ngg20 2
i3n8 g3 n8 2
i5tttv g5 tttv20 2
i5nnnn g5 nnnn5 2
f25ngg fake25 ngg20 2
f25grrt fake25 nngrrt20 3
f30ngg fake30 ngg20 2
f30grrt fake30 nngrrt20 3
EOF

# search <tag> <index> <search pam> <guides> <prune 0|1> <mm> <bDNA> <bRNA> [hook flags...]
search() {
  local tag=$1 idx=$2 spam=$3 guides=$4 prune=$5 mm=$6 bd=$7 br=$8; shift 8
  local flag=""; [ "${prune}" = 1 ] && flag="--prune-pam"
  ( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${IDX[${idx}]}/" "${P}/${spam}.txt" "${D}/${guides}" "${tag}" \
      -mm "${mm}" -bDNA "${bd}" -bRNA "${br}" -th "${TH:-1}" -r ${flag} "$@" > "${tag}.out" 2> "${tag}.err" ) \
    || fail "${tag}: search exited non-zero: $(tail -3 "${WORK}/${tag}.err")"
  [ ! -s "${WORK}/${tag}.err" ] || fail "${tag}: search wrote to stderr: $(head -3 "${WORK}/${tag}.err")"
}

n_cases=0; n_lines=0; n_ref=0
declare -A PARAMS
# tag | index | search pam | guides | prune | mm bd br | genome | bed | G K mb Gs L | checker args
while IFS='|' read -r tag idx spam guides prune budget genome bed geo extra; do
  [ -z "${tag// /}" ] && continue
  tag=${tag// /}; idx=${idx// /}; spam=${spam// /}; guides=${guides// /}; prune=${prune// /}
  genome=${genome// /}; bed=${bed// /}; extra=${extra//DATA/${D}}
  read -r mm bd br <<< "${budget}"; read -r G K mb gs L <<< "${geo}"
  [ -n "${QUICK:-}" ] && [ $((bd + br)) -ge 3 ] && [ "${mm}" -gt 1 ] && mm=$((mm - 1))
  PARAMS[${tag}]="${idx} ${spam} ${guides} ${prune} ${mm} ${bd} ${br}"
  B="${WORK}/${tag}.bed"; cutbed "${WORK}/${idx}.info.json" "${D}/${bed}" "${B}"  # the case's BED
  search "${tag}.U" "${idx}" "${spam}" "${guides}" "${prune}" "${mm}" "${bd}" "${br}"
  if [ -n "${SEARCHTST_REF:-}" ]; then  # without a BED: byte for byte the engine without the hooks
    ( cd "${WORK}" && "${SEARCHTST_REF}" "${IDX[${idx}]}/" "${D}/${guides}" "${mm}" "${bd}" "${br}" \
        "${P}/${spam}.txt" "${tag}.REF" r 1 "${mb}" -1 "${prune}" > "${tag}.REF.out" 2> "${tag}.REF.err" ) \
      || fail "${tag}: the reference searchTST failed"
    cmp -s "${WORK}/${tag}.U.targets.txt" "${WORK}/${tag}.REF.targets.txt" \
      || fail "${tag}: without a BED the output differs from SEARCHTST_REF's"
    n_ref=$((n_ref + 1))
  fi
  search "${tag}.S" "${idx}" "${spam}" "${guides}" "${prune}" "${mm}" "${bd}" "${br}" \
    --skip-bed "${B}" --leaf-dump "${WORK}/${tag}.dump.tsv"
  search "${tag}.K" "${idx}" "${spam}" "${guides}" "${prune}" "${mm}" "${bd}" "${br}" --keep-bed "${B}"
  CRISPRITZ_SKIP_PRUNE=0 search "${tag}.S0" "${idx}" "${spam}" "${guides}" "${prune}" "${mm}" "${bd}" "${br}" \
    --skip-bed "${B}"
  pflag=""; [ "${prune}" = 1 ] && pflag="--prune"
  side=3; case "${spam}" in tttv*) side=5;; esac
  "${PYTHON}" "${SCRIPT_DIR}/region_hooks/check_hooks.py" --label "${tag}" --ordered \
    --unskipped "${WORK}/${tag}.U.targets.txt" --skip "${WORK}/${tag}.S.targets.txt" \
    --skip-out "${WORK}/${tag}.S.out" --keep "${WORK}/${tag}.K.targets.txt" --keep-out "${WORK}/${tag}.K.out" \
    --skip-noprune "${WORK}/${tag}.S0.targets.txt" --bed "${B}" --pam-side "${side}" \
    --G "${G}" --K "${K}" --mb "${mb}" --gs "${gs}" --L "${L}" ${pflag} \
    --dump "${WORK}/${tag}.dump.tsv" --genome "${D}/${genome}" --leaves "${LEAVES[${idx}]}" ${extra} \
    || fail "${tag}: the region hooks are not exact (above)"
  n_cases=$((n_cases + 1)); n_lines=$((n_lines + $(grep -vc '^#' "${WORK}/${tag}.U.targets.txt")))
done <<'EOF'
d3_ngg_b00   | i3ngg  | ngg20    | guides_ngg20.txt         | 0 | 3 0 0 | g3     | g3.bed     | 20 3 2 20 3 |
d3_ngg_b22   | i3ngg  | ngg20    | guides_ngg20.txt         | 0 | 2 2 2 | g3     | g3.bed     | 20 3 2 20 3 |
d3_ngg21     | i3ngg  | ngg21    | guides_ngg21.txt         | 0 | 3 1 2 | g3     | g3.bed     | 20 3 2 21 3 | --min-old-dup 1
d3_ngg19     | i3ngg  | ngg19    | guides_ngg19.txt         | 0 | 3 2 1 | g3     | g3.bed     | 20 3 2 19 3 |
k8_ngg       | i3n8   | ngg20    | guides_ngg20.txt         | 1 | 2 2 2 | g3     | g3.bed     | 20 8 2 20 3 | --min-old-lost 1
k8_nrg       | i3n8   | nrg20    | guides_ngg20.txt         | 1 | 3 1 1 | g3     | g3.bed     | 20 8 2 20 3 | --min-old-lost 1
k8_nngrrt21  | i3n8   | nngrrt21 | guides_nngrrt21.txt      | 1 | 3 1 2 | g3     | g3.bed     | 20 8 2 21 6 | --min-old-lost 1 --min-old-dup 1
k8_n8        | i3n8   | n8       | guides_n8.txt            | 0 | 2 1 1 | g3     | g3.bed     | 20 8 2 20 8 |
d5_tttv_b00  | i5tttv | tttv20   | guides_tttv20.txt        | 0 | 3 0 0 | g5     | g5.bed     | 20 4 2 20 4 |
d5_tttv_b22  | i5tttv | tttv20   | guides_tttv20.txt        | 0 | 2 2 2 | g5     | g5.bed     | 20 4 2 20 4 |
d5_tttv23    | i5tttv | tttv23   | guides_tttv23.txt        | 0 | 3 1 2 | g5     | g5.bed     | 20 4 2 23 4 | --min-old-dup 1
k4_5_tttv    | i5nnnn | tttv20   | guides_tttv20.txt        | 1 | 2 2 2 | g5     | g5.bed     | 20 4 2 20 4 |
f25_ngg      | f25ngg | ngg20    | guides_fake_ngg20.txt    | 0 | 3 2 2 | fake25 | fake25.bed | 20 3 2 20 3 |
f25_grrt21   | f25grrt| nngrrt21 | guides_fake_nngrrt21.txt | 0 | 2 0 3 | fake25 | fake25.bed | 20 6 3 21 6 | --windows DATA/fake25.windows.tsv --hr27 none --min-dropped 0 --min-kept 0
f25_grrt22   | f25grrt| nngrrt22 | guides_fake_nngrrt22.txt | 0 | 2 0 3 | fake25 | fake25.bed | 20 6 3 22 6 | --windows DATA/fake25.windows.tsv --hr27 any --min-kept 0
f25_grrt23   | f25grrt| nngrrt23 | guides_fake_nngrrt23.txt | 0 | 2 0 3 | fake25 | fake25.bed | 20 6 3 23 6 | --windows DATA/fake25.windows.tsv --hr27 any --min-kept 0
f25s_grrt21  | f25grrt| nngrrt21 | guides_fake_nngrrt21.txt | 0 | 2 0 3 | fake25 | fake25.sep.bed | 20 6 3 21 6 | --windows DATA/fake25.windows.tsv --hr27 any --min-kept 0
f30_ngg      | f30ngg | ngg20    | guides_fake_ngg20.txt    | 0 | 3 2 2 | fake30 | fake30.bed | 20 3 2 20 3 |
f30_grrt22   | f30grrt| nngrrt22 | guides_fake_nngrrt22.txt | 0 | 2 0 3 | fake30 | fake30.bed | 20 6 3 22 6 | --windows DATA/fake30.windows.tsv --hr27 none --min-dropped 0 --min-kept 0
f30_grrt23   | f30grrt| nngrrt23 | guides_fake_nngrrt23.txt | 0 | 2 0 3 | fake30 | fake30.bed | 20 6 3 23 6 | --windows DATA/fake30.windows.tsv --hr27 none --min-dropped 0 --min-kept 0
f30s_grrt21  | f30grrt| nngrrt21 | guides_fake_nngrrt21.txt | 0 | 2 0 3 | fake30 | fake30.sep.bed | 20 6 3 21 6 | --windows DATA/fake30.windows.tsv --hr27 none --min-kept 0
EOF
echo "--- ${n_cases} searches: skip + keep = the search without a BED, by THE core; dumps exact (${n_lines} lines)"
[ -n "${SEARCHTST_REF:-}" ] && echo "--- ${n_ref} searches without a BED byte-identical to SEARCHTST_REF"

# --- indexes built by CRISPRitz 2.9.1 (format 1: the committed fixtures of
#     test/run_format1_test.sh; max_bulges from the command line; an index-0 artefact on
#     two contigs, which is no site and is never dropped, dumped or counted), searched by
#     searchTST directly with its 13th-15th arguments
F1="${SCRIPT_DIR}/fixtures/format1"
dsearch() {  # <tag> <index> <guides> <pam> <prune> [skip keep dump]: -mm 2 -bDNA 1 -bRNA 1, max_bulges 2
  local tag=$1 idx=$2 guides=$3 pam=$4 prune=$5; shift 5
  ( cd "${WORK}" && "${ST}" "${F1}/${idx}/" "${F1}/${guides}" 2 1 1 "${F1}/${pam}" "${tag}" r 1 2 -1 "${prune}" "$@" \
      > "${tag}.out" 2> "${tag}.err" ) || fail "${tag}: searchTST exited non-zero: $(tail -3 "${WORK}/${tag}.err")"
  [ ! -s "${WORK}/${tag}.err" ] || fail "${tag}: searchTST wrote to stderr: $(head -3 "${WORK}/${tag}.err")"
}
# (TTTV: every line comes from one site, which is dropped)
for c in "f1_ngg idx_NGG_2 guide3.txt pam_NGG.txt 0 3 20 3 3 24 1" "f1_k8_ngg idx_N8_2 guide3.txt pam_NGG.txt 1 8 20 3 3 24 1" \
         "f1_tttv idx_TTTV_2 guide5.txt pam_TTTV.txt 0 4 20 4 5 14 0"; do
  read -r tag idx guides pam prune K G L side every minkept <<< "${c}"
  ( cd "${WORK}" && "${PYTHON}" "${CZ}" index-info "${F1}/${idx}" > "${tag}.info.json" ) || fail "index-info ${idx}"
  "${PYTHON}" - "${F1}/genome" "${WORK}/${tag}.info.json" "${WORK}/${tag}.bed" "${every}" <<'PY'
import json, os, random, sys
rng = random.Random(11)
every = int(sys.argv[4])  # a BED base about every that many bases
contigs = set(json.load(open(sys.argv[2]))["contigs"])
with open(sys.argv[3], "w") as fo:
    for name in sorted(os.listdir(sys.argv[1])):
        seq = "".join(l.strip() for l in open(os.path.join(sys.argv[1], name)) if not l.startswith(">"))
        for p in range(len(seq)):
            if rng.random() < 1 / every and name[:-3] in contigs:
                fo.write("%s\t%d\t%d\n" % (name[:-3], p, p + 1))
PY
  dsearch "${tag}.U" "${idx}" "${guides}" "${pam}" "${prune}"
  # and one base: the first line's core start (so its site, at least, is dropped)
  "${PYTHON}" - "${WORK}/${tag}.U.targets.txt" "${side}" >> "${WORK}/${tag}.bed" <<'PY'
import sys
f = [l for l in open(sys.argv[1]) if not l.startswith("#")][0].split("\t")
d, r, p = f[1].count("-"), f[2].count("-"), int(f[4])
cs = p + d - r if (f[6] == "+") == (sys.argv[2] == "3") else p
print("%s\t%d\t%d" % (f[3], cs, cs + 1))
PY
  dsearch "${tag}.S" "${idx}" "${guides}" "${pam}" "${prune}" "${WORK}/${tag}.bed" "" "${WORK}/${tag}.dump.tsv"
  dsearch "${tag}.K" "${idx}" "${guides}" "${pam}" "${prune}" "" "${WORK}/${tag}.bed" ""
  CRISPRITZ_SKIP_PRUNE=0 dsearch "${tag}.S0" "${idx}" "${guides}" "${pam}" "${prune}" "${WORK}/${tag}.bed"
  pflag=""; [ "${prune}" = 1 ] && pflag="--prune"
  "${PYTHON}" "${SCRIPT_DIR}/region_hooks/check_hooks.py" --label "${tag}" --ordered \
    --unskipped "${WORK}/${tag}.U.targets.txt" --skip "${WORK}/${tag}.S.targets.txt" --skip-out "${WORK}/${tag}.S.out" \
    --keep "${WORK}/${tag}.K.targets.txt" --keep-out "${WORK}/${tag}.K.out" --skip-noprune "${WORK}/${tag}.S0.targets.txt" \
    --bed "${WORK}/${tag}.bed" --pam-side "${side}" --G "${G}" --K "${K}" --mb 2 --gs 20 --L "${L}" ${pflag} \
    --dump "${WORK}/${tag}.dump.tsv" --genome "${F1}/genome" --dup-ok --min-kept "${minkept}" \
    || fail "${tag}: the region hooks are not exact on a format-1 index (above)"
  n_cases=$((n_cases + 1)); n_lines=$((n_lines + $(grep -vc '^#' "${WORK}/${tag}.U.targets.txt")))
done
echo "--- format 1 (indexes built by 2.9.1): 3 searches, the same"

# --- 4 threads: the same lines (in another order across trees) and the same dump
for tag in d3_ngg_b22 k8_ngg f25_ngg; do
  read -r idx spam guides prune mm bd br <<< "${PARAMS[${tag}]}"
  TH=4 OMP_NUM_THREADS=4 search "${tag}.S4" "${idx}" "${spam}" "${guides}" "${prune}" "${mm}" "${bd}" "${br}" \
    --skip-bed "${WORK}/${idx}.bed" --leaf-dump "${WORK}/${tag}.dump4.tsv"
  cmp -s <(LC_ALL=C sort "${WORK}/${tag}.S.targets.txt") <(LC_ALL=C sort "${WORK}/${tag}.S4.targets.txt") \
    || fail "${tag}: 4 threads skip other lines than 1"
  cmp -s "${WORK}/${tag}.dump.tsv" "${WORK}/${tag}.dump4.tsv" || fail "${tag}: 4 threads dump other leaves than 1"
  [ "$(grep '^SKIPPED LEAVES' "${WORK}/${tag}.S.out")" = "$(grep '^SKIPPED LEAVES' "${WORK}/${tag}.S4.out")" ] \
    || fail "${tag}: 4 threads count other leaves than 1"
done
echo "--- 4 threads: the same lines, dump and counts as 1 (3 searches)"

# --- the rule of section 13.5 for a scan of the dense fake-indel windows alone (review
#     findings H-R2-7, HK-R1-1), site by site: in 25/27 and 30/30 windows, a PAM-right and
#     a PAM-left NNGRRT site 0 to 5 bases from the window edge its core reaches across
#     (between two dense windows), searched with 21- to 26-nt guides on a 20-nt index
#     (n = 27..32, -mm 1, 3 RNA bulges) and skipped with the dense windows in four forms:
#     without the separators, with the one after each (the haplotype-scan design's INDELS
#     BED), the one before, both. The sites whose dropped lines lie in their own window and
#     cover its ALT -- what such a scan misses -- must be the rule's, in all 48 skips. And
#     where a line's PAM runs across a separator (a PAM-less K=8 index searched with 8 N's
#     and no --prune-pam, whose flank is not tested; windows of 28, 29 and 30 bases, the
#     guide region next to the window's edge on its PAM side): the sides that lose lines
#     must be the rule's, in 12 skips
"${PYTHON}" "${SCRIPT_DIR}/region_hooks/hr27_bound.py" genome "${D}" > /dev/null
for gs in 24 25 26; do printf '%sGRRT 6\n' "$(n $((gs + 2)))" > "${P}/nngrrt${gs}.txt"; done
for g in hr25 hr30; do
  ( cd "${LIB}" && "${PYTHON}" "${CZ}" index-genome "${g}" "${D}/${g}/" "${P}/nngrrt20.txt" -bMax 3 -th 1 \
      > "${WORK}/${g}.index.log" 2>&1 ) || fail "index-genome ${g}: $(tail -3 "${WORK}/${g}.index.log")"
  IH="$(ls -d "${LIB}/genome_library/"*"_3_${g}")"
  for gs in 21 22 23 24 25 26; do
    for form in U win trail lead both; do
      hooks=(); [ "${form}" = U ] || hooks=("${D}/${g}.${form}.bed" "" "")
      ( cd "${WORK}" && "${ST}" "${IH}/" "${D}/guides_hr${gs}.txt" 1 0 3 "${P}/nngrrt${gs}.txt" "${g}.${gs}.${form}" \
          r 1 3 -1 0 ${hooks[@]+"${hooks[@]}"} > "${g}.${gs}.${form}.out" 2> "${g}.${gs}.${form}.err" ) \
        || fail "${g}.${gs}.${form}: searchTST exited non-zero: $(tail -3 "${WORK}/${g}.${gs}.${form}.err")"
      [ ! -s "${WORK}/${g}.${gs}.${form}.err" ] || fail "${g}.${gs}.${form}: searchTST wrote to stderr"
    done
  done
  "${PYTHON}" "${SCRIPT_DIR}/region_hooks/hr27_bound.py" check "${D}" "${g}" "${WORK}" 20 6 3 3 21 22 23 24 25 26 \
    || fail "${g}: the sites a scan of the dense windows alone misses are not the rule's (above)"
done
for g in st28 st29 st30; do
  ( cd "${LIB}" && "${PYTHON}" "${CZ}" index-genome "${g}" "${D}/${g}/" "${P}/n8.txt" -bMax 3 -th 1 \
      > "${WORK}/${g}.index.log" 2>&1 ) || fail "index-genome ${g}: $(tail -3 "${WORK}/${g}.index.log")"
  IH="$(ls -d "${LIB}/genome_library/"*"_3_${g}")"
  for form in U win trail lead both; do
    hooks=(); [ "${form}" = U ] || hooks=("${D}/${g}.${form}.bed" "" "")
    ( cd "${WORK}" && "${ST}" "${IH}/" "${D}/guides_st.txt" 0 3 0 "${P}/n8.txt" "${g}.${form}" \
        r 1 3 -1 0 ${hooks[@]+"${hooks[@]}"} > "${g}.${form}.out" 2> "${g}.${form}.err" ) \
      || fail "${g}.${form}: searchTST exited non-zero: $(tail -3 "${WORK}/${g}.${form}.err")"
    [ ! -s "${WORK}/${g}.${form}.err" ] || fail "${g}.${form}: searchTST wrote to stderr"
  done
  "${PYTHON}" "${SCRIPT_DIR}/region_hooks/hr27_bound.py" check-pam "${D}" "${g}" "${WORK}" 20 8 3 \
    || fail "${g}: the lines a scan of the dense windows alone misses are not the rule's (above)"
done
echo "--- section 13.5's rule for a scan of the dense fake-indel windows alone: what it misses is the rule's, in 60 skips"

# --- the refusals: exit 1, the reason on stderr, no targets file
I="${IDX[i3ngg]}"
direct() {  # <out> <skip> <keep> <dump>: searchTST with the hooks as its 13th-15th arguments
  ( cd "${WORK}" && "${ST}" "${I}/" "${D}/guides_ngg20.txt" 2 1 1 "${P}/ngg20.txt" "$1" r 1 2 -1 0 "$2" "$3" "$4" )
}
refused() {  # <label> <expected stderr text> <command...>
  local label=$1 want=$2; shift 2
  rm -f "${WORK}/ref.targets.txt"
  if "$@" > "${WORK}/ref.stdout" 2> "${WORK}/ref.err"; then fail "${label}: accepted"; fi
  grep -qF -- "${want}" "${WORK}/ref.err" || fail "${label}: stderr lacks '${want}': $(head -2 "${WORK}/ref.err")"
  [ ! -e "${WORK}/ref.targets.txt" ] || fail "${label}: a targets file was written"
}
bad() { printf "$2" > "${WORK}/bad_$1.bed"; refused "BED $1" "$3" direct ref "${WORK}/bad_$1.bed" "" ""; }
bad cols      'chrA\t10\t20\nchrA\t30\n'        "line 2: fewer than 3 tab-separated columns"
bad spaces    'chrA 10 20\n'                    "line 1: fewer than 3 tab-separated columns"
bad equal     'chrA\t10\t10\n'                  "line 1: the end 10 is not greater than the start 10"
bad reversed  'chrA\t10\t5\n'                   "the end 5 is not greater than the start 10"
bad float     'chrA\t10.5\t20\n'                "the start '10.5' is not a non-negative integer"
bad negative  'chrA\t-1\t20\n'                  "the start '-1' is not a non-negative integer"
bad space     'chrA\t 10\t20\n'                 "the start ' 10' is not a non-negative integer"
bad endtext   '# a comment\nchrA\t1\t2x\n'      "line 2: the end '2x' is not a non-negative integer"
bad huge      'chrA\t1\t1234567890123456789\n'  "the end '1234567890123456789' is not a non-negative integer"
bad nocontig  '\t1\t2\n'                        "line 1: an empty contig name"
bad unknown   'chrA\t1\t2\nchrZ\t1\t2\n'        "line 2: the contig 'chrZ' is not one of the index's"
bad fake      'fakechrA\t1\t2\n'                "the contig 'fakechrA' is not one of the index's"
refused "BED missing" "cannot be opened" direct ref "${WORK}/no_such.bed" "" ""
refused "a skip and a keep BED" "a skip BED or a keep BED, not both" direct ref "${D}/g3.bed" "${D}/g3.bed" ""
refused "CRISPRITZ_SKIP_BED with CRISPRITZ_KEEP_BED" "a skip BED or a keep BED, not both" \
  env CRISPRITZ_SKIP_BED="${WORK}/i3ngg.bed" CRISPRITZ_KEEP_BED="${WORK}/i3ngg.bed" \
  "${ST}" "${I}/" "${D}/guides_ngg20.txt" 2 1 1 "${P}/ngg20.txt" "${WORK}/ref" r 1 2 -1 0
refused "a dump with no BED" "has no skip BED" direct ref "" "" "${WORK}/d.tsv"
refused "a dump with a keep BED" "has a keep BED, not a skip BED" direct ref "" "${WORK}/i3ngg.bed" "${WORK}/d.tsv"
refused "a dump into no folder" "cannot be written" direct ref "${WORK}/i3ngg.bed" "" "${WORK}/no/such/dir/d.tsv"
# a contig of the genome that holds no site, so no tree: the TTTV index has none for chrB_0
printf 'chrB_0\t0\t5\n' > "${WORK}/b0.bed"
refused "a contig with no tree" "the contig 'chrB_0' is not one of the index's" \
  "${ST}" "${IDX[i5tttv]}/" "${D}/guides_tttv20.txt" 2 0 0 "${P}/tttv20.txt" "${WORK}/ref" r 1 2 -1 0 "${WORK}/b0.bed"
cz() { ( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${I}/" "${P}/ngg20.txt" "${D}/guides_ngg20.txt" ref -mm 2 -bDNA 1 -bRNA 1 -r "$@" ); }
refused "crispritz.py: a skip and a keep BED" "--skip-bed or --keep-bed, not both" cz --skip-bed "${WORK}/i3ngg.bed" --keep-bed "${WORK}/i3ngg.bed"
refused "crispritz.py: a dump without a skip BED" "--leaf-dump writes the leaves a --skip-bed BED drops" cz --leaf-dump "${WORK}/d.tsv"
refused "crispritz.py: no path" "--skip-bed takes one path" cz -th 1 --skip-bed
refused "crispritz.py: a flag twice" "given 2 times" cz --skip-bed "${WORK}/i3ngg.bed" --skip-bed "${WORK}/i3ngg.bed"
refused "crispritz.py: a missing BED" "no such file" cz --keep-bed "${WORK}/no_such.bed"
refused "crispritz.py: a malformed BED" "the end 10 is not greater than the start 10" cz --skip-bed "${WORK}/bad_equal.bed"
refused "crispritz.py: a BED on a brute-force search" "act on an indexed search" \
  bash -c "cd '${WORK}' && '${PYTHON}' '${CZ}' search '${D}/g3/' '${P}/ngg20.txt' '${D}/guides_ngg20.txt' ref -mm 2 -r --skip-bed '${WORK}/i3ngg.bed'"
echo "--- refusals: 12 malformed BEDs, a missing one, skip with keep (arguments, variables), a dump without a skip BED or into no folder, a contig with no tree; crispritz.py: 7"

# the later checks reuse d3_ngg_b22's search: its index, PAM, guides and budget
read -r _ _ _ _ MM BD BR <<< "${PARAMS[d3_ngg_b22]}"

# --- accepted: an empty BED (the search without one, byte for byte, and SKIPPED 0); the
#     BED with comments, empty lines, CRLF, extra columns, unsorted and overlapping
#     intervals (the clean BED's output); the dump and SKIPPED LEAVES of a search with no
#     dropped leaf
printf '# nothing\n\n' > "${WORK}/empty.bed"
search e.S i3ngg ngg20 guides_ngg20.txt 0 ${MM} ${BD} ${BR} --skip-bed "${WORK}/empty.bed" --leaf-dump "${WORK}/e.dump.tsv"
cmp -s "${WORK}/e.S.targets.txt" "${WORK}/d3_ngg_b22.U.targets.txt" || fail "an empty skip BED changed the output"
grep -qx "SKIPPED LEAVES	0	KEPT LEAVES	${LEAVES[i3ngg]}" "${WORK}/e.S.out" || fail "an empty skip BED: SKIPPED LEAVES"
grep -qx "LEAF DUMP	0" "${WORK}/e.S.out" && [ ! -s "${WORK}/e.dump.tsv" ] || fail "an empty skip BED: the dump"
search e.K i3ngg ngg20 guides_ngg20.txt 0 ${MM} ${BD} ${BR} --keep-bed "${WORK}/empty.bed"
[ "$(grep -vc '^#' "${WORK}/e.K.targets.txt")" = 0 ] || fail "an empty keep BED kept lines"
"${PYTHON}" - "${WORK}/i3ngg.bed" "${WORK}/messy.bed" <<'PY'
import random, sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1]) if not l.startswith("#")]
out = ["# a header", ""]
for c, s, e in rows:
    out.append("%s\t%s\t%s\tname\t0\t+" % (c, s, e))
    if int(e) - int(s) > 2:  # an overlapping piece of it
        out.append("%s\t%d\t%d" % (c, int(s) + 1, int(e) - 1))
random.Random(7).shuffle(out)
open(sys.argv[2], "w", newline="").write("\r\n".join(out) + "\r\n")
PY
search m.S i3ngg ngg20 guides_ngg20.txt 0 ${MM} ${BD} ${BR} --skip-bed "${WORK}/messy.bed" --leaf-dump "${WORK}/m.dump.tsv"
cmp -s "${WORK}/m.S.targets.txt" "${WORK}/d3_ngg_b22.S.targets.txt" && cmp -s "${WORK}/m.dump.tsv" "${WORK}/d3_ngg_b22.dump.tsv" \
  || fail "the BED with comments, CRLF, extra columns and overlaps gave other lines or another dump"
echo "--- accepted: an empty BED (the search without one; keep: nothing), a BED with comments, CRLF, extra columns, overlaps (the clean BED's lines)"

# --- the environment: crispritz.py does not pass CRISPRITZ_SKIP_BED, _KEEP_BED or
#     _LEAF_DUMP on (a search given no flag is the search without the hooks); searchTST
#     reads them where its argument is absent or empty
CRISPRITZ_SKIP_BED="${WORK}/i3ngg.bed" CRISPRITZ_LEAF_DUMP="${WORK}/env.dump.tsv" \
  search env.U i3ngg ngg20 guides_ngg20.txt 0 ${MM} ${BD} ${BR}
cmp -s "${WORK}/env.U.targets.txt" "${WORK}/d3_ngg_b22.U.targets.txt" && ! grep -q 'SKIPPED LEAVES' "${WORK}/env.U.out" \
  && [ ! -e "${WORK}/env.dump.tsv" ] || fail "crispritz.py passed a CRISPRITZ_ variable on to searchTST"
( cd "${WORK}" && CRISPRITZ_SKIP_BED="${WORK}/i3ngg.bed" "${ST}" "${I}/" "${D}/guides_ngg20.txt" "${MM}" "${BD}" "${BR}" "${P}/ngg20.txt" env.S r 1 2 -1 0 > env.S.out ) \
  || fail "searchTST with CRISPRITZ_SKIP_BED failed"
( cd "${WORK}" && CRISPRITZ_KEEP_BED="${WORK}/i3ngg.bed" "${ST}" "${I}/" "${D}/guides_ngg20.txt" "${MM}" "${BD}" "${BR}" "${P}/ngg20.txt" env.K r 1 2 -1 0 "" "" "" > env.K.out ) \
  || fail "searchTST with CRISPRITZ_KEEP_BED failed"
cmp -s "${WORK}/env.S.targets.txt" "${WORK}/d3_ngg_b22.S.targets.txt" && cmp -s "${WORK}/env.K.targets.txt" "${WORK}/d3_ngg_b22.K.targets.txt" \
  || fail "searchTST read CRISPRITZ_SKIP_BED / CRISPRITZ_KEEP_BED otherwise than its arguments"
echo "--- environment: crispritz.py does not pass the CRISPRITZ_ variables on; searchTST reads them as its arguments"

# --- a searchTST without the hooks (before 2.10.0 ignores the 13th-15th arguments):
#     crispritz.py search --skip-bed must stop rather than return every site
OLDROOT="${WORK}/oldroot"; mkdir -p "${OLDROOT}/bin" "${OLDROOT}/opt/crispritz"
cp "${CZ}" "${OLDROOT}/bin/"
cat > "${OLDROOT}/opt/crispritz/searchTST" <<SH
#!/bin/bash
exec "${ST}" "\${@:1:12}"
SH
chmod +x "${OLDROOT}/opt/crispritz/searchTST"
refused "crispritz.py on a searchTST without the hooks" "did not apply the region BED" \
  bash -c "cd '${WORK}' && '${PYTHON}' '${OLDROOT}/bin/crispritz.py' search '${I}/' '${P}/ngg20.txt' '${D}/guides_ngg20.txt' old -mm 2 -bDNA 1 -bRNA 1 -r --skip-bed '${WORK}/i3ngg.bed'"
echo "--- a searchTST without the hooks: crispritz.py search --skip-bed stops"

# --- the capability: crispritz.py version lists the hooks, and SEARCH_HOOKS states them
#     in the form CRISPRme reads VERSION and INDEX_FORMAT in
"${PYTHON}" "${CZ}" version > "${WORK}/version.out"
grep -qx "search hooks: skip-bed keep-bed leaf-dump" "${WORK}/version.out" || fail "crispritz.py version: $(cat "${WORK}/version.out")"
"${PYTHON}" - "${CZ}" <<'PY' || fail "SEARCH_HOOKS is not readable as CRISPRme reads crispritz.py's VERSION"
import re, sys
vals = {}
for line in open(sys.argv[1]):
    m = re.match(r"""(VERSION|INDEX_FORMAT|SEARCH_HOOKS)\s*=\s*["']?([^"'\s#]+)""", line)
    if m and m.group(1) not in vals:
        vals[m.group(1)] = m.group(2)
assert vals.get("SEARCH_HOOKS") == "skip-bed,keep-bed,leaf-dump", vals
PY
echo "--- capability: crispritz.py version lists skip-bed keep-bed leaf-dump; SEARCH_HOOKS"
echo "region-hooks test: PASSED (${n_cases} searches by THE core, ${n_lines} lines)"
