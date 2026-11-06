#!/usr/bin/env bash
#
# The four-allele symbol Q, end to end through crispritz.py (docs/INDEX_FORMAT.md
# section 4, docs/MULTIALLELIC_ENRICHMENT.md section 2.2).
#
# add-variants writes Q where a position's reference base and SNP ALTs cover A, C, G
# and T, and buildTST, searchTST and searchBruteForce read it as any base (1111); N
# stays a gap (0000) that matches nothing. Until the fifth symbol, such a position was
# written N (a record listing three SNP ALTs: every target through it lost) or, after
# the four-base rule, without its rarest ALT.
#
# four_allele/make_genome.py writes contigs and VCFs whose split and multi-ALT SNP
# records put Q in guide regions (making a reference mismatch a match), in PAMs (a
# PAM that exists only through Q), beyond a guide's far end (the pad a DNA bulge
# reads), next to N runs and on contig ends, a run of 30 Q's, ALTs on reference N
# (which stay N), and a background with Q, IUPAC codes and N runs. Then:
#
#   1. ENRICH  : add-variants with each enricher port; the two enriched genomes are
#                identical, Q exactly where the reference base and ALTs cover four
#                bases, and an N exactly where the reference has one.
#   2. INDEX   : index-genome (NGG at -bMax 0 and 2, NRG, a PAM-less K=8 flank, 5'
#                TTTV and a 5' PAM-less K=4 flank); a .bin holding Q is flagged
#                four_allele in its header (index-info), the variant-free contig's is
#                not, and buildTST says so on stdout.
#   3. SEARCH  : every search, with and without bulges and --prune-pam, equals the
#                independent oracle (test/boundary/oracle.py, Q = any base, N = gap)
#                as a set of rows; every 0-bulge one equals searchBruteForce, so every
#                target the brute-force scan finds is found by searchTST; every row is
#                its contig's own bases (check_rows.py), and no row holds an N; every
#                planted perfect site is found; --prune-pam equals the dedicated index.
#   4. GAPS    : searchBruteForce with a mismatch budget of 120 (more than the 101 an N
#                used to cost) reports no window with an N in its guide, as an index
#                never holds one.
#   5. BYTES   : a genome byte that is no base, IUPAC code, N or Q is a gap too: the PAM
#                scans of buildTST and searchBruteForce used to skip it, reading every
#                later PAM one base off, so a perfect site after it was lost.
#
# Any stage writing to stderr fails the test (stderr is fatal to a CRISPRme stage).
#
# Configuration: CXX, CXXFLAGS (e.g. ASan), BUILDTST, SEARCHTST, SEARCHBF, ENRICHER
# (prebuilt binaries), PYTHON, KEEP_WORKDIR.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC="${REPO_ROOT}/sourceCode"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:-}"
PYTHON="${PYTHON:-python3}"
PYTHON="$(command -v "${PYTHON}")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

WORK="$(mktemp -d)"
cleanup() { [ -n "${KEEP_WORKDIR:-}" ] || rm -rf "${WORK}"; }
trap cleanup EXIT
[ -n "${KEEP_WORKDIR:-}" ] && echo "work dir: ${WORK}"
fail() { echo "FAIL: $*" >&2; exit 1; }
quiet() {  # <label> <log> <err>: the stage wrote nothing to stderr
  [ ! -s "$3" ] || fail "$1 wrote to stderr: $(head -3 "$3")"
}

# --- binaries, laid out as an installed package: <root>/bin/crispritz.py and
#     <root>/opt/crispritz/<binary> (crispritz.py finds them at bin/../opt/crispritz)
ROOT="${WORK}/root"
OPT="${ROOT}/opt/crispritz"
mkdir -p "${ROOT}/bin" "${OPT}/Python_Scripts/Enrichment"
cp "${REPO_ROOT}/crispritz.py" "${ROOT}/bin/crispritz.py"
cp "${SRC}/Python_Scripts/Enrichment/enricher.py" "${OPT}/Python_Scripts/Enrichment/"
build() {  # <name> <env override> <compile command...>
  local name=$1 pre=$2; shift 2
  if [ -n "${pre}" ]; then cp "${pre}" "${OPT}/${name}"; else echo "--- compiling ${name} ---"; "$@"; fi
}
build buildTST "${BUILDTST:-}" "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} \
  "${SRC}/CRISPR-Cas-Tree/mainParallel.cpp" -o "${OPT}/buildTST"
build searchTST "${SEARCHTST:-}" "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} \
  "${SRC}/CRISPR-Cas-Tree/searchOnTST.cpp" "${SRC}/CRISPR-Cas-Tree/detailedOutput.cpp" \
  "${SRC}/CRISPR-Cas-Tree/convert.cpp" -I "${SRC}/CRISPR-Cas-Tree/include" -o "${OPT}/searchTST"
build searchBruteForce "${SEARCHBF:-}" bash -c "cd '${SRC}/CRISPRofiler' && '${CXX}' -std=c++11 -O2 -fopenmp ${CXXFLAGS} \
  main.cpp profiling.cpp guide_searching.cpp pam_searching.cpp pre_computation.cpp reading.cpp analysis.cpp \
  -o '${OPT}/searchBruteForce'"
ENR="${OPT}/Python_Scripts/Enrichment/enricher"
build Python_Scripts/Enrichment/enricher "${ENRICHER:-}" "${CXX}" -std=c++11 -O2 ${CXXFLAGS} \
  "${SRC}/Python_Scripts/Enrichment/enricher.cpp" -o "${ENR}" -lz
CZ="${ROOT}/bin/crispritz.py"

# --- the genome
G="${WORK}/g"
"${PYTHON}" "${SCRIPT_DIR}/four_allele/make_genome.py" "${G}" > "${WORK}/gen.log"
echo "--- $(cat "${WORK}/gen.log")"

echo "########## 1. ENRICH ##########"
# a PATH with no other enricher on it (an installed CRISPRitz puts one in the conda
# bin), so add-variants runs the port it is given: the compiled one, found beside
# crispritz.py as in a package; then enricher.py, once the compiled one is moved away
mkdir -p "${WORK}/pybin" && ln -s "${PYTHON}" "${WORK}/pybin/python3"
CLEAN_PATH="${WORK}/pybin:/usr/bin:/bin"
enrich() {  # <dir>
  mkdir -p "$1"
  ( cd "$1" && PATH="${CLEAN_PATH}" python3 "${CZ}" add-variants "${G}/vcf/" "${G}/ref/" > add.log 2> add.err ) \
    || fail "add-variants in $1: $(tail -3 "$1/add.err")"
  quiet "add-variants in $1" "$1/add.log" "$1/add.err"
}
enrich "${WORK}/enr_cpp"
mv "${ENR}" "${WORK}/enricher.moved"
enrich "${WORK}/enr_py"
mv "${WORK}/enricher.moved" "${ENR}"
GE="${WORK}/enr_cpp/variants_genome/SNPs_genome/ref_enriched"
diff -r "${GE}" "${WORK}/enr_py/variants_genome/SNPs_genome/ref_enriched" > /dev/null \
  || fail "the two enricher ports wrote different enriched genomes"
grep -q "each is coded Q, which the search reads as any base" "${WORK}/enr_cpp/add.log" \
  || fail "add-variants printed no four-allele note"
"${PYTHON}" - "${G}/ref" "${G}/vcf" "${GE}" <<'PY' || fail "the enriched genome is not the expected one"
import gzip, os, sys
ref_dir, vcf_dir, enr_dir = sys.argv[1:]
BIT = {"A": 1, "C": 2, "G": 4, "T": 8}
CODE = {3: "M", 5: "R", 9: "W", 6: "S", 10: "Y", 12: "K", 14: "B", 13: "D", 11: "H", 7: "V", 15: "Q"}
n_q = n_n = n_seq = 0
for name in sorted(os.listdir(ref_dir)):
    chrom = name[:-3]
    seq = "".join(open(os.path.join(ref_dir, name)).read().split("\n")[1:])
    want = list(seq)
    vcf = os.path.join(vcf_dir, "fourallele.%s.vcf.gz" % chrom)
    if os.path.exists(vcf):
        mask = {}
        for line in gzip.open(vcf, "rt"):
            if line.startswith("#"):
                continue
            f = line.split("\t")
            p = int(f[1]) - 1
            m = mask.get(p, BIT.get(seq[p], 0))
            for a in f[4].split(","):
                m |= BIT.get(a, 0)
            mask[p] = m
        for p, m in mask.items():
            if seq[p] != "N" and m in CODE:
                want[p] = CODE[m]
    want = "".join(want)
    got = "".join(open(os.path.join(enr_dir, chrom + ".enriched.fa")).read().split("\n")[1:])
    assert got == want, (chrom, [i for i in range(len(want)) if got[i:i + 1] != want[i]][:5])
    assert [i for i, c in enumerate(got) if c == "N"] == [i for i, c in enumerate(seq) if c == "N"], chrom
    n_q += got.count("Q"); n_n += got.count("N"); n_seq += 1
assert n_q > 1000, n_q
print("  ok: %d contigs enriched alike by both ports: %d Q, exactly where the reference base and "
      "ALTs cover four bases; %d N, exactly the reference's" % (n_seq, n_q, n_n))
PY

# --- PAM files
PAMS="${WORK}/pams"; mkdir -p "${PAMS}"
N20="NNNNNNNNNNNNNNNNNNNN"
printf '%s 3\n' "${N20}NGG" > "${PAMS}/NGG.txt"
printf '%s 3\n' "${N20}NRG" > "${PAMS}/NRG.txt"
printf '%s 8\n' "${N20}NNNNNNNN" > "${PAMS}/N8.txt"
printf '%s -4\n' "TTTV${N20}" > "${PAMS}/TTTV.txt"
printf '%s -4\n' "NNNN${N20}" > "${PAMS}/NNNN5.txt"
GUIDE="GACTTAGCCTAAGGTCTACG"
GUIDE2="GACTTAGCCTAAGGTCTAGC"
printf '%sNNN\n%sNNN\n' "${GUIDE}" "${GUIDE2}" > "${WORK}/guides3.txt"
printf 'NNNN%s\nNNNN%s\n' "${GUIDE}" "${GUIDE2}" > "${WORK}/guides5.txt"

echo "########## 2. INDEX ##########"
LIB="${WORK}/lib"; mkdir -p "${LIB}"
idxdir() { ls -d "${LIB}/genome_library/"*"_$2_$1"; }  # <name> <mb>
index() {  # <name> <pam> <mb>
  ( cd "${LIB}" && "${PYTHON}" "${CZ}" index-genome "$1" "${GE}/" "${PAMS}/$2.txt" -bMax "$3" -th 1 \
      > "${WORK}/$1.index.log" 2> "${WORK}/$1.index.err" ) || fail "index-genome $1: $(tail -3 "${WORK}/$1.index.err")"
  quiet "index-genome $1" "${WORK}/$1.index.log" "${WORK}/$1.index.err"
  grep -q "^FOUR-ALLELE (Q): [0-9]* target(s) hold Q" "${WORK}/$1.index.log" \
    || fail "index-genome $1: no four-allele note on stdout"
  ( cd "${WORK}" && "${PYTHON}" "${CZ}" index-info "$(idxdir "$1" "$3")" > "$1.info.json" 2> "$1.info.err" ) \
    || fail "index-info $1: $(head -3 "${WORK}/$1.info.err")"
  quiet "index-info $1" "${WORK}/$1.info.json" "${WORK}/$1.info.err"
  # a .bin is flagged exactly when one of its leaves stores Q: its guide region, the pad
  # bases before the first one off the contig or on an N, or its flank (the sites of
  # docs/INDEX_FORMAT.md section 3.1, enumerated here from the FASTA)
  "${PYTHON}" - "${WORK}/$1.info.json" "${GE}" "${PAMS}/$2.txt" "$3" <<'PY' || fail "index-info $1: the four-allele flags are wrong"
import json, os, sys
info, genome, pam_file, mb = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
BITS = {"A": 1, "C": 2, "G": 4, "T": 8, "R": 5, "Y": 10, "S": 6, "W": 9, "K": 12, "M": 3,
        "B": 14, "D": 13, "H": 11, "V": 7, "N": 15, "Q": 15}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M",
        "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N", "Q": "Q"}
motif, limit = open(pam_file).read().split()
limit = int(limit)
side, K = ("5", -limit) if limit < 0 else ("3", limit)
pam = motif[:K] if side == "5" else motif[-K:]
Gi = len(motif) - K
pamless = set(pam) == {"N"}
want = {}
for name in os.listdir(genome):
    with open(os.path.join(genome, name)) as fh:
        chrom = fh.readline()[1:].split()[0]
        seq = "".join(l.strip() for l in fh)
    n, holds = len(seq), False
    for T in (seq, "".join(COMP[c] for c in reversed(seq))):
        at = lambda x: T[x] if 0 <= x < n else None
        for g0 in range(n - Gi + 1):
            guide = T[g0:g0 + Gi]
            if "N" in guide:
                continue
            if side == "3":
                flank = [at(x) for x in range(g0 + Gi, g0 + Gi + K)]
                pad = [at(x) for x in range(g0 - 1, g0 - 1 - mb, -1)]
            else:
                flank = [at(x) for x in range(g0 - 1, g0 - 1 - K, -1)]
                pad = [at(x) for x in range(g0 + Gi, g0 + Gi + mb)]
            if pamless:
                if flank[0] is None:  # the flank base next to the guide must be on the contig
                    continue
            elif None in flank or any(b == "N" or not BITS[b] & BITS[p]
                                      for b, p in zip(flank, pam[::-1] if side == "5" else pam)):
                continue
            valid = 0
            while valid < mb and pad[valid] not in (None, "N"):
                valid += 1
            if "Q" in guide or "Q" in pad[:valid] or "Q" in flank:
                holds = True
                break
        if holds:
            break
    want[chrom] = holds
files = json.load(open(info))["files"]
assert all(f["format"] == 2 and "error" not in f for f in files.values()), files
got = {}
for name, f in files.items():
    chrom = name.split("_", 1)[1].rsplit("_", 1)[0]
    got[chrom] = got.get(chrom, False) or f["four_allele"]
bad = [c for c in got if got[c] != want[c]]
assert not bad, [(c, got[c], want[c]) for c in bad]
assert sum(got.values()) >= 20 and not got.get("chr4actlp", False) and not got.get("chr4aQnp", False), got
print("    %d contigs' .bin flagged four_allele, %d not, each as its leaves hold Q or not"
      % (sum(got.values()), len(got) - sum(got.values())))
PY
}
index q3a NGG 0
index q3b NGG 2
index q3c NRG 1
index q3d N8 2
index q5a TTTV 1
index q5b NNNN5 1
echo "  ok: 6 indexes; a .bin is flagged four_allele in its header exactly when it holds Q"

echo "########## 3. SEARCH ##########"
body() { tail -n +2 "$1" | LC_ALL=C sort; }
n_cmp=0; n_rows=0; n_bf=0; n_q_rows=0
search() {  # <index name> <index pam> <mb> <search pam> <prune 0|1> <mm> <bDNA> <bRNA>
  local name=$1 ipam=$2 mb=$3 spam=$4 prune=$5 mm=$6 bd=$7 br=$8 side=3 flag="" oflag=""
  [ "${spam}" = TTTV ] && side=5
  [ "${prune}" = 1 ] && flag="--prune-pam" && oflag="--prune"
  local tag="${name}.${spam}.p${prune}.mm${mm}.b${bd}${br}" idx; idx="$(idxdir "${name}" "${mb}")/"
  ( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${idx}" "${PAMS}/${spam}.txt" "${WORK}/guides${side}.txt" "${tag}" \
      -mm "${mm}" -bDNA "${bd}" -bRNA "${br}" -th 1 -r ${flag} > "${tag}.log" 2> "${tag}.err" ) \
    || fail "${tag}: search exited non-zero: $(tail -3 "${WORK}/${tag}.err")"
  quiet "${tag}" "${WORK}/${tag}.log" "${WORK}/${tag}.err"
  grep -q "^INDEX FORMAT: 2 x[0-9]* ([0-9]* holding Q, the four-allele symbol)" "${WORK}/${tag}.log" \
    || fail "${tag}: searchTST did not report the files holding Q"
  "${PYTHON}" "${SCRIPT_DIR}/boundary/oracle.py" --genome "${GE}" --guides "${WORK}/guides${side}.txt" \
    --pam "${PAMS}/${spam}.txt" --index-pam "${PAMS}/${ipam}.txt" --mb "${mb}" \
    --mm "${mm}" --bdna "${bd}" --brna "${br}" ${oflag} > "${WORK}/${tag}.oracle"
  diff <(body "${WORK}/${tag}.targets.txt") "${WORK}/${tag}.oracle" > "${WORK}/${tag}.diff" \
    || { head -12 "${WORK}/${tag}.diff"; fail "${tag}: indexed search != oracle ($(grep -c '^<' "${WORK}/${tag}.diff" || true) extra, $(grep -c '^>' "${WORK}/${tag}.diff" || true) missing)"; }
  "${PYTHON}" "${SCRIPT_DIR}/boundary/check_rows.py" "${GE}" "${WORK}/${tag}.targets.txt" > "${WORK}/${tag}.chk" \
    || { head -5 "${WORK}/${tag}.chk"; fail "${tag}: $(tail -1 "${WORK}/${tag}.chk")"; }
  ! tail -n +2 "${WORK}/${tag}.targets.txt" | cut -f3 | grep -q '[Nn]' || fail "${tag}: a row holds an N"
  local rows; rows="$(wc -l < "${WORK}/${tag}.oracle" | tr -d ' ')"
  [ "${rows}" -gt 0 ] || fail "${tag}: no rows (vacuous comparison)"
  n_cmp=$((n_cmp + 1)); n_rows=$((n_rows + rows))
  n_q_rows=$((n_q_rows + $(tail -n +2 "${WORK}/${tag}.targets.txt" | cut -f3 | grep -c Q || true)))
  if [ "${bd}" = 0 ] && [ "${br}" = 0 ]; then
    # every target the brute-force scan finds on the enriched FASTA, searchTST finds
    ( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${GE}/" "${PAMS}/${spam}.txt" "${WORK}/guides${side}.txt" "${tag}.bf" \
        -mm "${mm}" -th 1 -r > "${tag}.bf.log" 2> "${tag}.bf.err" ) || fail "${tag}: brute-force search failed"
    quiet "${tag}: brute force" "${WORK}/${tag}.bf.log" "${WORK}/${tag}.bf.err"
    diff <(body "${WORK}/${tag}.targets.txt") <(body "${WORK}/${tag}.bf.targets.txt") > "${WORK}/${tag}.bf.diff" \
      || { head -12 "${WORK}/${tag}.bf.diff"; fail "${tag}: indexed 0-bulge search != searchBruteForce ($(grep -c '^<' "${WORK}/${tag}.bf.diff" || true) only indexed, $(grep -c '^>' "${WORK}/${tag}.bf.diff" || true) only brute force)"; }
    n_bf=$((n_bf + 1))
  fi
}
for s in "0 0 0" "4 0 0"; do search q3a NGG 0 NGG 0 ${s}; done
for s in "0 0 0" "3 0 0" "1 1 0" "1 0 1" "2 1 1" "1 2 2"; do search q3b NGG 2 NGG 0 ${s}; done
for s in "0 0 0" "2 1 1"; do search q3c NRG 1 NRG 0 ${s}; done
for s in "0 0 0" "3 0 0" "1 1 0" "1 0 1" "2 1 1" "1 2 2"; do search q3d N8 2 NGG 1 ${s}; done
for s in "0 0 0" "2 1 1"; do search q3d N8 2 NRG 1 ${s}; done
for s in "0 0 0" "3 0 0" "1 1 1"; do search q5a TTTV 1 TTTV 0 ${s}; done
for s in "0 0 0" "1 1 1"; do search q5b NNNN5 1 TTTV 1 ${s}; done
[ "${n_q_rows}" -gt 1000 ] || fail "only ${n_q_rows} rows hold Q (vacuous test)"
echo "  ok: ${n_cmp} searches == oracle (${n_rows} rows, ${n_q_rows} of them through Q, none holding N); ${n_bf} 0-bulge searches == searchBruteForce"

# --prune-pam off the PAM-less indexes == the dedicated index of the same -bMax
for s in mm0.b00 mm3.b00 mm1.b10 mm1.b01 mm2.b11 mm1.b22; do
  cmp -s <(body "${WORK}/q3d.NGG.p1.${s}.targets.txt") <(body "${WORK}/q3b.NGG.p0.${s}.targets.txt") \
    || fail "--prune-pam NGG off the K=8 index != the NGG index at ${s}"
done
for s in mm0.b00 mm1.b11; do
  cmp -s <(body "${WORK}/q5b.TTTV.p1.${s}.targets.txt") <(body "${WORK}/q5a.TTTV.p0.${s}.targets.txt") \
    || fail "--prune-pam TTTV off the 5' K=4 index != the TTTV index at ${s}"
done
echo "  ok: --prune-pam off the PAM-less K=8 and 5' K=4 indexes == the dedicated NGG and TTTV indexes"

# every planted perfect site is found, at its position and strand, as a perfect match
"${PYTHON}" - "${G}/sites.tsv" "${WORK}/q3b.NGG.p0.mm0.b00.targets.txt" "${WORK}/q5a.TTTV.p0.mm0.b00.targets.txt" <<'PY' \
  || fail "a planted perfect site was not found"
import sys
sites, t3, t5 = sys.argv[1:]
found = set()
for path in (t3, t5):
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if not line.startswith("#") and c[1].startswith(("GACTTAGCCTAAGGTCTACG", "NNNNGACTTAGCCTAAGGTCTACG")) and c[7] == "0":
            found.add((c[3], c[6], int(c[4])))
want = [(c, st, int(p)) for c, side, st, p in (l.split() for l in open(sites))]
missing = [w for w in want if w not in found]
assert not missing, missing
print("  ok: all %d planted perfect sites found (Q in the guide, in the PAM, beside N runs, at contig ends)" % len(want))
PY

echo "########## 4. GAPS in the brute-force scan ##########"
( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${GE}/" "${PAMS}/NGG.txt" "${WORK}/guides3.txt" gaps \
    -mm 120 -th 1 -r > gaps.log 2> gaps.err ) || fail "brute force at -mm 120 failed"
quiet "brute force at -mm 120" "${WORK}/gaps.log" "${WORK}/gaps.err"
rows="$(tail -n +2 "${WORK}/gaps.targets.txt" | wc -l | tr -d ' ')"
[ "${rows}" -gt 1000 ] || fail "brute force at -mm 120: only ${rows} rows"
! tail -n +2 "${WORK}/gaps.targets.txt" | cut -f3 | grep -q '[Nn]' \
  || fail "brute force at -mm 120 reports $(tail -n +2 "${WORK}/gaps.targets.txt" | cut -f3 | grep -c '[Nn]') windows holding an N"
"${PYTHON}" "${SCRIPT_DIR}/boundary/check_rows.py" "${GE}" "${WORK}/gaps.targets.txt" > /dev/null \
  || fail "brute force at -mm 120: a row is not its contig's bases"
echo "  ok: searchBruteForce at -mm 120: ${rows} rows, none with an N in its guide"

echo "########## 5. OTHER BYTES in the genome ##########"
# an 'X' amid A/T (no PAM within reach, so no stored window holds it), then a perfect site
mkdir -p "${WORK}/gx"
"${PYTHON}" - "${WORK}/gx/chrX1.fa" "${GUIDE}" <<'PY'
import random, sys
random.seed(3)
seq = "".join(random.choice("ACGT") for _ in range(60)) + "AT" * 20 + "X" + "TA" * 20
seq += "TTTTT" + sys.argv[2] + "TGG" + "".join(random.choice("ACGT") for _ in range(60))
open(sys.argv[1], "w").write(">chrX1\n%s\n" % seq)
print(seq.index(sys.argv[2]), file=open(sys.argv[1] + ".site", "w"))
PY
site="$(cat "${WORK}/gx/chrX1.fa.site")"; rm "${WORK}/gx/chrX1.fa.site"
( cd "${LIB}" && "${PYTHON}" "${CZ}" index-genome gx "${WORK}/gx/" "${PAMS}/NGG.txt" -bMax 1 -th 1 > "${WORK}/gx.index.log" 2> "${WORK}/gx.index.err" ) \
  || fail "index-genome gx: $(tail -3 "${WORK}/gx.index.err")"
quiet "index-genome gx" "${WORK}/gx.index.log" "${WORK}/gx.index.err"
( cd "${WORK}" && "${PYTHON}" "${CZ}" search "$(idxdir gx 1)/" "${PAMS}/NGG.txt" "${WORK}/guides3.txt" gx \
    -mm 0 -bDNA 0 -bRNA 0 -th 1 -r > gx.log 2> gx.err ) || fail "search gx: $(tail -3 "${WORK}/gx.err")"
quiet "search gx" "${WORK}/gx.log" "${WORK}/gx.err"
( cd "${WORK}" && "${PYTHON}" "${CZ}" search "${WORK}/gx/" "${PAMS}/NGG.txt" "${WORK}/guides3.txt" gx.bf \
    -mm 0 -th 1 -r > gx.bf.log 2> gx.bf.err ) || fail "brute force gx: $(tail -3 "${WORK}/gx.bf.err")"
quiet "brute force gx" "${WORK}/gx.bf.log" "${WORK}/gx.bf.err"
for f in gx gx.bf; do
  got="$(tail -n +2 "${WORK}/${f}.targets.txt" | awk -F'\t' '$3 == "'"${GUIDE}"'TGG" {print $5 $7}')"
  [ "${got}" = "${site}+" ] || fail "${f}: the perfect site at ${site} after the 'X' was not found there (got '${got}')"
done
echo "  ok: a perfect site after an 'X' in the genome is found at its position, indexed and by brute force"

echo ""
echo "four-allele test: PASSED"
