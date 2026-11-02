#!/usr/bin/env bash
#
# Index format 2 (docs/INDEX_FORMAT.md): the header and trailer buildTST writes, the
# reader's refusals, and the corners the boundary suite does not reach.
#
#   1. HEADER     : every .bin buildTST writes starts with the magic "CZI\xfe", format 2,
#                   and the build parameters, and ends with the magic; crispritz.py
#                   index-info reads them (and reads a format-1 fixture as format 1).
#   2. MIXED      : a directory holding format-1 and format-2 files searches as the union
#                   of the two.
#   3. MAX BULGES : a format-2 index searched with another max_bulges on the command line
#                   decodes coordinates with its own, with a note on stdout (stderr is
#                   fatal to CRISPRme).
#   4. REFUSALS   : a cut-short file, a PAM side other than the index's, an unknown future
#                   format version, a header guide length other than the body's, and
#                   -bMax above 31 at build time each fail with an error naming the
#                   reason, not a crash or a silent wrong result.
#   5. GUIDE LENGTH : a search guide longer or shorter than the index's is searched, with
#                   a note on stdout naming both lengths (docs/INDEX_FORMAT.md 5.1); the
#                   boundary suite checks its rows.
#
# Configuration: CXX, CXXFLAGS, BUILDTST, SEARCHTST (prebuilt binaries), PYTHON,
# KEEP_WORKDIR.
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

BT="${BUILDTST:-${WORK}/buildTST}"
if [ -z "${BUILDTST:-}" ]; then
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC_DIR}/mainParallel.cpp" -o "${BT}"
fi
ST="${SEARCHTST:-${WORK}/searchTST}"
if [ -z "${SEARCHTST:-}" ]; then
  "${CXX}" -std=c++11 -O2 -fopenmp ${CXXFLAGS} "${SRC_DIR}/searchOnTST.cpp" \
    "${SRC_DIR}/detailedOutput.cpp" "${SRC_DIR}/convert.cpp" -I "${SRC_DIR}/include" -o "${ST}"
fi
hex() { od -An -tx1 | tr -d ' \n'; }

cd "${WORK}"
mkdir -p v2 && ( cd v2 && for fa in "${FIX}"/genome/*.fa; do "${BT}" "${fa}" "${FIX}/pam_NGG.txt" 1 2 > /dev/null; done )
search() {  # <index dir> <pam file> <max bulges> <out>; sets RC
  set +e
  "${ST}" "$1/" "${FIX}/guide3.txt" 2 1 1 "$2" "$4" r 1 "$3" > "$4.log" 2> "$4.err"
  RC=$?
  set -e
}

echo "########## 1. HEADER and trailer ##########"
for f in v2/*.bin; do
  # magic, version 2, header length 28 + 3, flags 0 (3' PAM, constrained), guide 20, PAM 3, mb 2, "NGG"
  [ "$(head -c 31 "${f}" | hex)" = "435a49fe020000001f00000000000000140000000300000002000000""4e4747" ] \
    || fail "${f}: unexpected header $(head -c 31 "${f}" | hex)"
  [ "$(tail -c 4 "${f}" | hex)" = "435a49fe" ] || fail "${f}: no end-of-index marker"
done
mkdir -p cli/bin && cp "${REPO_ROOT}/crispritz.py" cli/bin/
"${PYTHON}" cli/bin/crispritz.py index-info v2 > info2.json
"${PYTHON}" cli/bin/crispritz.py index-info "${FIX}/idx_NGG_2" > info1.json
"${PYTHON}" - info2.json info1.json <<'PY'
import json, sys
v2, v1 = (json.load(open(p)) for p in sys.argv[1:])
assert v2["formats"] == ["2"] and len(v2["files"]) == 6, v2
for f in v2["files"].values():
    assert f == dict(f, format=2, pam_at_start=False, pamless=False, guide_len=20, pam_len=3,
                     max_bulges=2, pam="NGG", complete=True) and f["num_leaves"] > 0, f
assert v1["formats"] == ["1"] and all(f["format"] == 1 for f in v1["files"].values()), v1
PY
echo "  ok: $(ls v2 | wc -l | tr -d ' ') .bin files carry the format-2 header and trailer; index-info reads formats 2 and 1"

echo "########## 2. MIXED directory ##########"
mkdir -p mixed p1 p2
cp "${FIX}"/idx_NGG_2/NGG_f[123]_1.bin mixed/ && cp "${FIX}"/idx_NGG_2/NGG_f[123]_1.bin p1/
cp v2/NGG_f[456]_1.bin mixed/ && cp v2/NGG_f[456]_1.bin p2/
for d in mixed p1 p2; do search "${d}" "${FIX}/pam_NGG.txt" 2 "s_${d}"; [ "${RC}" = 0 ] || fail "${d}: exit ${RC}"; done
grep -q "INDEX FORMAT: 1 (CRISPRitz <= 2.9.1) x3 2 x3" s_mixed.log || fail "mixed: formats not reported"
cmp -s <(tail -n +2 s_mixed.targets.txt | sort) <(cat <(tail -n +2 s_p1.targets.txt) <(tail -n +2 s_p2.targets.txt) | sort) \
  || fail "mixed directory != union of its parts"
echo "  ok: mixed directory == union of its format-1 and format-2 parts ($(tail -n +2 s_mixed.targets.txt | wc -l | tr -d ' ') rows)"

echo "########## 3. MAX BULGES from the header ##########"
search v2 "${FIX}/pam_NGG.txt" 2 s_mb2
search v2 "${FIX}/pam_NGG.txt" 3 s_mb3
[ "${RC}" = 0 ] && [ ! -s s_mb3.err ] || fail "max_bulges 3 on a max_bulges 2 index: exit ${RC}, stderr '$(cat s_mb3.err)'"
grep -q "NOTE: the index was built with max_bulges=2 (command line: 3)" s_mb3.log || fail "no max_bulges note"
cmp -s s_mb2.targets.txt s_mb3.targets.txt || fail "max_bulges 3 on a max_bulges 2 index changed the rows"
echo "  ok: command-line max_bulges 3 on a -bMax 2 index: same rows, a note on stdout"

echo "########## 4. REFUSALS ##########"
refused() {  # <label> <expected stderr fragment>
  [ "${RC}" = 1 ] && grep -q "$2" "$3" || fail "$1: expected exit 1 and '$2', got exit ${RC}: $(cat "$3")"
  echo "  ok: $1 -> $(head -1 "$3" | cut -c1-110)"
}
mkdir -p cut && cp v2/*.bin cut/ && f=cut/NGG_f1_1.bin && head -c $(( $(wc -c < "${f}") - 10 )) "${f}" > "${f}.t" && mv "${f}.t" "${f}"
search cut "${FIX}/pam_NGG.txt" 2 s_cut;   refused "cut-short file" "is cut short" s_cut.err
printf 'NNNNNNNNNNNNNNNNNNNNNNN -3\n' > n3_5p.txt
search v2 n3_5p.txt 2 s_side;              refused "5' search PAM on a 3' index" "was built for a PAM at the 3' end" s_side.err
mkdir -p v3 && cp v2/NGG_f1_1.bin v3/ && printf '\003' | dd of=v3/NGG_f1_1.bin bs=1 seek=4 conv=notrunc 2> /dev/null
search v3 "${FIX}/pam_NGG.txt" 2 s_v3;     refused "format version 3" "is index format 3" s_v3.err
mkdir -p glen && cp v2/NGG_f1_1.bin glen/ && printf '\025' | dd of=glen/NGG_f1_1.bin bs=1 seek=16 conv=notrunc 2> /dev/null
search glen "${FIX}/pam_NGG.txt" 2 s_glen; refused "header guide length 21, body 20" "header gives guide length 21, its body 20" s_glen.err
mkdir -p b32 && set +e && ( cd b32 && "${BT}" "${FIX}/genome/f1.fa" "${FIX}/pam_NGG.txt" 1 32 > /dev/null 2> ../b32.err ); RC=$?; set -e
refused "buildTST -bMax 32" "max bulges must be 0 to 31" b32.err
[ -z "$(ls b32)" ] || fail "buildTST -bMax 32 wrote files"

echo "########## 5. GUIDE LENGTH other than the index's ##########"
for g in 21 19; do
  set +e
  "${ST}" v2/ "${FIX}/guide${g}.txt" 2 1 1 "${FIX}/pam_NGG${g}.txt" "s_g${g}" r 1 2 > "s_g${g}.log" 2> "s_g${g}.err"
  RC=$?
  set -e
  [ "${RC}" = 0 ] && [ ! -s "s_g${g}.err" ] || fail "${g}-nt guides: exit ${RC}, stderr '$(cat "s_g${g}.err")'"
  [ "$(grep -c "^NOTE: the index was built for 20-nt guides and the search guides are ${g} nt" "s_g${g}.log")" = 1 ] \
    || fail "${g}-nt guides: no single guide-length note on stdout"
  [ "$(tail -n +2 "s_g${g}.targets.txt" | wc -l)" -gt 0 ] || fail "${g}-nt guides: no rows"
  echo "  ok: ${g}-nt guides on the 20-nt index: $(tail -n +2 "s_g${g}.targets.txt" | wc -l | tr -d ' ') rows, one note on stdout"
done
search v2 "${FIX}/pam_NGG.txt" 2 s_g20
! grep -q "^NOTE: the index was built for" s_g20.log || fail "a note for a guide of the index's own length"

echo "index-format test: PASSED"
