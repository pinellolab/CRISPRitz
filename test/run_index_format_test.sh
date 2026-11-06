#!/usr/bin/env bash
#
# Index format 2 (docs/INDEX_FORMAT.md): the header and trailer buildTST writes, the
# reader's refusals, and the corners the boundary suite does not reach.
#
#   1. HEADER     : every .bin buildTST writes starts with the magic "CZI\xfe", format 2,
#                   and the build parameters, and ends with the magic; crispritz.py
#                   index-info reads them (and reads a format-1 fixture as format 1),
#                   and its INDEX_FORMAT is the format buildTST writes.
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
#   6. INDEX-INFO : index-info refuses (exit 1, an "error" entry) each of 18 .bin files
#                   damaged in one way that searchTST refuses, for the same reason (and,
#                   not as root, an unreadable one), and a directory whose .bin files
#                   are named for more than one PAM (a macOS "._" file beside the
#                   index), and accepts the good ones; the reference reader of
#                   docs/INDEX_FORMAT.md section 7 returns what index-info does for every
#                   file.
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
"${PYTHON}" - info2.json info1.json "${REPO_ROOT}/crispritz.py" "${SRC_DIR}/include/indexFormat.h" <<'PY'
import json, re, sys
v2, v1 = (json.load(open(p)) for p in sys.argv[1:3])
assert v2["formats"] == ["2"] and len(v2["files"]) == 6 and v2["pams"] == v1["pams"] == ["NGG"], v2
for f in v2["files"].values():
    assert f == dict(f, format=2, pam_at_start=False, pamless=False, four_allele=False, guide_len=20,
                     pam_len=3, max_bulges=2, pam="NGG", complete=True) and f["num_leaves"] > 0, f
    assert "error" not in f, f
assert v1["formats"] == ["1"] and all(f["format"] == 1 and "error" not in f
                                      for f in v1["files"].values()), v1
# INDEX_FORMAT, read from crispritz.py as CRISPRme reads VERSION, is the format buildTST
# writes (CZI_FORMAT_VERSION, and the files above), and index-info reports it
src = open(sys.argv[3]).read()
fmt = [int(m.group(1)) for m in re.finditer(r"^INDEX_FORMAT\s*=\s*(\d+)\s*$", src, re.M)]
czi = re.findall(r"CZI_FORMAT_VERSION = (\d+);", open(sys.argv[4]).read())
assert fmt == [2] and czi == ["2"] and v2["writes_format"] == v1["writes_format"] == 2, (fmt, czi)
PY
echo "  ok: $(ls v2 | wc -l | tr -d ' ') .bin files carry the format-2 header and trailer; index-info reads formats 2 and 1; INDEX_FORMAT = 2 = what buildTST writes"

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

echo "########## 6. INDEX-INFO refuses what searchTST refuses ##########"
# Each case is one .bin (v2/NGG_f1_1.bin, or for the last one the 2.9.1-built fixture's)
# damaged in one way, alone in a directory of its own. searchTST must refuse it (exit 1,
# its reason on stderr), and so must index-info (exit 1, the file and the reason on
# stderr, an "error" entry in its JSON).
"${PYTHON}" - "${WORK}" "${FIX}/idx_NGG_2/NGG_f1_1.bin" <<'PY'
import os, struct, sys
work, fmt1 = sys.argv[1:]
data = open(os.path.join(work, "v2", "NGG_f1_1.bin"), "rb").read()
def put(case, d):
    os.makedirs(os.path.join(work, "bad", case))
    open(os.path.join(work, "bad", case, "NGG_f1_1.bin"), "wb").write(bytes(d))
def at(off, val, d=data):  # an int32 of the header (or body) set to val
    d = bytearray(d)
    struct.pack_into("=i", d, off, val)
    return d
put("header_len_40", at(8, 40))      # header_len (the reviewer's case)
put("flag_bit_3", at(12, 8))         # a reserved flag bit (bit 2 marks the four-allele symbol Q)
put("guide_len_0", at(16, 0))
put("pam_len_0", at(20, 0))
put("max_bulges_32", at(24, 32))
put("max_bulges_-1", at(24, -1))
put("version_1", at(4, 1))           # the magic, then a version other than 2
put("version_3", at(4, 3))
put("header_cut", data[:20])         # the file ends inside the fixed header
put("pam_cut", data[:29])            # ... inside the header's PAM
put("cut_10", data[:-10])            # the end marker gone with the tree's end
put("no_end_marker", data[:-4] + b"\0\0\0\0")
put("header_and_marker_only", data[:31] + data[-4:])
put("guide_len_21", at(16, 21))      # header guide length 21, body 20
d = bytearray(data); d[28:31] = b"NAG"
put("pam_NAG", d)                    # header PAM other than the file name's
put("leaves_0", at(31, 0))           # numLeaves, at header_len
put("tiny", data[:2])
put("format1_leaves_0", at(0, 0, open(fmt1, "rb").read()))
PY
n_ref=0
while IFS='|' read -r case engine info; do
  [ -n "${case}" ] || continue
  search "bad/${case}" "${FIX}/pam_NGG.txt" 2 "s_bad_${case}"
  [ "${RC}" = 1 ] && grep -qF -- "${engine}" "s_bad_${case}.err" \
    || fail "${case}: searchTST: expected exit 1 and '${engine}', got exit ${RC}: $(cat "s_bad_${case}.err")"
  set +e
  "${PYTHON}" cli/bin/crispritz.py index-info "bad/${case}" > "ii_${case}.json" 2> "ii_${case}.err"
  RC=$?
  set -e
  [ "${RC}" = 1 ] && grep -qF -- "NGG_f1_1.bin: ${info}" "ii_${case}.err" \
    || fail "${case}: index-info: expected exit 1 and 'NGG_f1_1.bin: ${info}', got exit ${RC}: $(cat "ii_${case}.err")"
  "${PYTHON}" -c 'import json, sys; f = json.load(open(sys.argv[1]))["files"]["NGG_f1_1.bin"]; assert sys.argv[2] in f["error"], f' \
    "ii_${case}.json" "${info}" || fail "${case}: index-info's JSON has no such error entry"
  n_ref=$((n_ref + 1))
  echo "  ok: ${case}: searchTST and index-info refuse it ($(head -1 "ii_${case}.err" | sed 's/.*NGG_f1_1.bin: //' | cut -c1-70))"
done <<'EOF'
header_len_40|has a malformed index header (version 2, flags 0, PAM length 3, header length 40)|malformed index header (version 2, flags 0, PAM length 3, header length 40,
flag_bit_3|has a malformed index header (version 2, flags 8,|malformed index header (version 2, flags 8,
guide_len_0|has a malformed index header (version 2, flags 0, PAM length 3, header length 31)|malformed index header (version 2, flags 0, PAM length 3, header length 31, guide length 0, max bulges 2)
pam_len_0|has a malformed index header (version 2, flags 0, PAM length 0,|malformed index header (version 2, flags 0, PAM length 0,
max_bulges_32|has a malformed index header (version 2, flags 0, PAM length 3, header length 31)|malformed index header (version 2, flags 0, PAM length 3, header length 31, guide length 20, max bulges 32)
max_bulges_-1|has a malformed index header (version 2, flags 0, PAM length 3, header length 31)|malformed index header (version 2, flags 0, PAM length 3, header length 31, guide length 20, max bulges -1)
version_1|has a malformed index header (version 1,|malformed index header (version 1,
version_3|is index format 3, which is newer than this CRISPRitz reads|index format 3, newer than this CRISPRitz reads
header_cut|has a cut-short index header|cut-short index header
pam_cut|has a cut-short index header|cut-short index header
cut_10|is cut short: it does not end in the end-of-index marker|cut short: it does not end in the end-of-index marker
no_end_marker|is cut short: it does not end in the end-of-index marker|cut short: it does not end in the end-of-index marker
header_and_marker_only|is not a valid CRISPRitz index (.bin) file (numLeaves=-28747197)|cut short: no leaf count and guide length after its header
guide_len_21|has a malformed index: its header gives guide length 21, its body 20|malformed index: its header gives guide length 21, its body 20
pam_NAG|holds a 3-nt PAM/flank region ('NAG') but is named for 'NGG'|holds the PAM/flank region 'NAG' but is named for 'NGG'
leaves_0|is not a valid CRISPRitz index (.bin) file (numLeaves=0)|not a valid CRISPRitz index (.bin) file (numLeaves=0)
tiny|is too short to be a CRISPRitz index (.bin) file|too short to be a CRISPRitz index (.bin) file
format1_leaves_0|is not a valid CRISPRitz index (.bin) file (numLeaves=0)|not a valid CRISPRitz index (.bin) file (numLeaves=0)
EOF
[ "${n_ref}" = 18 ] || fail "${n_ref} damaged files tested, expected 18"
# a .bin that cannot be opened (root reads any file, so not as root)
if [ "$(id -u)" != 0 ]; then
  mkdir -p unreadable && cp v2/NGG_f1_1.bin unreadable/ && chmod 000 unreadable/NGG_f1_1.bin
  search unreadable "${FIX}/pam_NGG.txt" 2 s_unreadable
  [ "${RC}" = 1 ] && grep -qF "NGG_f1_1.bin' cannot be opened" s_unreadable.err \
    || fail "an unreadable .bin: searchTST: exit ${RC}: $(cat s_unreadable.err)"
  set +e
  "${PYTHON}" cli/bin/crispritz.py index-info unreadable > /dev/null 2> ii_unreadable.err
  RC=$?
  set -e
  chmod 644 unreadable/NGG_f1_1.bin
  [ "${RC}" = 1 ] && grep -qF "NGG_f1_1.bin: cannot be opened" ii_unreadable.err \
    || fail "an unreadable .bin: index-info: exit ${RC}: $(cat ii_unreadable.err)"
  echo "  ok: unreadable: searchTST and index-info refuse it ($(sed 's/.*NGG_f1_1.bin: //' ii_unreadable.err))"
else
  echo "  (an unreadable .bin: not tested as root)"
fi
# .bin files named for more than one PAM: a macOS "._" AppleDouble file copied along with
# an index (its first int, read as a format-1 leaf count, is positive, so the file alone
# passes). searchTST reads every .bin as an index of the PAM that names the first one it
# lists; it refused such a directory ("this index stores a 1-nt PAM/flank region ('.')")
# on Linux amd64 and arm64. It is not run here: listed the other way round, it would take
# the 163-byte file's first int, 118,883,584, for a leaf count and allocate as many leaves.
mkdir -p appledouble && cp v2/NGG_f1_1.bin appledouble/
"${PYTHON}" -c 'import sys; open(sys.argv[1], "wb").write(b"\x00\x05\x16\x07\x00\x02\x00\x00Mac OS X        " + bytes(139))' \
  appledouble/._NGG_f1_1.bin
set +e
"${PYTHON}" cli/bin/crispritz.py index-info appledouble > ii_ad.json 2> ii_ad.err
RC=$?
set -e
[ "${RC}" = 1 ] && grep -qF "the .bin files are named for more than one PAM ('.', 'NGG')" ii_ad.err \
  || fail "a ._ AppleDouble file beside the index: index-info: exit ${RC}: $(cat ii_ad.err)"
"${PYTHON}" -c 'import json, sys; d = json.load(open(sys.argv[1])); assert d["pams"] == [".", "NGG"] and d["files"]["._NGG_f1_1.bin"]["format"] == 1, d' \
  ii_ad.json || fail "a ._ AppleDouble file beside the index: index-info's JSON"
echo "  ok: ._NGG_f1_1.bin (AppleDouble) beside NGG_f1_1.bin: index-info refuses the directory (.bin files named for '.' and 'NGG')"
# what index-info accepts, searchTST searches: the good directories, both formats
for d in v2 mixed "${FIX}/idx_NGG_2"; do
  "${PYTHON}" cli/bin/crispritz.py index-info "${d}" > /dev/null 2> ii_good.err || fail "index-info refuses ${d}: $(cat ii_good.err)"
  [ ! -s ii_good.err ] || fail "index-info wrote to stderr on ${d}"
done
# docs/INDEX_FORMAT.md section 7's reference reader returns what index-info does for every
# file above
n_doc="$("${PYTHON}" - "${REPO_ROOT}/docs/INDEX_FORMAT.md" "${PYTHON}" cli/bin/crispritz.py \
  v2 mixed "${FIX}/idx_NGG_2" bad/* <<'PY'
import json, re, subprocess, sys
doc, py, cz, dirs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
block = re.search(r"^## 7\..*?^```python\n(.*?)^```", open(doc).read(), re.S | re.M).group(1)
ns = {}
exec(compile(block, "INDEX_FORMAT.md section 7", "exec"), ns)
n = 0
for d in dirs:
    got = json.loads(subprocess.run([py, cz, "index-info", d], capture_output=True, text=True).stdout)
    for name, entry in got["files"].items():
        ref = ns["read_index_header"]("%s/%s" % (d, name))
        assert ref == entry, (d, name, ref, entry)
        n += 1
print(n)
PY
)" || fail "docs/INDEX_FORMAT.md section 7's reader differs from index-info"
[ "${n_doc}" = 36 ] || fail "the section 7 reader was compared on ${n_doc} files, expected 36"
echo "  ok: ${n_ref} damaged files refused by both, for the same reason; the good directories accepted; docs/INDEX_FORMAT.md section 7's reader == index-info on ${n_doc} files"

echo "index-format test: PASSED"
