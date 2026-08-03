#!/usr/bin/env bash
# Verify the C++ enricher is byte-identical to the reference enricher.py.
#
# Builds a small synthetic chromosome FASTA + gzipped VCF containing a mix of
# variant types (biallelic SNP, multiallelic SNP, non-PASS -> skipped,
# insertion, deletion, multiallelic indel), then runs BOTH enricher.py and the
# compiled `enricher` on identical inputs in separate dirs (mode "yes") and
# diffs all four outputs.
#
# Prints "BYTE-IDENTICAL: PASS" only if every output diff is empty.
#
# Usage: verify_enricher.sh [path-to-enricher-binary]
#        (defaults to ./enricher, then /tmp/enricher_opt)

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENRICHER_PY="$SCRIPT_DIR/enricher.py"

# locate C++ binary
if [ "${1:-}" != "" ]; then
    ENRICHER_BIN="$1"
elif [ -x "$SCRIPT_DIR/enricher" ]; then
    ENRICHER_BIN="$SCRIPT_DIR/enricher"
elif [ -x "./enricher" ]; then
    ENRICHER_BIN="$(pwd)/enricher"
elif [ -x "/tmp/enricher_opt" ]; then
    ENRICHER_BIN="/tmp/enricher_opt"
else
    echo "ERROR: cannot find enricher binary. Build it or pass its path." >&2
    exit 1
fi
echo "Using C++ binary: $ENRICHER_BIN"
echo "Using reference : $ENRICHER_PY"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
echo "Work dir: $WORK"

# ---------------------------------------------------------------------------
# synthetic chromosome FASTA (a few hundred bp, mixed case to exercise upper())
# ---------------------------------------------------------------------------
CHR_FA="$WORK/chrTest.fa"
python3 - "$CHR_FA" <<'PY'
import sys, random
random.seed(1234)
seq = ''.join(random.choice('ACGTacgt') for _ in range(400))
with open(sys.argv[1], 'w') as f:
    f.write(">chrTest\n")
    # wrap at 60 cols like real fasta
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")
PY

# uppercased reference sequence, for building valid REF fields
REFSEQ="$(python3 - "$CHR_FA" <<'PY'
import sys
lines=open(sys.argv[1]).read().splitlines()
print(''.join(lines[1:]).upper())
PY
)"

# ---------------------------------------------------------------------------
# synthetic VCF: proper ## header + #CHROM, 5 samples, mixed variant types.
# REF bases are taken from the actual sequence so indel replacement matches.
# ---------------------------------------------------------------------------
VCF_PLAIN="$WORK/variants.vcf"
python3 - "$VCF_PLAIN" "$REFSEQ" <<'PY'
import sys
out = sys.argv[1]
ref = sys.argv[2]
def b(pos1):  # 1-based ref base
    return ref[pos1-1]

lines = []
lines.append("##fileformat=VCFv4.2")
lines.append('##FILTER=<ID=PASS,Description="All filters passed">')
lines.append('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">')
lines.append('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
samples = ["S1","S2","S3","S4","S5"]
lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))

def row(pos, rsid, r, a, info, gts, filt="PASS"):
    return "\t".join([
        "chrTest", str(pos), rsid, r, a, ".", filt, info, "GT", *gts])

# 1) biallelic PASS SNP at pos 50
p=50
lines.append(row(p, "rs100", b(p), ("T" if b(p)!="T" else "A"),
                 "AC=2;AF=0.25", ["0|1","0|0","1|1","0|0","1|0"]))

# 2) multiallelic SNP at pos 90 (two alt alleles, both single-base)
p=90
alts=[x for x in "ACGT" if x!=b(p)][:2]
lines.append(row(p, "rs200,rs201", b(p), ",".join(alts),
                 "AC=1,1;AF=0.1,0.2", ["0|1","0|2","0|0","2|1","1|0"]))

# 3) non-PASS variant at pos 120 (must be skipped)
p=120
lines.append(row(p, "rs300", b(p), ("G" if b(p)!="G" else "C"),
                 "AC=1;AF=0.05", ["0|1","0|0","0|0","0|0","0|0"], filt="LowQual"))

# 4) simple insertion at pos 150 (REF 1bp, ALT 2bp)
p=150
lines.append(row(p, "rs400", b(p), b(p)+"A",
                 "AC=1;AF=0.3", ["0|1","0|0","1|0","0|0","0|0"]))

# 5) simple deletion at pos 200 (REF 2bp, ALT 1bp)
p=200
r2=ref[p-1:p+1]
lines.append(row(p, "rs500", r2, r2[0],
                 "AC=1;AF=0.15", ["0|0","0|1","0|0","1|1","0|0"]))

# 6) multiallelic indel at pos 250 (a SNP alt + an insertion alt)
p=250
snpalt=("T" if b(p)!="T" else "A")
insalt=b(p)+"GG"
lines.append(row(p, "rs600,rs601", b(p), snpalt+","+insalt,
                 "AC=1,1;AF=0.12,0.22", ["0|1","0|2","1|2","0|0","2|0"]))

open(out,"w").write("\n".join(lines)+"\n")
PY

VCF_GZ="$WORK/variants.vcf.gz"
gzip -c "$VCF_PLAIN" > "$VCF_GZ"

# ---------------------------------------------------------------------------
# run reference (python) in its own dir
# ---------------------------------------------------------------------------
PYDIR="$WORK/py"
CPPDIR="$WORK/cpp"
mkdir -p "$PYDIR/SNPs_genome" "$CPPDIR/SNPs_genome"

run_one() {
    local dir="$1"; shift
    ( cd "$dir" && "$@" "$VCF_GZ" "$CHR_FA" "resultTest" "yes" "$VCF_PLAIN" ) \
        > "$dir/run.log" 2>&1
    return $?
}

echo "Running enricher.py ..."
run_one "$PYDIR" python3 "$ENRICHER_PY" || { echo "python run FAILED"; cat "$PYDIR/run.log"; exit 1; }
echo "Running C++ enricher ..."
run_one "$CPPDIR" "$ENRICHER_BIN" || { echo "C++ run FAILED"; cat "$CPPDIR/run.log"; exit 1; }

# ---------------------------------------------------------------------------
# diff all four outputs
# ---------------------------------------------------------------------------
FAIL=0
diff_out() {
    local rel="$1"
    if [ ! -e "$PYDIR/$rel" ] && [ ! -e "$CPPDIR/$rel" ]; then
        echo "  SKIP (neither produced): $rel"
        return
    fi
    if diff "$PYDIR/$rel" "$CPPDIR/$rel" > /dev/null 2>&1; then
        echo "  OK   : $rel"
    else
        echo "  DIFF : $rel"
        diff "$PYDIR/$rel" "$CPPDIR/$rel" | head -40
        FAIL=1
    fi
}

ENRICHED="SNPs_genome/resultTest_enriched/chrTest.enriched.fa"
JSON="SNPs_genome/my_dict_chrTest.json"
LOG="SNPs_genome/logchrTest.txt"
FAKE="fake_variants.vcf_chrTest/fakechrTest.fa"

echo "Comparing outputs:"
diff_out "$ENRICHED"
diff_out "$JSON"
diff_out "$LOG"
diff_out "$FAKE"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "BYTE-IDENTICAL: PASS"
    exit 0
else
    echo "BYTE-IDENTICAL: FAIL"
    exit 1
fi
