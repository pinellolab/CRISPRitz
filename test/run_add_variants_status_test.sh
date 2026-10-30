#!/usr/bin/env bash
#
# `crispritz.py add-variants` must fail when an enricher fails.
#
# add-variants runs one enricher per chromosome in a multiprocessing Pool. Nothing
# read the workers' results, so an enricher that crashed, was killed (OOM) or ran
# out of disk left add-variants at exit 0, with that chromosome's enriched genome,
# SNP dict, indel log and contigs missing or truncated. CRISPRme's index build
# trusts that exit status and goes on to index whatever is there.
#
# A stub `enricher` on PATH (which add-variants prefers over the shipped ports)
# stands in for the real one, so each outcome is exact:
#
#   1. every chromosome succeeds         -> exit 0, no ERROR line
#   2. one chromosome exits 3            -> exit 1, ERROR names it and "exit code 3"
#   3. one chromosome is killed (KILL)   -> exit 1, ERROR names it and the signal
#   4. the shipped enricher on real input -> exit 0 and the enriched FASTA exists
#
# Each run is under a timeout: a worker that exits instead of returning its status
# leaves pool.join() waiting forever, which must read as a failure, not a hang.
#
# Usage: test/run_add_variants_status_test.sh

set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CRISPRITZ="$REPO/crispritz.py"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

# two chromosomes, one VCF each (the VCF name must carry ".chrN.")
mkdir -p "$WORK/genome" "$WORK/vcf" "$WORK/bin"
for c in chrA chrB; do
    python3 - "$WORK/genome/$c.fa" "$c" <<'PY'
import random, sys
random.seed(sys.argv[2])
seq = "".join(random.choice("ACGT") for _ in range(300))
with open(sys.argv[1], "w") as f:
    f.write(">%s\n%s\n" % (sys.argv[2], seq))
PY
    ref=$(sed -n 2p "$WORK/genome/$c.fa" | cut -c100)
    alt=$([ "$ref" = A ] && echo C || echo A)
    printf '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n%s\t100\trs1\t%s\t%s\t.\tPASS\tAF=0.5\tGT\t0|1\n' \
        "$c" "$ref" "$alt" | gzip > "$WORK/vcf/test.$c.vcf.gz"
done

# stub enricher: argv = <vcf> <fasta> <genome dir name> <yes|no> <vcf dir>
cat > "$WORK/bin/enricher" <<'SH'
#!/usr/bin/env bash
case "$2" in
    *chrB*) case "${STUB_CHRB:-ok}" in
                exit3) exit 3 ;;
                kill)  kill -KILL $$ ;;
            esac ;;
esac
exit 0
SH
chmod +x "$WORK/bin/enricher"

# a hang must fail the test, not stall it: GNU timeout (exit 124) where present,
# else a perl alarm (SIGALRM, exit 142)
if command -v timeout > /dev/null 2>&1; then TIMEOUT="timeout 120"
elif command -v gtimeout > /dev/null 2>&1; then TIMEOUT="gtimeout 120"
else TIMEOUT="perl -e alarm(shift);exec(@ARGV) 120"; fi
hung() { [ "$RC" = 124 ] || [ "$RC" = 142 ]; }

run() {  # $1 label, $2 PATH prefix ("" = none), rest: env assignments
    local label=$1 prefix=$2; shift 2
    rm -rf "$WORK/run_$label"; mkdir -p "$WORK/run_$label"
    # --debug: resolve the shipped ports in this checkout's sourceCode/, not a conda prefix
    ( cd "$WORK/run_$label" && env ${prefix:+PATH="$prefix:$PATH"} "$@" \
        $TIMEOUT python3 "$CRISPRITZ" add-variants "$WORK/vcf" "$WORK/genome" -th 2 --debug \
        > out.txt 2> err.txt; echo $? > rc.txt )
    RC=$(cat "$WORK/run_$label/rc.txt"); ERR=$(cat "$WORK/run_$label/err.txt")
}

run ok "$WORK/bin" STUB_CHRB=ok
[ "$RC" = 0 ] || fail "all enrichers succeeded but add-variants exited $RC: $ERR"
case "$ERR" in *ERROR*) fail "all enrichers succeeded but stderr reports: $ERR" ;; esac
echo "  ok  every chromosome succeeds -> exit 0"

run exit3 "$WORK/bin" STUB_CHRB=exit3
hung && fail "add-variants hung on a failed enricher (timeout)"
[ "$RC" = 1 ] || fail "an enricher exited 3 but add-variants exited $RC"
case "$ERR" in *"enriching test.chrB.vcf.gz failed (exit code 3)"*) ;;
    *) fail "the ERROR line does not name chrB and its exit code: $ERR" ;; esac
case "$ERR" in *chrA*) fail "chrA succeeded but is reported: $ERR" ;; esac
echo "  ok  one chromosome exits 3 -> exit 1, reported by name and code"

run kill "$WORK/bin" STUB_CHRB=kill
hung && fail "add-variants hung on a killed enricher (timeout)"
[ "$RC" = 1 ] || fail "an enricher was killed but add-variants exited $RC"
case "$ERR" in *"enriching test.chrB.vcf.gz failed (signal SIGKILL)"*) ;;
    *) fail "the ERROR line does not name chrB and the signal: $ERR" ;; esac
echo "  ok  one chromosome is killed -> exit 1, reported by name and signal"

# the shipped enricher (no stub on PATH), on the same real input
if ! command -v enricher > /dev/null 2>&1; then
    run real ""
    [ "$RC" = 0 ] || fail "the shipped enricher on valid input: add-variants exited $RC: $ERR"
    enr=$(find "$WORK/run_real" -name 'chrA.enriched.fa' | head -1)
    [ -s "$enr" ] || fail "the shipped enricher wrote no chrA.enriched.fa"
    echo "  ok  the shipped enricher on valid input -> exit 0, enriched FASTA written"
else
    echo "  --  skipped the shipped-enricher run: an 'enricher' is already on PATH"
fi

echo "add-variants exit status test: PASSED"
