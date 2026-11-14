#!/usr/bin/env bash
#
# `crispritz.py add-variants` must fail when an enricher fails, and when a chromosome it
# copies holds Q; and only then.
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
# and, with the shipped enricher, the chromosomes add-variants copies (those without a
# VCF):
#
#   5. one of them holds a 'q'           -> exit 1 before any enrichment, ERROR names
#                                           it and the position (Q is the four-allele
#                                           symbol, which the search reads as any base)
#   6. a README, a .DS_Store, a .fai and -> exit 0, nothing on stderr, every chromosome
#      a macOS ._chrC.fa (AppleDouble)      enriched or copied, the stray files ignored
#      beside the chromosomes               (2.9.1 wrote 'cp: cannot stat' on stderr,
#                                           copied ._chrC.fa as a chromosome or, with a
#                                           stray file listed first, copied none)
#   7. a genome folder with no FASTA     -> exit 1, ERROR names the folder
#
# Each run is under a timeout: a worker that exits instead of returning its status
# leaves pool.join() waiting forever, which must read as a failure, not a hang. The
# shipped enricher runs with a PATH holding only this python3 and the system tools, as
# add-variants prefers an enricher on PATH (an installed CRISPRitz puts one there).
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
if command -v timeout > /dev/null 2>&1; then TIMEOUT="$(command -v timeout) 120"
elif command -v gtimeout > /dev/null 2>&1; then TIMEOUT="$(command -v gtimeout) 120"
else TIMEOUT="$(command -v perl) -e alarm(shift);exec(@ARGV) 120"; fi
hung() { [ "$RC" = 124 ] || [ "$RC" = 142 ]; }

# the PATH of the shipped enricher's runs: this python3 (enricher.py runs under it) and
# the system tools, no enricher
mkdir -p "$WORK/pybin" && ln -s "$(command -v python3)" "$WORK/pybin/python3"
SHIPPED="$WORK/pybin:/usr/bin:/bin"

run() {  # $1 label, $2 PATH, $3 genome dir, rest: env assignments
    local label=$1 path=$2 genome=$3; shift 3
    rm -rf "$WORK/run_$label"; mkdir -p "$WORK/run_$label"
    # --debug: resolve the shipped ports in this checkout's sourceCode/, not a conda prefix
    ( cd "$WORK/run_$label" && env PATH="$path" "$@" \
        $TIMEOUT python3 "$CRISPRITZ" add-variants "$WORK/vcf" "$genome" -th 2 --debug \
        > out.txt 2> err.txt; echo $? > rc.txt )
    RC=$(cat "$WORK/run_$label/rc.txt"); ERR=$(cat "$WORK/run_$label/err.txt")
}

run ok "$WORK/bin:$PATH" "$WORK/genome" STUB_CHRB=ok
[ "$RC" = 0 ] || fail "all enrichers succeeded but add-variants exited $RC: $ERR"
case "$ERR" in *ERROR*) fail "all enrichers succeeded but stderr reports: $ERR" ;; esac
echo "  ok  every chromosome succeeds -> exit 0"

run exit3 "$WORK/bin:$PATH" "$WORK/genome" STUB_CHRB=exit3
hung && fail "add-variants hung on a failed enricher (timeout)"
[ "$RC" = 1 ] || fail "an enricher exited 3 but add-variants exited $RC"
case "$ERR" in *"enriching test.chrB.vcf.gz failed (exit code 3)"*) ;;
    *) fail "the ERROR line does not name chrB and its exit code: $ERR" ;; esac
case "$ERR" in *chrA*) fail "chrA succeeded but is reported: $ERR" ;; esac
echo "  ok  one chromosome exits 3 -> exit 1, reported by name and code"

run kill "$WORK/bin:$PATH" "$WORK/genome" STUB_CHRB=kill
hung && fail "add-variants hung on a killed enricher (timeout)"
[ "$RC" = 1 ] || fail "an enricher was killed but add-variants exited $RC"
case "$ERR" in *"enriching test.chrB.vcf.gz failed (signal SIGKILL)"*) ;;
    *) fail "the ERROR line does not name chrB and the signal: $ERR" ;; esac
echo "  ok  one chromosome is killed -> exit 1, reported by name and signal"

# the shipped enricher (no stub on PATH), on the same real input
run real "$SHIPPED" "$WORK/genome"
[ "$RC" = 0 ] || fail "the shipped enricher on valid input: add-variants exited $RC: $ERR"
enr=$(find "$WORK/run_real" -name 'chrA.enriched.fa' | head -1)
[ -s "$enr" ] || fail "the shipped enricher wrote no chrA.enriched.fa"
echo "  ok  the shipped enricher on valid input -> exit 0, enriched FASTA written"

# a chromosome without a VCF, which add-variants copies: chrC, holding a 'q' at 150
mkdir -p "$WORK/genome_q"
cp "$WORK/genome/chrA.fa" "$WORK/genome/chrB.fa" "$WORK/genome_q/"
python3 - "$WORK/genome_q/chrC.fa" <<'PY'
import random, sys
random.seed("chrC")
seq = [random.choice("ACGT") for _ in range(300)]
seq[149] = "q"
with open(sys.argv[1], "w") as f:
    f.write(">chrC\n%s\n%s\n" % ("".join(seq[:120]), "".join(seq[120:])))
PY
run q "$SHIPPED" "$WORK/genome_q"
hung && fail "add-variants hung on a reference holding q (timeout)"
[ "$RC" = 1 ] || fail "a copied chromosome holds a 'q' but add-variants exited $RC: $ERR"
case "$ERR" in *"ERROR: chrC holds 'Q' at 150: a reference FASTA must not hold Q, which is reserved for four-allele SNP positions"*) ;;
    *) fail "the ERROR line does not name chrC and position 150: $ERR" ;; esac
[ ! -e "$WORK/run_q/variants_genome" ] || fail "add-variants enriched before refusing chrC"
echo "  ok  a copied chromosome holds a 'q' -> exit 1 before any enrichment, chrC and 150 named"

# stray files beside the chromosomes: none is a chromosome, whatever the listing order
mkdir -p "$WORK/genome_stray"
cp "$WORK/genome/chrA.fa" "$WORK/genome/chrB.fa" "$WORK/genome_stray/"
sed 's/q/A/' "$WORK/genome_q/chrC.fa" > "$WORK/genome_stray/chrC.fa"
echo "chrA..chrC, one FASTA each" > "$WORK/genome_stray/README"
printf '\0\0\0\1Bud1' > "$WORK/genome_stray/.DS_Store"
printf 'chrA\t300\t6\t300\t301\n' > "$WORK/genome_stray/chrA.fa.fai"
# what macOS tar puts beside chrC.fa for its extended attributes, a quarantine one here
python3 -c 'import sys; open(sys.argv[1], "wb").write(b"\x00\x05\x16\x07\x00\x02\x00\x00Mac OS X"
    b"        \x00\x02ATTR\x00\x00com.apple.quarantine\x000083;00000000;Safari;")' "$WORK/genome_stray/._chrC.fa"
run stray "$SHIPPED" "$WORK/genome_stray"
hung && fail "add-variants hung on a genome folder with stray files (timeout)"
[ "$RC" = 0 ] || fail "stray files in the genome folder: add-variants exited $RC: $ERR"
[ -z "$ERR" ] || fail "stray files in the genome folder: add-variants wrote to stderr: $ERR"
out="$WORK/run_stray/variants_genome/SNPs_genome/genome_stray_enriched"
[ "$(ls -A "$out" | tr '\n' ' ')" = "chrA.enriched.fa chrB.enriched.fa chrC.enriched.fa " ] \
    || fail "the enriched genome is not chrA..chrC: $(ls -A "$out" | tr '\n' ' ')"
cmp -s "$out/chrC.enriched.fa" "$WORK/genome_stray/chrC.fa" || fail "chrC (no VCF) was not copied as it is"
for c in chrA chrB; do
    ! cmp -s "$out/$c.enriched.fa" "$WORK/genome_stray/$c.fa" || fail "$c was not enriched"
done
echo "  ok  a README, a .DS_Store, a .fai and a ._chrC.fa beside the chromosomes -> exit 0, nothing on stderr, ignored"

mkdir -p "$WORK/genome_none"
cp "$WORK/genome_stray/README" "$WORK/genome_stray/chrA.fa.fai" "$WORK/genome_none/"
run none "$SHIPPED" "$WORK/genome_none"
hung && fail "add-variants hung on a genome folder with no FASTA (timeout)"
[ "$RC" = 1 ] || fail "a genome folder with no FASTA: add-variants exited $RC: $ERR"
[ "$ERR" = "ERROR: $(cd "$WORK/genome_none" && pwd -P) holds no .fa or .fasta file" ] \
    || fail "a genome folder with no FASTA: not the ERROR line expected: $ERR"
echo "  ok  a genome folder with no FASTA -> exit 1, ERROR names it"

echo "add-variants exit status test: PASSED"
