#!/usr/bin/env bash
# Verify the C++ enricher is byte-identical to the reference enricher.py.
#
# Builds a small synthetic chromosome FASTA + gzipped VCF containing a mix of
# variant types (biallelic SNP, multiallelic SNP, non-PASS -> skipped,
# insertion, deletion, multiallelic indel, and a site with TWO carried indel
# alleles), then runs BOTH enricher.py and the
# compiled `enricher` on identical inputs in separate dirs (mode "yes") and
# diffs all four outputs.
#
# Cases 9-16 cover the SNP+indel COMBINATIONS: SNPs inside an indel's flanking
# window (same haplotype / other haplotype / other sample, window edges), a
# deletion spanning SNPs, indels with overlapping windows, multiallelic records
# mixing SNP and indel ALTs, unphased / missing / haploid genotypes, indels at
# the chromosome start (window floor) and end, a REF that does not match the
# FASTA, and same-POS SNP/indel orders beyond case 8. Case 17 covers multiallelic
# SNP records whose ALTs include a non-A/C/G/T member ('*', 'N', a lowercase base)
# or whose FASTA base is N. Cases 18-30 cover multi-allelic sites SPLIT into one
# record per ALT, as `bcftools norm -m -any` writes them: two and three ALTs in
# every record order, the same sites as one record, the all-four-bases rule, split
# records that are not adjacent, haploid / unphased / missing genotypes, SNP + indel
# at one position and inside an indel's window, a split record beside a multi-ALT
# record, a non-PASS split record, the same ALT twice, AF spellings and AF_<x>
# fields, a reference IUPAC base and N, '*', the chromosome edges, and split records
# phased so that one haplotype slot holds two ALTs. Cases 31-33 cover each record's
# INFO/AF read by its key wherever it stands, the same REF>ALT on two records merged
# into one dict entry, and ALTs on a reference N and on reference bytes that are no
# base or IUPAC code (X, '-'). Besides the four outputs, the file
# set and the run log (stdout+stderr, e.g. the skip warnings and the four-base note)
# must match.
#
# Agreement is not correctness, so the outputs are also checked against what they
# must be: the enriched genome against an independent reconstruction from the FASTA
# and the VCF (all four bases are Q, the search's any-base symbol, and an N stays an
# N, as does a byte that is no base or IUPAC code: the enriched genome holds an N
# exactly where the reference does -- the N-count invariant), the dict at cases 24
# and 31-33 (merged carriers, AF by key, the records at a reference N, X or '-', a
# record listing three SNP ALTs) and the indel log's AF at case 31.
# A reference holding 'Q' must be refused by both ports alike.
#
# Prints "BYTE-IDENTICAL: PASS" only if every check passes.
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
# synthetic chromosome FASTA (2,800 bp, mixed case to exercise upper()). The
# first 400 bp -- all that cases 1-8 touch -- are the same draws as when the
# chromosome was 400 bp long, the first 1,800 -- cases 9-17 -- the same as when it
# was 1,800 bp long, and the first 2,400 -- cases 18-30 -- the same as when it was
# 2,400 bp long; the rest holds cases 31-33 (case 14 follows the end).
# ---------------------------------------------------------------------------
CHR_FA="$WORK/chrTest.fa"
python3 - "$CHR_FA" <<'PY'
import sys, random
random.seed(1234)
seq = [random.choice('ACGTacgt') for _ in range(2800)]
# case 15: a hard-masked N run and an IUPAC base in the reference itself
for i in range(1610, 1621):
    seq[i - 1] = 'N'
seq[1650 - 1] = 'R'
# case 26: another reference IUPAC base, for split records
seq[2200 - 1] = 'Y'
# case 33: bytes that are no base or IUPAC code (RepeatMasker -x masks with X)
seq[2420 - 1] = 'X'
seq[2425 - 1] = '-'
seq = ''.join(seq)
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
lines.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
samples = ["S1","S2","S3","S4","S5"]
lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))

def row(pos, rsid, r, a, info, gts, filt="PASS", fmt="GT"):
    return "\t".join([
        "chrTest", str(pos), rsid, r, a, ".", filt, info, fmt, *gts])

def other(pos1, k=0):  # the k-th base that differs from the reference base
    return [x for x in "ACGT" if x != b(pos1)][k]

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

# 7) TWO CARRIED INDEL ALLELES at one site (pos 330): REF T -> T+A and T+AA, the
#    shape of the real chr22:42762462 T -> TA / TAA / TAAA site. Case 6 above is a
#    SNP alt plus ONE insertion, so its `indels` list has a single entry and the
#    per-allele loop runs once -- it cannot catch an allele being built on top of
#    the previous one. This case makes that loop run twice, with S5 carrying both
#    alleles so the carrier scan is exercised too.
p=330
lines.append(row(p, "rs700,rs701", b(p), b(p)+"A," + b(p)+"AA",
                 "AC=2,1;AF=0.3,0.1", ["0|1","0|2","1|1","0|0","2|1"]))

# 8) A SNP and an INDEL at the SAME POS, as SEPARATE BIALLELIC RECORDS (pos 300).
#    This is NOT case 6: there the SNP alt and the insertion alt share one
#    multiallelic record, while here they are two rows -- which is how every
#    CRISPRme-supplied VCF is normalised (0 of 1,998,332 real chr22 records carry
#    a comma in ALT). Measured on a real 120,000-record chr22 slice: the two ports
#    disagree at exactly the 176 positions of this shape and nowhere else. The
#    Python reference keeps the SNP's IUPAC code (C/A -> M); the C++ port writes
#    the plain reference base, so the SNP disappears from the enriched genome --
#    and a SNP co-occurring with an indel is precisely what the SNP+indel feature
#    exists to find. None of the cases above could catch it.
p=300
lines.append(row(p, "rs800", b(p), ("T" if b(p)!="T" else "A"),
                 "AC=1;AF=0.2", ["0|1","0|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs801", b(p), b(p)+"CC",
                 "AC=1;AF=0.1", ["0|0","0|1","0|0","0|0","0|0"]))

# ---- cases 9-16: SNP + indel COMBINATIONS ------------------------------------
# The enricher cuts each indel's contig from the PLAIN reference: a 25 bp left
# flank, REF->ALT at offset 25, a 27 bp right flank, i.e. the 1-based window
# [POS-25, POS+26+len(REF)]. A SNP in that window goes into the enriched genome
# and the dict but never into the contig (CRISPRme overlays it later, in
# overlay_indel_snps.py, and judges cis/trans from the genotypes), so here the
# two ports must agree on all of: where the SNP lands, what the contig holds, and
# which carriers each output lists.

# 9) SNPs INSIDE an indel's flanking window: carriers on the SAME haplotype as
#    the indel (S1 1|0 + 1|0), the OTHER haplotype (S2 0|1 + 1|0) and a
#    DIFFERENT sample (S4); SNPs on the first/last window base and one base
#    outside each edge; around an insertion (450) and a 4-bp deletion (550).
p=450
lines.append(row(p, "rs900", b(p), b(p)+"AG", "AC=4;AF=0.4", ["1|0","0|1","1|1","0|0","0|0"]))
for s, gts in ((424, ["1|0","0|0","0|0","0|1","0|0"]), (425, ["1|0","1|0","0|0","0|0","0|0"]),
               (440, ["1|0","1|0","0|0","0|1","0|0"]), (451, ["0|1","0|1","1|0","0|0","0|0"]),
               (465, ["1|0","0|0","0|1","1|1","0|0"]), (477, ["1|0","1|0","0|0","0|0","0|0"]),
               (478, ["1|0","0|0","0|0","0|0","1|0"])):
    lines.append(row(s, "rsW%d" % s, b(s), other(s), "AC=1;AF=0.1", gts))
p=550
lines.append(row(p, "rs901", ref[p-1:p+3], b(p), "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
for s, gts in ((524, ["1|0","0|0","0|0","0|1","0|0"]), (525, ["1|0","1|0","0|0","0|0","0|0"]),
               (540, ["1|0","1|0","0|0","0|1","0|0"]), (554, ["1|0","0|1","0|0","0|0","0|0"]),
               (570, ["0|1","1|0","0|0","1|1","0|0"]), (580, ["1|0","0|0","0|0","0|0","0|0"]),
               (581, ["1|0","0|0","0|0","0|0","0|1"])):
    lines.append(row(s, "rsW%d" % s, b(s), other(s), "AC=1;AF=0.1", gts))

# 10) A DELETION SPANNING SNP positions: a 6-bp deletion at 650 with SNPs at its
#     anchor (650), inside the deleted span (651, 653, 655) and just after it
#     (656); and a 3-bp deletion at 700 with an insertion starting inside its
#     span (701) plus a SNP at 702.
p=650
lines.append(row(p, "rs1000", ref[p-1:p+5], b(p), "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
for s, gts in ((650, ["0|1","0|0","1|0","0|0","0|0"]), (651, ["0|1","0|1","0|0","0|0","0|0"]),
               (653, ["0|1","1|0","0|0","1|0","0|0"]), (655, ["1|0","0|0","0|0","0|0","0|0"]),
               (656, ["1|0","0|1","0|0","0|0","0|0"])):
    lines.append(row(s, "rsX%d" % s, b(s), other(s), "AC=1;AF=0.1", gts))
p=700
lines.append(row(p, "rs1001", ref[p-1:p+2], b(p), "AC=1;AF=0.1", ["1|0","0|0","0|0","0|0","0|0"]))
lines.append(row(701, "rs1002", b(701), b(701)+"T", "AC=2;AF=0.2", ["0|1","1|0","0|0","0|0","0|0"]))
lines.append(row(702, "rs1003", b(702), other(702), "AC=2;AF=0.2", ["0|1","0|0","1|0","0|0","0|0"]))

# 11) INDELS WITH OVERLAPPING WINDOWS: two insertions and a deletion at the SAME
#     POS (800) as separate records, a SNP between them (805), a 3-bp deletion at
#     810 that S1 carries on the SAME haplotype as the 800 insertion (S3: a
#     different sample), and an insertion at 820 on S1's OTHER haplotype. Each
#     contig is cut from the reference, so none holds a neighbouring indel.
lines.append(row(800, "rs1100", b(800), b(800)+"A", "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
lines.append(row(800, "rs1101", b(800), b(800)+"AAT", "AC=1;AF=0.1", ["0|0","1|0","0|0","0|0","0|0"]))
lines.append(row(800, "rs1102", ref[799:801], b(800), "AC=1;AF=0.1", ["0|0","0|0","0|0","0|1","0|0"]))
lines.append(row(805, "rs1103", b(805), other(805), "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
lines.append(row(810, "rs1104", ref[809:812], b(810), "AC=2;AF=0.2", ["1|0","0|0","0|1","0|0","0|0"]))
lines.append(row(820, "rs1105", b(820), b(820)+"GG", "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))

# 12) MULTIALLELIC records beyond cases 6-7: the indel ALT listed BEFORE the SNP
#     ALT (900); indel-only, two deletions + an insertion (960); REF > 1 with a
#     deletion and a same-length SNP-like ALT, which the enricher treats as an
#     indel contig and does not IUPAC-code (1020); two SNP ALTs + an insertion
#     with ONE AF value for three ALTs and two rsIDs (1080); 11 ALTs, so GT
#     indices reach 10 and 11 ("0|11" must not count as allele 1) (1140); a SNP
#     ALT next to a symbolic <DEL> (1200).
p=900
lines.append(row(p, "rs1200,rs1201", b(p), b(p)+"T,"+other(p), "AC=1,1;AF=0.1,0.2",
                 ["1|2","0|1","2|0","0|0","0|0"]))
p=960
r3=ref[p-1:p+2]
lines.append(row(p, "rs1202", r3, r3[0]+","+r3[:2]+","+r3+"T", "AC=1,1,1;AF=0.1,0.1,0.1",
                 ["1|2","0|3","2|3","1|1","0|0"]))
p=1020
r2=ref[p-1:p+1]
lines.append(row(p, "rs1203,rs1204", r2, r2[0]+","+other(p)+r2[1], "AC=1,1;AF=0.1,0.1",
                 ["1|2","0|2","2|2","0|0","0|0"]))
p=1080
lines.append(row(p, "rs1205,rs1206", b(p), other(p,0)+","+other(p,1)+","+b(p)+"CA", "AC=1,1,1;AF=0.25",
                 ["1|2","3|0","2|3","0|0","0|0"]))
p=1140
alts=[b(p)+"A"*k for k in range(1,11)]+[other(p)]
lines.append(row(p, "rs1207", b(p), ",".join(alts), "AC=1;AF="+",".join(["0.01"]*11),
                 ["0|11","10|0","1|11","2|10","0|0"]))
p=1200
lines.append(row(p, "rs1208", b(p), other(p)+",<DEL>", "AC=1,1;AF=0.1,0.1", ["0|1","2|0","1|2","0|0","0|0"]))

# 13) GENOTYPE SHAPES on SNP+indel combinations (cases 1-8 are all phased
#     diploid): unphased, missing, partially missing and haploid calls on an
#     insertion + a same-POS SNP + a window SNP (1260/1270); on a deletion + a
#     SNP in its deleted span (1300/1301); GT:DP FORMAT on an indel and a SNP
#     (1340/1350); multi-index calls (1/2, 2|., ./2, 2, 0/2) on a multiallelic
#     SNP+indel record (1400).
p=1260
lines.append(row(p, "rs1300", b(p), b(p)+"GT", "AC=3;AF=0.3", ["0/1","./.",".|1","1","1/1"]))
lines.append(row(p, "rs1301", b(p), other(p), "AC=2;AF=0.2", ["1/0",".|.","1|.","0","0|1"]))
lines.append(row(1270, "rs1302", b(1270), other(1270), "AC=2;AF=0.2", ["./1","1/1","0/0",".","1"]))
p=1300
lines.append(row(p, "rs1303", ref[p-1:p+2], b(p), "AC=3;AF=0.3", ["0/1","1/1","./.","1","0|1"]))
lines.append(row(1301, "rs1304", b(1301), other(1301), "AC=3;AF=0.3", ["1/0","./1","1","0|0","1|1"]))
p=1340
lines.append(row(p, "rs1305", b(p), b(p)+"A", "AC=3;AF=0.3", ["0|1:12","./.:0","1/1:30",".:5","1:7"],
                 fmt="GT:DP"))
lines.append(row(1350, "rs1306", b(1350), other(1350), "AC=3;AF=0.3", ["1|0:9","0|1:3","./.:0","1:4","0/1:8"],
                 fmt="GT:DP"))
p=1400
lines.append(row(p, "rs1307", b(p), other(p)+","+b(p)+"TT", "AC=1,4;AF=0.1,0.4", ["1/2","2|.","./2","2","0/2"]))

# 15) A REF THAT DOES NOT MATCH THE FASTA: the SNP IUPAC code is built from the
#     FASTA base, not the VCF REF (1460); an insertion whose REF base is not the
#     FASTA base at POS but occurs later in the window, so the first-occurrence
#     substitution lands there (1500); a 3-bp deletion whose REF occurs nowhere in
#     the window, so the contig stays reference (1540); lowercase REF/ALT on an
#     insertion (1580) and a lowercase ALT on a biallelic SNP (1585); a biallelic
#     SNP on a FASTA N (1615) next to an insertion whose window holds the N run
#     (1600); a biallelic SNP and an insertion on a FASTA IUPAC base 'R' (1650).
p=1460
lines.append(row(p, "rs1500", other(p,0), other(p,1), "AC=1;AF=0.1", ["0|1","1|0","0|0","0|0","0|0"]))
p=1500
wrong=[x for x in "ACGT" if x != b(p) and x in ref[p:p+27]][0]
lines.append(row(p, "rs1501", wrong, wrong+"TT", "AC=1;AF=0.1", ["1|0","0|0","0|0","0|0","0|0"]))
p=1540
win=ref[p-26:p+29]
bad=next(x+y+z for x in "ACGT" for y in "ACGT" for z in "ACGT" if x+y+z not in win)
lines.append(row(p, "rs1502", bad, bad[0], "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))
p=1580
lines.append(row(p, "rs1503", b(p).lower(), b(p).lower()+"at", "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
lines.append(row(1585, "rs1504", b(1585), other(1585).lower(), "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))
lines.append(row(1600, "rs1505", b(1600), b(1600)+"A", "AC=1;AF=0.1", ["1|0","0|0","0|0","0|0","0|0"]))
lines.append(row(1615, "rs1506", b(1615), "A", "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))
lines.append(row(1650, "rs1507", b(1650), "C", "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))
lines.append(row(1650, "rs1508", b(1650), b(1650)+"T", "AC=1;AF=0.1", ["0|0","1|0","0|0","0|0","0|0"]))

# 16) SAME-POS SNP and indel records beyond case 8: the indel record FIRST (1700);
#     a SNP + a 3-bp deletion (1730); two SNP records -- a split multiallelic SNP --
#     + an insertion (1760). Up to 2.9.1 both ports kept only the LAST SNP record's
#     ALT there (each record rebuilt the code from the reference base and overwrote
#     the dict key); both now code both ALTs and list both in the dict (cases 18-30).
p=1700
lines.append(row(p, "rs1600", b(p), b(p)+"TG", "AC=2;AF=0.2", ["0|1","1|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs1601", b(p), other(p), "AC=2;AF=0.2", ["1|0","0|0","0|1","0|0","0|0"]))
p=1730
lines.append(row(p, "rs1602", b(p), other(p), "AC=3;AF=0.3", ["0|1","0|0","1|1","0|0","0|0"]))
lines.append(row(p, "rs1603", ref[p-1:p+2], b(p), "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
p=1760
lines.append(row(p, "rs1604", b(p), other(p,0), "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs1604", b(p), other(p,1), "AC=1;AF=0.1", ["0|0","1|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs1605", b(p), b(p)+"CC", "AC=1;AF=0.1", ["0|0","0|0","1|0","0|0","0|0"]))

# 17) MULTIALLELIC SNP records whose single-character ALTs include a member that is
#     not A/C/G/T: the '*' of a GATK spanning deletion (1228, and with two real
#     ALTs at 1230), an 'N' (1232), a lowercase base (1234); a FASTA N carrying
#     two ALTs (1618, in the planted N run); and a split SNP record followed by a
#     multiallelic '*' record at the same POS (1229). The A/C/G/T subset is coded
#     and the other member adds nothing, so the real SNP is kept (enricher.py used
#     to leave the base unchanged here, enricher.cpp coded the subset).
p=1228
lines.append(row(p, "rs1700", b(p), other(p)+",*", "AC=2,1;AF=0.2,0.1", ["0|1","0|2","1|2","0|0","0|0"]))
p=1229
lines.append(row(p, "rs1701", b(p), other(p,0), "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs1702", b(p), other(p,1)+",*", "AC=1,1;AF=0.1,0.1", ["0|0","1|2","0|0","0|0","0|0"]))
p=1230
lines.append(row(p, "rs1703", b(p), other(p,0)+","+other(p,1)+",*", "AC=1,1,1;AF=0.1,0.1,0.1",
                 ["1|2","0|3","0|0","0|0","0|0"]))
p=1232
lines.append(row(p, "rs1704", b(p), other(p)+",N", "AC=1,1;AF=0.1,0.1", ["0|1","2|0","0|0","0|0","0|0"]))
p=1234
lines.append(row(p, "rs1705", b(p), other(p,0)+","+other(p,1).lower(), "AC=1,1;AF=0.1,0.1",
                 ["0|1","0|2","0|0","0|0","0|0"]))
p=1618
lines.append(row(p, "rs1706", b(p), "A,C", "AC=1,1;AF=0.1,0.1", ["0|1","0|2","0|0","0|0","0|0"]))

# ---- cases 18-30: multi-allelic sites SPLIT into one record per ALT --------------
# `bcftools norm -m -any` writes a site with several ALTs as one record per ALT, its
# genotypes recoded so 1 is that record's ALT. Each record must ADD its ALT to the
# position's code and dict entry (2.9.1 and earlier rebuilt the code from the
# reference base and overwrote the entry, so the last record won). A position whose
# reference base (A, C, G or T) and ALTs cover all four bases is Q, which the search
# reads as any base (N, the IUPAC letter for the four, is a gap to buildTST and
# searchTST).
import itertools
G1 = ["0|1","0|0","0|0","0|0","0|0"]

# 18) Two ALTs in both orders (1830, 1834); three ALTs -- all four bases -- in all six
#     orders (1840..1865), Q whatever the order and the AFs (the middle ALT is the
#     rarest); the same sites as ONE record (1870: two ALTs, 1875: three), whose codes
#     the split sites equal; three ALTs whose AFs tie (1880: Q too).
p=1830
lines.append(row(p, "rs1800", b(p), other(p,0), "AC=1;AF=0.1", G1))
lines.append(row(p, "rs1801", b(p), other(p,1), "AC=2;AF=0.2", ["0|0","1|0","0|1","0|0","0|0"]))
p=1834
lines.append(row(p, "rs1802", b(p), other(p,1), "AC=2;AF=0.2", ["0|0","1|0","0|1","0|0","0|0"]))
lines.append(row(p, "rs1803", b(p), other(p,0), "AC=1;AF=0.1", G1))
afs3 = ["0.3", "0.002", "0.05"]
gts3 = [["0|1","0|0","1|0","0|0","0|0"], ["0|0","1|0","0|0","0|0","0|0"], ["0|0","0|0","0|1","1|1","0|0"]]
for k, perm in enumerate(itertools.permutations(range(3))):
    p = 1840 + 5*k
    for j in perm:
        lines.append(row(p, "rs%d" % (1810 + 3*k + j), b(p), other(p,j), "AC=1;AF=" + afs3[j], gts3[j]))
p=1870
lines.append(row(p, "rs1830", b(p), other(p,0)+","+other(p,1), "AC=1,2;AF=0.1,0.2",
                 ["0|1","2|0","0|2","0|0","0|0"]))
p=1875
lines.append(row(p, "rs1831", b(p), ",".join(other(p,j) for j in range(3)), "AC=1,1,3;AF=0.3,0.002,0.05",
                 ["0|1","2|0","1|3","3|3","0|0"]))
p=1880
for j in (2, 0, 1):
    lines.append(row(p, "rs%d" % (1832 + j), b(p), other(p,j), "AC=1;AF=0.01", G1))

# 19) SPLIT RECORDS NOT ADJACENT: interleaved with another position's records and with
#     an insertion at the same POS (1900, 1905).
lines.append(row(1900, "rs1840", b(1900), other(1900,0), "AC=1;AF=0.1", G1))
lines.append(row(1905, "rs1841", b(1905), other(1905,0), "AC=1;AF=0.1", ["1|0","0|0","0|0","0|0","0|0"]))
lines.append(row(1900, "rs1842", b(1900), b(1900)+"TT", "AC=1;AF=0.1", ["0|0","0|0","1|0","0|0","0|0"]))
lines.append(row(1905, "rs1843", b(1905), other(1905,2), "AC=1;AF=0.1", ["0|0","0|1","0|0","0|0","0|0"]))
lines.append(row(1900, "rs1844", b(1900), other(1900,1), "AC=1;AF=0.2", ["0|0","1|0","0|0","0|0","0|0"]))

# 20) GENOTYPE SHAPES on split records: haploid, unphased, missing, partially missing
#     (1930: two ALTs), and GT:DP on three ALTs (1935).
lines.append(row(1930, "rs1850", b(1930), other(1930,0), "AC=2;AF=0.2", ["0/1","./.",".|1","1","0"]))
lines.append(row(1930, "rs1851", b(1930), other(1930,1), "AC=2;AF=0.2", ["1/0","1|.","0|0",".","1"]))
for j, gts in enumerate((["0|1:9","./.:0","1/1:30",".:5","1:7"], ["1|0:3","0|0:2","./1:4","0:6","0:1"],
                         ["0|0:1",".|1:2","0|0:3","1:8","0/1:5"])):
    lines.append(row(1935, "rs%d" % (1852 + j), b(1935), other(1935,j), "AC=2;AF=0.%d" % (j + 1), gts,
                     fmt="GT:DP"))

# 21) SNP + INDEL AT ONE POSITION with split SNPs: three SNP ALTs, an insertion and a
#     4-bp deletion at 1980; split SNP sites inside those indels' windows (1990: two
#     ALTs, 2005: three). The contigs are cut from the reference, so they stay those
#     of the indels alone.
p=1980
lines.append(row(p, "rs1860", b(p), other(p,0), "AC=1;AF=0.1", G1))
lines.append(row(p, "rs1861", b(p), b(p)+"GA", "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
lines.append(row(p, "rs1862", b(p), other(p,1), "AC=1;AF=0.3", ["0|0","1|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs1863", ref[p-1:p+3], b(p), "AC=1;AF=0.2", ["0|0","0|0","1|0","0|0","0|0"]))
lines.append(row(p, "rs1864", b(p), other(p,2), "AC=1;AF=0.05", ["0|0","0|0","0|1","0|0","0|0"]))
for j in (1, 0):
    lines.append(row(1990, "rs%d" % (1865 + j), b(1990), other(1990,j), "AC=1;AF=0.1",
                     [["1|0","0|0","0|0","0|0","0|0"], ["0|0","0|1","0|0","0|0","0|0"]][j]))
for j in (0, 2, 1):
    lines.append(row(2005, "rs%d" % (1867 + j), b(2005), other(2005,j), "AC=1;AF=0.%d" % (j + 1),
                     [["1|0","0|0","0|0","0|0","0|0"], ["0|0","1|0","0|0","0|0","0|0"],
                      ["0|0","0|0","1|1","0|0","0|0"]][j]))

# 22) A RECORD LISTING SEVERAL ALTS and a split record at ONE position, both orders
#     (2050: the multi-ALT record first; 2055: the split record first).
p=2050
lines.append(row(p, "rs1870", b(p), other(p,0)+","+other(p,1), "AC=1,1;AF=0.1,0.2", ["0|1","2|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs1871", b(p), other(p,2), "AC=1;AF=0.05", ["0|0","0|0","1|0","0|0","0|0"]))
p=2055
lines.append(row(p, "rs1872", b(p), other(p,2), "AC=1;AF=0.3", ["0|0","0|0","1|0","0|0","0|0"]))
lines.append(row(p, "rs1873", b(p), other(p,0)+","+other(p,1), "AC=1,1;AF=0.1,0.002", ["0|1","2|0","0|0","0|0","0|0"]))

# 23) A NON-PASS split record at a multi-allelic position adds nothing (2100).
p=2100
lines.append(row(p, "rs1880", b(p), other(p,0), "AC=1;AF=0.1", G1))
lines.append(row(p, "rs1881", b(p), other(p,1), "AC=5;AF=0.5", ["1|1","1|1","1|0","0|0","0|0"], filt="LowQual"))
lines.append(row(p, "rs1882", b(p), other(p,2), "AC=1;AF=0.2", ["0|0","0|1","0|0","0|0","0|0"]))

# 24) THE SAME ALT ON TWO RECORDS: one entry with the carriers of both, and the rsID and
#     AF of the more frequent record (2130; up to 2.9.1 the later record's entry
#     replaced the earlier).
p=2130
lines.append(row(p, "rs1890", b(p), other(p,0), "AC=1;AF=0.1", G1))
lines.append(row(p, "rs1891", b(p), other(p,1), "AC=1;AF=0.2", ["0|0","0|1","0|0","0|0","0|0"]))
lines.append(row(p, "rs1892", b(p), other(p,0), "AC=1;AF=0.3", ["0|0","0|0","1|0","0|0","0|0"]))

# 25) AF SPELLINGS at four-base sites: an exponent, '.', no AF value (2160); an ALT on
#     two records, whose one entry takes the rsID and AF of the more frequent (2170).
#     Both are Q whatever the AFs (the four-base rule ranked the ALTs by them).
p=2160
lines.append(row(p, "rs1900", b(p), other(p,0), "AC=1;AF=2e-1", G1))
lines.append(row(p, "rs1901", b(p), other(p,1), "AC=1;AF=.", ["0|0","0|1","0|0","0|0","0|0"]))
lines.append(row(p, "rs1902", b(p), other(p,2), "AC=1", ["0|0","0|0","0|1","0|0","0|0"]))
p=2170
lines.append(row(p, "rs1903", b(p), other(p,0), "AC=1;AF=0.001", G1))
lines.append(row(p, "rs1904", b(p), other(p,1), "AC=1;AF=0.01", ["0|0","0|1","0|0","0|0","0|0"]))
lines.append(row(p, "rs1905", b(p), other(p,2), "AC=1;AF=0.1", ["0|0","0|0","0|1","0|0","0|0"]))
lines.append(row(p, "rs1906", b(p), other(p,0), "AC=1;AF=0.5", ["0|0","0|0","0|0","0|1","0|0"]))

# 26) A REFERENCE IUPAC base with two split ALTs (2200, planted 'Y'): it adds no base
#     of its own, as for one record; a reference N with four split ALTs (1614, in the
#     planted N run): it stays N, as in case 33.
lines.append(row(2200, "rs1910", "Y", "A", "AC=1;AF=0.1", G1))
lines.append(row(2200, "rs1911", "Y", "G", "AC=1;AF=0.2", ["0|0","0|1","0|0","0|0","0|0"]))
for j, (alt, af) in enumerate((("T", "0.2"), ("A", "0.4"), ("G", "0.01"), ("C", "0.3"))):
    lines.append(row(1614, "rs%d" % (1912 + j), "N", alt, "AC=1;AF=" + af, G1))

# 27) A '*' split record after a SNP record (2230): listed in the dict, not coded.
p=2230
lines.append(row(p, "rs1920", b(p), other(p,0), "AC=1;AF=0.1", G1))
lines.append(row(p, "rs1921", b(p), "*", "AC=1;AF=0.2", ["0|0","1|0","0|0","0|0","0|0"]))

# 29) FOUR-BASE sites whose records carry AF_<x> fields: no AF, only per-dataset AF_<x>
#     and AF_max as a sites-only merge (CRISPRme's mega) writes them (2260), and an
#     AF_pop above every AF at its position (2270). Both are Q whatever these fields
#     say: the four-base rule ranked the ALTs by AF and every AF_<x>, 2.10.0 ranks
#     nothing.
p=2260
for j, info in enumerate(("AF_one=0.01;AF_max=0.01", "AF_two=0.0002;AF_max=0.0002",
                          "AF_one=0.2;AF_two=0.3;AF_max=0.3")):
    lines.append(row(p, "rs%d" % (1930 + j), b(p), other(p,j), info, G1))
p=2270
for j, info in enumerate(("AC=1;AF=0.001;AF_pop=0.4", "AC=1;AF=0.01", "AC=1;AF=0.002")):
    lines.append(row(p, "rs%d" % (1933 + j), b(p), other(p,j), info, G1))

# 30) SPLIT RECORDS PHASED INDEPENDENTLY, as HGDP + 1000 Genomes is: one haplotype slot
#     holds two ALTs of one position -- S1 0|1 and 0|1 (chr22:10661629, NA18869), S2 1|0
#     and 1|0, S3 more ALT copies than slots (1|1 and 0|1), S4 one ALT per slot (2300);
#     three ALTs on one slot at a four-base position (2305). Every entry keeps its own
#     record's genotypes verbatim, in both ports: resolving the conflict is CRISPRme's
#     (docs/MULTIALLELIC_ENRICHMENT.md, section 6.3).
p=2300
lines.append(row(p, "rs1940", b(p), other(p,0), "AC=5;AF=0.5", ["0|1","1|0","1|1","1|0","0|0"]))
lines.append(row(p, "rs1941", b(p), other(p,1), "AC=4;AF=0.4", ["0|1","1|0","0|1","0|1","0|0"]))
p=2305
for j, af in enumerate(("0.3", "0.001", "0.2")):
    lines.append(row(p, "rs%d" % (1942 + j), b(p), other(p,j), "AC=1;AF=" + af, G1))

# 31) INFO/AF READ BY ITS KEY in every record: the file's first record (case 1) has AF
#     second in INFO, as up to 2.9.1 every record's AF was looked for there; here AF
#     stands first (2330), after a gnomAD-style AF_afr (2335), third on an insertion
#     (2340, the indel log's AF) and nowhere, beside AF_max (2345: no AF).
lines.append(row(2330, "rs1950", b(2330), other(2330), "AF=0.3;AC=1", G1))
lines.append(row(2335, "rs1951", b(2335), other(2335), "DP=7;AC=1;AF_afr=0.5;AF=0.15", G1))
lines.append(row(2340, "rs1952", b(2340), b(2340)+"GT", "AN=10;AC=1;AF=0.05", G1))
lines.append(row(2345, "rs1953", b(2345), other(2345), "AF_max=0.2;AC=1", G1))

# 32) THE SAME REF>ALT ON TWO SPLIT RECORDS with carriers on both, one of them on the
#     same slot and one on the other (2360), and on a record listing several ALTs then
#     a split record (2365): one entry each, genotypes in split form (1 = that ALT).
p=2360
lines.append(row(p, "rs1960", b(p), other(p,0), "AC=2;AF=0.1", ["0|1","0|0","0|1","0|0","0|0"]))
lines.append(row(p, "rs1961", b(p), other(p,1), "AC=1;AF=0.2", ["0|0","1|0","0|0","0|0","0|0"]))
lines.append(row(p, "rs1962", b(p), other(p,0), "AC=3;AF=0.4", ["1|0","0|0","0|1","1|1","0|0"]))
p=2365
lines.append(row(p, "rs1963", b(p), other(p,0)+","+other(p,1), "AC=1,2;AF=0.1,0.2",
                 ["1|2","0|2","0|0","0|0","0|0"]))
lines.append(row(p, "rs1964", b(p), other(p,1), "AC=1;AF=0.05", ["0|0","0|0","1|0","0|0","0|0"]))

# 33) ALTS ON A REFERENCE N (the planted N run 1610-1620): three split ALTs (1612) and
#     a record listing all four (1616) leave the N; and on reference bytes that are no
#     base or IUPAC code, which the search reads as N: two split ALTs (the X at 2420)
#     and a record listing two (the '-' at 2425) leave them; the four-base sites
#     elsewhere are Q, here one record listing three SNP ALTs (2380).
for j, alt in enumerate("ACG"):
    lines.append(row(1612, "rs%d" % (1970 + j), "N", alt, "AC=1;AF=0.1", G1))
lines.append(row(1616, "rs1973", "N", "A,C,G,T", "AC=1,1,1,1;AF=0.1,0.1,0.1,0.1",
                 ["0|1","2|3","0|4","0|0","0|0"]))
lines.append(row(2420, "rs1975", "A", "C", "AC=1;AF=0.1", G1))
lines.append(row(2420, "rs1976", "A", "G", "AC=1;AF=0.2", ["0|0","1|0","0|0","0|0","0|0"]))
lines.append(row(2425, "rs1977", "A", "C,T", "AC=1,1;AF=0.1,0.2", ["0|1","2|0","0|0","0|0","0|0"]))
p=2380
lines.append(row(p, "rs1974", b(p), ",".join(other(p,j) for j in range(3)), "AC=1,1,1;AF=0.1,0.2,0.3",
                 ["0|1","0|2","3|0","0|0","0|0"]))

# 14) THE CHROMOSOME EDGES (records listed last; the enricher does not need
#     sorted input): an insertion at POS 1 and a 2-bp deletion at POS 25 have
#     POS-26 < 0, so BOTH ports must skip them with the same stdout warning; an
#     insertion at POS 26 (window start exactly 0) with SNPs at 26 and at 1; at
#     the END, a deletion whose REF runs past the last base (cannot match), a
#     deletion of the last two bases, and a SNP + an insertion on the last base
#     (windows truncated by the chromosome end). 28) A second SNP record at POS 1
#     and at the last base makes each a split multi-allelic site.
L=len(ref)
lines.append(row(1, "rs1400", b(1), other(1), "AC=1;AF=0.1", ["0|1","0|0","0|0","0|0","0|0"]))
lines.append(row(1, "rs1401", b(1), b(1)+"A", "AC=2;AF=0.2", ["1|0","0|1","0|0","0|0","0|0"]))
lines.append(row(25, "rs1402", ref[24:26], b(25), "AC=1;AF=0.1", ["1|0","0|0","0|0","0|0","0|0"]))
lines.append(row(26, "rs1403", b(26), b(26)+"CC", "AC=2;AF=0.2", ["0|1","1|0","0|0","0|0","0|0"]))
lines.append(row(26, "rs1404", b(26), other(26), "AC=2;AF=0.2", ["1|0","0|0","0|1","0|0","0|0"]))
lines.append(row(L-2, "rs1405", ref[L-3:L]+"AA", b(L-2), "AC=1;AF=0.1", ["0|0","0|0","0|0","1|0","0|0"]))
lines.append(row(L-1, "rs1406", ref[L-2:L], b(L-1), "AC=1;AF=0.1", ["0|0","0|0","0|0","0|0","1|0"]))
lines.append(row(L, "rs1407", b(L), other(L), "AC=1;AF=0.1", ["0|0","0|0","1|0","0|0","0|0"]))
lines.append(row(L, "rs1408", b(L), b(L)+"GGG", "AC=2;AF=0.2", ["0|0","0|0","0|0","1|1","0|0"]))
lines.append(row(1, "rs1409", b(1), other(1,1), "AC=1;AF=0.3", ["0|0","1|0","0|0","0|0","0|0"]))
lines.append(row(L, "rs1410", b(L), other(L,1), "AC=1;AF=0.3", ["0|0","0|0","0|1","0|0","0|0"]))

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

# every file either port wrote, not only the four above: an extra or a missing
# output must fail too
( cd "$PYDIR"  && find . -type f ! -name run.log | sort ) > "$WORK/files_py.txt"
( cd "$CPPDIR" && find . -type f ! -name run.log | sort ) > "$WORK/files_cpp.txt"
if diff "$WORK/files_py.txt" "$WORK/files_cpp.txt" > /dev/null 2>&1; then
    echo "  OK   : same file set ($(wc -l < "$WORK/files_py.txt" | tr -d ' ') files)"
else
    echo "  DIFF : file set"
    diff "$WORK/files_py.txt" "$WORK/files_cpp.txt"
    FAIL=1
fi

# stdout+stderr: the skip warnings of case 14 must match as well; only the
# elapsed seconds on the DONE line may differ
sed 's/^DONE IN .*/DONE IN/' "$PYDIR/run.log"  > "$WORK/log_py.txt"
sed 's/^DONE IN .*/DONE IN/' "$CPPDIR/run.log" > "$WORK/log_cpp.txt"
if diff "$WORK/log_py.txt" "$WORK/log_cpp.txt" > /dev/null 2>&1; then
    echo "  OK   : run log (stdout+stderr; $(grep -c WARNING "$WORK/log_py.txt") skip warning(s))"
else
    echo "  DIFF : run log"
    diff "$WORK/log_py.txt" "$WORK/log_cpp.txt" | head -20
    FAIL=1
fi

# ---------------------------------------------------------------------------
# agreement is not correctness: the outputs against what they must be
# ---------------------------------------------------------------------------
echo "Checking outputs against an independent reconstruction:"
if python3 - "$CHR_FA" "$VCF_PLAIN" "$PYDIR" "$CPPDIR" "$ENRICHED" "$JSON" "$LOG" <<'PY'
import json, sys
fa, vcf, pydir, cppdir, enriched, jsn, log = sys.argv[1:]
ref = "".join(open(fa).read().splitlines()[1:]).upper()
BIT = {"A": 1, "C": 2, "G": 4, "T": 8}
CODE = {3: "M", 5: "R", 9: "W", 6: "S", 10: "Y", 12: "K", 14: "B", 13: "D", 11: "H", 7: "V", 15: "Q"}
# every PASS record with a one-base REF codes its one-character A/C/G/T ALTs at POS,
# added to the reference base when that is A/C/G/T; a reference N, or a byte that is no
# base or IUPAC code, stays as it is; one base leaves the reference letter
mask = {}
for line in open(vcf):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[6] not in ("PASS", ".") or len(f[3]) != 1:
        continue
    pos = int(f[1]) - 1
    m = mask.get(pos, BIT.get(ref[pos], 0))
    for alt in f[4].split(","):
        if len(alt) == 1:
            m |= BIT.get(alt, 0)
    mask[pos] = m
want = list(ref)
for pos, m in mask.items():
    if ref[pos] in "ACGTRYSWKMBDHV" and m in CODE:
        want[pos] = CODE[m]
want = "".join(want)
bad = 0
for port, d in (("enricher.py", pydir), ("enricher.cpp", cppdir)):
    got = "".join(open("%s/%s" % (d, enriched)).read().splitlines()[1:])
    diff = [i + 1 for i in range(max(len(got), len(want))) if got[i:i + 1] != want[i:i + 1]]
    nref = [i for i, c in enumerate(ref) if c == "N"]
    ngot = [i for i, c in enumerate(got) if c == "N"]
    print("  %s  : %s enriched genome == the reconstruction (%d bases, %d coded, %d Q)%s"
          % ("OK  " if not diff else "DIFF", port, len(want), len(mask),
             want.count("Q"), "" if not diff else ": differs at %s" % diff[:10]))
    print("  %s  : %s N-count invariant: %d N in the reference, %d in the enriched genome%s"
          % ("OK  " if nref == ngot else "DIFF", port, len(nref), len(ngot),
             ", at the same positions" if nref == ngot else ""))
    bad += bool(diff) + (nref != ngot)
assert want.count("Q") >= 12 and "N" not in [want[p] for p in mask if ref[p] != "N"]
assert [want[2419], want[2424]] == ["X", "-"] and 2419 in mask and 2424 in mask
# the dict at the cases that pin merged alleles, per-record AF, the records at a
# reference N, X or '-' and a record listing three SNP ALTs; the indel log's AF read by
# key
d = json.load(open("%s/%s" % (pydir, jsn)))
def entries(p):
    v = d.get("chrTest,%d" % p)
    return [] if v is None else v.split("$")
def other(p, k):
    return [x for x in "ACGT" if x != ref[p - 1]][k]
r = lambda p: ref[p - 1]
expect = {
    1612: ("listed at a reference N",
           ["S1:0|1;N,A;rs1970;0.1", "S1:0|1;N,C;rs1971;0.1", "S1:0|1;N,G;rs1972;0.1"]),
    1616: ("listed at a reference N",
           ["S1:0|1;N,A;rs1973;0.1", "S2:2|3;N,C;rs1973;0.1", "S2:2|3;N,G;rs1973;0.1",
            "S3:0|4;N,T;rs1973;0.1"]),
    2420: ("listed at a reference X", ["S1:0|1;A,C;rs1975;0.1", "S2:1|0;A,G;rs1976;0.2"]),
    2425: ("listed at a reference '-'", ["S1:0|1;A,C;rs1977;0.1", "S2:2|0;A,T;rs1977;0.2"]),
    2130: ("carriers merged",
           ["S1:0|1,S3:1|0;%s,%s;rs1892;0.3" % (r(2130), other(2130, 0)),
            "S2:0|1;%s,%s;rs1891;0.2" % (r(2130), other(2130, 1))]),
    2330: ("AF by key", ["S1:0|1;%s,%s;rs1950;0.3" % (r(2330), other(2330, 0))]),
    2335: ("AF by key", ["S1:0|1;%s,%s;rs1951;0.15" % (r(2335), other(2335, 0))]),
    2345: ("AF by key", ["S1:0|1;%s,%s;rs1953;" % (r(2345), other(2345, 0))]),
    2360: ("carriers merged",
           ["S1:1|1,S3:0|1,S4:1|1;%s,%s;rs1962;0.4" % (r(2360), other(2360, 0)),
            "S2:1|0;%s,%s;rs1961;0.2" % (r(2360), other(2360, 1))]),
    2365: ("carriers merged",
           ["S1:1|2;%s,%s;rs1963;0.1" % (r(2365), other(2365, 0)),
            "S1:0|1,S2:0|1,S3:1|0;%s,%s;rs1963;0.2" % (r(2365), other(2365, 1))]),
    2380: ("a record listing three SNP ALTs",
           ["S1:0|1;%s,%s;rs1974;0.1" % (r(2380), other(2380, 0)),
            "S2:0|2;%s,%s;rs1974;0.2" % (r(2380), other(2380, 1)),
            "S3:3|0;%s,%s;rs1974;0.3" % (r(2380), other(2380, 2))]),
}
for p, (what, e) in sorted(expect.items()):
    ok = entries(p) == e
    bad += not ok
    print("  %s  : dict at %d (%s)" % ("OK  " if ok else "DIFF", p,
                                       what if ok else "$".join(entries(p)) or "no such key"))
rows = [x.split("\t") for x in open("%s/%s" % (pydir, log)).read().splitlines()[1:]]
af = [x[3] for x in rows if x[4].startswith("chrTest_2340_")]
bad += af != ["0.05"]
print("  %s  : indel log AF at 2340 by key (%s)" % ("OK  " if af == ["0.05"] else "DIFF", af))
sys.exit(1 if bad else 0)
PY
then :; else FAIL=1; fi

# a reference holding Q, the four-allele symbol, is refused by both ports alike
QDIR="$WORK/refQ"
mkdir -p "$QDIR/py/SNPs_genome" "$QDIR/cpp/SNPs_genome"
sed '3s/^\(.\{10\}\)./\1q/' "$CHR_FA" > "$QDIR/chrTest.fa"
for port in py cpp; do
    if [ "$port" = py ]; then cmd=(python3 "$ENRICHER_PY"); else cmd=("$ENRICHER_BIN"); fi
    ( cd "$QDIR/$port" && "${cmd[@]}" "$VCF_GZ" "$QDIR/chrTest.fa" "resultTest" "yes" "$VCF_PLAIN" ) \
        > "$QDIR/$port.out" 2> "$QDIR/$port.err"
    echo $? > "$QDIR/$port.rc"
done
QMSG="ERROR: chrTest holds 'Q' at 71: a reference FASTA must not hold Q, which is reserved for four-allele SNP positions"
if [ "$(cat "$QDIR/py.rc")" = 1 ] && [ "$(cat "$QDIR/cpp.rc")" = 1 ] \
   && [ "$(cat "$QDIR/py.err")" = "$QMSG" ] && [ "$(cat "$QDIR/cpp.err")" = "$QMSG" ]; then
    echo "  OK   : a reference holding 'q' is refused by both ports (exit 1, the same message)"
else
    echo "  DIFF : a reference holding 'q': exit $(cat "$QDIR/py.rc") / $(cat "$QDIR/cpp.rc")"
    cat "$QDIR/py.err" "$QDIR/cpp.err"
    FAIL=1
fi

echo
if [ "$FAIL" -eq 0 ]; then
    echo "BYTE-IDENTICAL: PASS"
    exit 0
else
    echo "BYTE-IDENTICAL: FAIL"
    exit 1
fi
