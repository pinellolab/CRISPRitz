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
# or whose FASTA base is N. Besides the four outputs, the file set and the run log
# (stdout+stderr, e.g. the skip warnings) must match.
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
# synthetic chromosome FASTA (1,800 bp, mixed case to exercise upper()). The
# first 400 bp -- all that cases 1-8 touch -- are the same draws as when the
# chromosome was 400 bp long; the rest holds cases 9-16.
# ---------------------------------------------------------------------------
CHR_FA="$WORK/chrTest.fa"
python3 - "$CHR_FA" <<'PY'
import sys, random
random.seed(1234)
seq = [random.choice('ACGTacgt') for _ in range(1800)]
# case 15: a hard-masked N run and an IUPAC base in the reference itself
for i in range(1610, 1621):
    seq[i - 1] = 'N'
seq[1650 - 1] = 'R'
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
#     + an insertion (1760). For the latter BOTH ports keep only the LAST SNP
#     record's ALT in the enriched genome and the dict (each record rebuilds the
#     code from the reference base and overwrites the dict key); that is
#     enricher.py's long-standing behaviour, pinned here only for port agreement.
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

# 14) THE CHROMOSOME EDGES (records listed last; the enricher does not need
#     sorted input): an insertion at POS 1 and a 2-bp deletion at POS 25 have
#     POS-26 < 0, so BOTH ports must skip them with the same stdout warning; an
#     insertion at POS 26 (window start exactly 0) with SNPs at 26 and at 1; at
#     the END, a deletion whose REF runs past the last base (cannot match), a
#     deletion of the last two bases, and a SNP + an insertion on the last base
#     (windows truncated by the chromosome end).
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

echo
if [ "$FAIL" -eq 0 ]; then
    echo "BYTE-IDENTICAL: PASS"
    exit 0
else
    echo "BYTE-IDENTICAL: FAIL"
    exit 1
fi
