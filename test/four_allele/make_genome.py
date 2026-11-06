#!/usr/bin/env python3
"""Synthetic genome and VCF for the four-allele test (test/run_four_allele_test.sh).

Writes a reference FASTA per contig and a VCF per contig whose SNP records give some
positions all four bases between the reference base and the ALTs -- split into one
record per ALT, as `bcftools norm -m -any` writes a panel, or one record listing three
ALTs. `crispritz.py add-variants` codes each such position Q, the symbol buildTST,
searchTST and searchBruteForce read as any base (N stays a gap, which matches nothing).
Two or three bases give the usual IUPAC codes, and ALTs at a reference N leave the N.

The contigs plant the test guide G next to an NGG PAM (3') and after a TTTV PAM (5'),
written in target-strand orientation:

  Qg   a guide base whose reference base mismatches G, made a match by Q (at several
       distances from the PAM, one and six Q's in one guide);
  Qp   a PAM that exists only through Q: NGG with its GG (or one G) on Q, NGG's N on Q,
       TTTV with a T or the V on Q;
  Qd   Q beyond the guide's far end (the pad a DNA-bulged alignment reads), and Q there
       right before an N run;
  Qe   Q on a contig's first and last bases, inside a site at the contig end;
  Qr   a run of 30 Q's (every window there matches any guide);
  Qn   ALTs at reference N positions inside an N run beside a site (the N stays);
  bg   a 20-kb background: random bases, 3% four-allele positions, 4% two- or
       three-base ones, N runs, and G planted with 0-4 mismatches next to NGG and
       after TTTV;
  ctl  a contig with no variant (its index files are not flagged for Q);
  Ng   G with one base an N, next to NGG: a guide region holding exactly one N, which
       no index stores and the brute force must drop whatever the mismatch budget.

Each case is written on the + strand and, reverse-complemented with its VCF, on the -
strand. sites.tsv lists every planted site whose enriched sequence is a perfect match
of G with its PAM (3' or 5'), for the test to check each is found.

usage: make_genome.py <out_dir> [seed]
writes <out_dir>/ref/<contig>.fa, <out_dir>/vcf/fourallele.<contig>.vcf.gz and
<out_dir>/sites.tsv
"""
import gzip
import os
import random
import sys

G = "GACTTAGCCTAAGGTCTACG"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
OTHER = {"A": "C", "C": "A", "G": "T", "T": "G"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def main():
    out = sys.argv[1]
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 2027)

    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    cases = []  # (name, sequence, {0-based pos: ALT bases}, [(site start, side)])

    def case(name, parts):
        """parts: (text, alts) pairs. text is reference bases; alts maps an offset in
        text to the ALT bases a record set adds there ("*": all but the reference
        base, i.e. four bases in all; or an explicit string of bases). A site is marked
        by a part ("SITE3", text) or ("SITE5", text): its text starts at the site."""
        seq, alts, sites = "", {}, []
        for part in parts:
            if part[0] in ("SITE3", "SITE5"):
                sites.append((len(seq), part[0][-1]))
                part = part[1:]
            text, a = part
            for off, bases in a.items():
                ref = text[off]
                alts[len(seq) + off] = "".join(b for b in "ACGT" if b != ref) if bases == "*" else bases
            seq += text
        cases.append((name, seq, alts, sites))

    pam3 = "TGG"
    # Qg: the guide base at distance k from the far (5') end mismatches, Q makes it match
    for k in (0, 4, 9, 14, 19):
        g = list(G)
        g[k] = OTHER[g[k]]
        case("Qg%d" % k, [(rnd(30), {}), ("SITE3", "".join(g) + pam3, {k: "*"}), (rnd(30), {})])
    g = list(G)
    for k in (1, 3, 6, 11, 15, 18):
        g[k] = OTHER[g[k]] if k % 2 else g[k]
    case("Qg6", [(rnd(30), {}), ("SITE3", "".join(g) + pam3, {k: "*" for k in (1, 3, 6, 11, 15, 18)}),
                 (rnd(30), {})])
    # Qp: the PAM only through Q
    for label, pam, qs in (("GG", "TCA", (1, 2)), ("G2", "TGA", (2,)), ("G1", "TAG", (1,)), ("N", "AGG", (0,))):
        case("Qp" + label, [(rnd(30), {}), ("SITE3", G + pam, {20 + q: "*" for q in qs}), (rnd(30), {})])
    for label, pam, qs in (("T", "TCTA", (1,)), ("V", "TTTT", (3,)), ("TT", "AATC", (0, 1))):
        case("Qp5" + label, [(rnd(30), {}), ("SITE5", pam + G, {q: "*" for q in qs}), (rnd(30), {})])
    # Qd: Q beyond the guide's far end, and Q there before an N run
    case("Qd", [(rnd(30), {}), (rnd(3), {0: "*", 2: "*"}), ("SITE3", G + pam3, {}), (rnd(30), {})])
    case("QdN", [(rnd(30), {}), ("NNNNNNNN", {}), ("A", {0: "*"}), ("SITE3", G + pam3, {}), (rnd(30), {})])
    case("Qd5", [(rnd(30), {}), ("SITE5", "TTTA" + G, {}), (rnd(2), {1: "*"}), ("NNNNNN", {}), (rnd(30), {})])
    # Qe: Q on the contig's first and last bases, in a site at each end
    case("Qe", [("SITE3", G + "TGG", {0: "*", 22: "*"})])
    case("Qe5", [("SITE5", "TTTC" + G, {0: "*", 23: "*"})])
    # Qr: a run of 30 Q's
    case("Qr", [(rnd(40), {}), (rnd(30), {k: "*" for k in range(30)}), (rnd(40), {})])
    # Qn: ALTs on reference N beside a site: the N stays, the site's pad stops there
    case("Qn", [(rnd(30), {}), ("NNNNNNNNNN", {2: "AC", 4: "ACG", 6: "ACGT"}), ("SITE3", G + pam3, {}),
                (rnd(30), {})])
    # bg: a background with sprinkled codes, N runs and planted guides
    bg, alts = list(rnd(20000)), {}
    for start in (3000, 11000, 17000):
        for i in range(start, start + 40 + rng.randrange(60)):
            bg[i] = "N"
    sites = []
    for i in range(0, 19900, 400):
        if "N" in bg[max(i - 1, 0):i + 30]:
            continue
        g = list(G)
        for k in rng.sample(range(20), rng.randrange(5)):
            g[k] = OTHER[g[k]]
        if (i // 400) % 2:
            bg[i:i + 23] = list("".join(g) + "TGG")
        else:
            bg[i:i + 24] = list("TTTA" + "".join(g))
    for i in range(20000):
        if bg[i] == "N":
            if rng.random() < 0.05:
                alts[i] = "".join(rng.sample("ACGT", rng.randrange(1, 5)))
            continue
        r = rng.random()
        if r < 0.03:
            alts[i] = "".join(b for b in "ACGT" if b != bg[i])
        elif r < 0.07:
            alts[i] = "".join(rng.sample([b for b in "ACGT" if b != bg[i]], rng.randrange(1, 3)))
    cases.append(("bg", "".join(bg), alts, []))
    case("ctl", [(rnd(40), {}), ("SITE3", G + pam3, {}), (rnd(40), {})])
    # Ng: one N in a guide region beside a PAM. searchBruteForce counted an N as 101
    # mismatches, so a budget of 101 or more reported the window; the N runs above give
    # none (no PAM where a guide reaches exactly one base into a run)
    case("Ng", [(rnd(30), {}), (G[:7] + "N" + G[8:] + pam3, {}), (rnd(30), {})])

    os.makedirs(os.path.join(out, "ref"), exist_ok=True)
    os.makedirs(os.path.join(out, "vcf"), exist_ok=True)
    samples = ["S1", "S2", "S3"]
    n_rec = n_quad = 0
    with open(os.path.join(out, "sites.tsv"), "w") as tab:
        for name, seq, alts, sites in cases:
            for strand in "+-":
                chrom = "chr4a%s%s" % (name, "p" if strand == "+" else "m")
                L = len(seq)
                s = seq if strand == "+" else revcomp(seq)
                a = alts if strand == "+" else {L - 1 - p: "".join(COMP[b] for b in bases)
                                                 for p, bases in alts.items()}
                with open(os.path.join(out, "ref", chrom + ".fa"), "w") as fh:
                    fh.write(">%s\n" % chrom)
                    for i in range(0, L, 60):
                        fh.write(s[i:i + 60] + "\n")
                recs = []
                for p in sorted(a):
                    bases = [b for b in a[p] if b != s[p]]
                    ref = s[p]
                    if ref != "N" and len(set(bases) | {ref}) == 4:
                        n_quad += 1
                    if p % 7 == 3 and len(bases) > 1:  # one record listing every ALT
                        gts = ["0|1", "%d|0" % len(bases), "0|0"]
                        recs.append((p, ref, ",".join(bases), ",".join(["1"] * len(bases)), gts))
                    else:  # split, as bcftools norm -m -any writes it
                        for k, b in enumerate(bases):
                            gts = [["0|1", "0|0", "0|0"], ["0|0", "1|0", "0|0"], ["0|0", "0|0", "1|1"]][k % 3]
                            recs.append((p, ref, b, "1", gts))
                if not recs:
                    continue
                vcf = os.path.join(out, "vcf", "fourallele.%s.vcf.gz" % chrom)
                with gzip.open(vcf, "wt") as fh:
                    fh.write("##fileformat=VCFv4.2\n##contig=<ID=%s,length=%d>\n" % (chrom, L))
                    fh.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n')
                    fh.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
                    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(samples))
                    for p, ref, alt, ac, gts in recs:
                        fh.write("\t".join([chrom, str(p + 1), "v%d" % (p + 1), ref, alt, ".", "PASS",
                                            "AC=%s;AF=%s" % (ac, ",".join(["0.1"] * len(ac.split(",")))),
                                            "GT"] + gts) + "\n")
                        n_rec += 1
                for start, side in sites:
                    width = 23 if side == "3" else 24
                    pos = start if strand == "+" else L - start - width
                    tab.write("%s\t%s\t%s\t%d\n" % (chrom, side, strand, pos))
    print("%d contigs, %d VCF records, %d four-allele positions" % (len(cases) * 2, n_rec, n_quad))


if __name__ == "__main__":
    main()
