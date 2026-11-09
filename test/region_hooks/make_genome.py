#!/usr/bin/env python3
"""Genomes, guides and BEDs for the region-hooks test (test/run_region_hooks_test.sh).

usage: make_genome.py <out_dir> [seed]

Writes, under <out_dir>:

  g3/*.fa       the 3'-PAM genome, one contig a file:
                  chrA      12 kb of random bases with IUPAC codes (about 4%) and Q, and a
                            site of the test guide (perfect, 2 mismatches, a DNA or an RNA
                            bulge; in a context that is NGG, NRG and NNGRRT) every 130
                            bases on alternating strands, some next to an N run on the pad
                            side (pad_valid < mb)
                  chrB_<j>  short contigs with a site at each end, on both strands: their
                            pads, or their 8-wide flanks, run off the contig (FLANK_OFF)
                  chrUn_KI270302v1_random  a contig whose name holds underscores
                  chrE      sites, and no interval of the BED
                  chrF      sites, and one interval over the whole contig
  g5/*.fa       the same for a 5' PAM (TTTV)
  fake25/fakechrA.fa, fake30/fakechrA.fa
                fake-indel genomes laid out as add-variants writes them: windows of LEFT
                reference bases, an ALT, RIGHT reference bases, one N between two windows;
                LEFT/RIGHT 25/27 (CRISPRitz up to 2.9.1) and 30/30 (2.10.0). Every third
                window (0, 3, 6, ...) is dense: in the BED, whole. Window 7 opens with the
                guide region of a site whose core, for a 22-nt guide and NNGRRT on a 20-nt
                index (28 bases), reaches back across the separator N into window 6:
                a 2-RNA-bulge alignment of it lies in window 7 and, with LEFT = 25 (not
                30), covers its ALT (docs/INDEX_FORMAT.md section 13.5)
  <g>.bed       each genome's region BED: one base next to or on each planted site's core,
                at offsets that cycle through its pad, the edges of its guide region and
                its PAM, and the distal flank of an 8-wide window; random single bases;
                three long intervals; chrF whole; nothing on chrE. The fake genomes:
                their dense windows, whole
  <g>.windows.tsv  the fake genomes' windows: contig, start, end, left, right, dense
  guides_*.txt  the search guides, one file a search PAM file of the runner
"""
import os
import random
import sys

G3 = "GACTTAGCCTAAGGTCTACG"  # the 3' test guide (test/boundary's)
P3 = "TGGAGTCA"  # 3' context: NGG, NRG and NNGRRT, and an 8-wide flank
G5 = "CTGAGCTTACGGATCAGCTA"  # the 5' test guide
P5 = "GCAGTTTA"  # 5' context: TTTA (TTTV) next to the guide
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N", "Q": "Q"}
CHANGE = {"A": "C", "C": "A", "G": "T", "T": "G"}
CODES = "RYSWKMBDHV"
NAMES = ("m0", "m2", "dna1", "rna1", "m0")
# The BED base planted at each site, cycling: (anchor, delta) on the site's own strand
# and orientation, read outward from the PAM: anchor 'g' is the guide region's far
# (distal) base, 'p' the PAM's first (proximal) base; delta > 0 moves away from the
# guide on the pad side ('g') or the PAM side ('p'), delta < 0 into the guide region.
# They land on the pad (g+1, g+2), beyond it (g+3, g+5, g+9), the guide region's edges
# (g0, g-1, p-1), the PAM's bases (p0..p+2), a longer PAM's (p+3, p+5) and the distal
# flank of an 8-wide window (p+6, p+7), and beyond (p+8, p+10).
REL = [("g", 1), ("g", 2), ("p", 7), ("g", 0), ("p", 0), ("p", 3), ("g", 3), ("p", -1),
       ("p", 2), ("g", -1), ("p", 6), ("p", 5), ("p", 1), ("g", 5), ("p", 8), ("g", 9),
       ("p", 10)]


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def main():
    out = sys.argv[1]
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 20261109)

    def rnd(n, codes=0.0, q=0.0):
        s = []
        for _ in range(n):
            x = rng.random()
            s.append(rng.choice(CODES) if x < codes else "Q" if x < codes + q else rng.choice("ACGT"))
        return "".join(s)

    def genome(side):
        guide, ctx = (G3, P3) if side == "3" else (G5, P5)

        def variant(name):
            """the guide with its edits counted from its far (distal) end"""
            at = (lambda k: k) if side == "3" else (lambda k: len(guide) - 1 - k)
            g = list(guide)
            if name == "m2":
                for k in (3, 11):
                    g[at(k)] = CHANGE[g[at(k)]]
            elif name == "dna1":
                g.insert(at(5) + (0 if side == "3" else 1), "T")
            elif name == "rna1":
                del g[at(5)]
            return "".join(g)

        def site(name, distal, t=8):
            """(target-strand string, offset of its guide region, guide length): the
            guide with t bases of the PAM-side context and the given distal bases"""
            g = variant(name)
            if side == "3":
                return distal + g + ctx[:t], len(distal), len(g)
            return ctx[len(ctx) - t:] + g + distal, t, len(g)

        counter = [0]

        def place(chars, pos, s, goff, glen, strand, bed, contig):
            """write site s at + strand position pos, and one BED base by REL"""
            seq = s if strand == "+" else revcomp(s)
            chars[pos:pos + len(seq)] = list(seq)
            g0 = pos + goff if strand == "+" else pos + len(s) - goff - glen  # + strand start
            g1 = g0 + glen
            pam_right = (strand == "+") == (side == "3")
            # far (distal) base of the guide region and the PAM's first base, + strand
            far, pam0, out_pad, out_pam = (g0, g1, -1, 1) if pam_right else (g1 - 1, g0 - 1, 1, -1)
            anchor, delta = REL[counter[0] % len(REL)]
            counter[0] += 1
            b = far + out_pad * delta if anchor == "g" else pam0 + out_pam * delta
            if 0 <= b < len(chars):
                bed.append((contig, b, b + 1))

        contigs, bed = [], []
        # chrA
        L = 12000
        chars = list(rnd(L, 0.04, 0.005))
        n = 0
        for pos in range(60, L - 80, 130):
            name = NAMES[n % len(NAMES)]
            strand = "+-"[n % 2]
            distal = rnd(4)
            if n % 9 == 4:  # an N run 0, 1 or 2 bases beyond the guide's far end
                k = n % 3
                distal = ("NNNN" + rnd(k)) if side == "3" else (rnd(k) + "NNNN")
            s, goff, glen = site(name, distal)
            place(chars, pos, s, goff, glen, strand, bed, "chrA")
            n += 1
        for p in range(L):
            if rng.random() < 1.0 / 70:
                bed.append(("chrA", p, p + 1))
        for a in (1500, 5200, 9100):
            bed.append(("chrA", a, a + rng.randrange(150, 450)))
        contigs.append(("chrA", "".join(chars)))

        # chrB_<j>: a site at each end, one per strand; even j: the PAM-side context cut
        # to t < 8 bases at the contig end (an 8-wide flank off the contig), odd j: the
        # far end of the guide 0 or 1 base from the contig end (the pad off it)
        for j in range(10):
            name = NAMES[j % len(NAMES)]
            mid = rnd(5 + j, 0.05)
            if j % 2 == 0:
                s, goff, glen = site(name, rnd(4), 3 + j % 5)
                pam_end_first = side == "5"  # the site whose PAM side faces the contig start
            else:
                s, goff, glen = site(name, rnd(j % 4 // 2))
                pam_end_first = side == "3"
            # first site: its PAM side (even j) or far side (odd j) at the contig start
            first_strand = "+" if pam_end_first else "-"
            second_strand = "-" if first_strand == "+" else "+"
            first = s if first_strand == "+" else revcomp(s)
            second = s if second_strand == "+" else revcomp(s)
            chars = list(first + mid + second)
            name_c = "chrB_%d" % j
            place(chars, 0, s, goff, glen, first_strand, bed, name_c)
            place(chars, len(first) + len(mid), s, goff, glen, second_strand, bed, name_c)
            contigs.append((name_c, "".join(chars)))

        # chrUn_KI270302v1_random, chrE, chrF
        for name_c, L, step in (("chrUn_KI270302v1_random", 1800, 200), ("chrE", 1500, 170),
                                ("chrF", 1500, 170)):
            chars = list(rnd(L, 0.03))
            dropped = []
            for i, pos in enumerate(range(40, L - 80, step)):
                s, goff, glen = site(NAMES[i % len(NAMES)], rnd(4))
                place(chars, pos, s, goff, glen, "+-"[i % 2], dropped if name_c != "chrUn_KI270302v1_random" else bed, name_c)
            contigs.append((name_c, "".join(chars)))
        bed.append(("chrF", 0, len(contigs[-1][1])))
        return contigs, bed

    def write(dirname, contigs):
        d = os.path.join(out, dirname)
        os.makedirs(d, exist_ok=True)
        for name, seq in contigs:
            with open(os.path.join(d, name + ".fa"), "w") as fh:
                fh.write(">%s\n" % name)
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")

    def write_bed(name, bed):
        with open(os.path.join(out, name), "w") as fh:
            fh.write("# region-hooks test BED\n")
            for c, s, e in sorted(bed):
                fh.write("%s\t%d\t%d\n" % (c, s, e))

    for side in "35":
        contigs, bed = genome(side)
        write("g" + side, contigs)
        write_bed("g%s.bed" % side, bed)

    # fake-indel genomes: the same windows cut at 25/27 and at 30/30
    n_win = 60
    wins = [[rnd(30, 0.05), rng.choice("ACGT") + rnd(rng.randrange(0, 6)), rnd(30, 0.05)]
            for _ in range(n_win)]
    ngg_guides = []
    for i in (2, 9, 11, 18, 19, 27, 35, 42, 44, 52):  # NGG sites; 9, 18, 27, 42 in dense windows
        g = rnd(20)
        wins[i][0] = wins[i][0][:5] + g + "TGG" + wins[i][0][28:]
        ngg_guides.append(g)
    gr = "ACGTTGCAGTCCATGACTGA"  # window 7's guide region
    for LEFT, RIGHT, tag in ((25, 27, "fake25"), (30, 30, "fake30")):
        seq, rows = [], []
        pos = 0
        for i, (left, alt, right) in enumerate(wins):
            lft = left[30 - LEFT:]
            if i == 7:  # NNGRRT (AAGAAT) at bases 20-25: its last base the ALT's when LEFT = 25
                if LEFT == 25:
                    lft, alt = gr + "AAGAA", "T" + alt[1:]
                else:
                    lft = gr + "AAGAAT" + lft[26:]
            w = lft + alt + right[:RIGHT]
            if i:
                seq.append("N")
                pos += 1
            rows.append(("fakechrA", pos, pos + len(w), LEFT, RIGHT, 1 if i % 3 == 0 else 0))
            seq.append(w)
            pos += len(w)
        write(tag, [("fakechrA", "".join(seq))])
        with open(os.path.join(out, tag + ".windows.tsv"), "w") as fh:
            for r in rows:
                fh.write("%s\t%d\t%d\t%d\t%d\t%d\n" % r)
        write_bed(tag + ".bed", [(c, s, e) for c, s, e, _, _, d in rows if d])

    def gfile(name, guides, pamlen, side):
        with open(os.path.join(out, name), "w") as fh:
            for g in guides:
                fh.write((g + "N" * pamlen if side == "3" else "N" * pamlen + g) + "\n")

    r3, r5 = [rnd(20), rnd(20)], [rnd(20), rnd(20)]
    gfile("guides_ngg20.txt", [G3] + r3, 3, "3")
    gfile("guides_ngg21.txt", ["T" + G3, "A" + r3[0]], 3, "3")
    gfile("guides_ngg19.txt", [G3[1:], r3[1][1:]], 3, "3")
    gfile("guides_nngrrt21.txt", ["T" + G3, "C" + r3[1]], 6, "3")
    gfile("guides_n8.txt", [G3], 8, "3")
    gfile("guides_tttv20.txt", [G5] + r5, 4, "5")
    gfile("guides_tttv23.txt", [G5 + "GAC", r5[0] + "TTA"], 4, "5")
    gfile("guides_fake_ngg20.txt", ngg_guides, 3, "3")
    gfile("guides_fake_nngrrt21.txt", ["C" + gr], 6, "3")  # n = 27
    gfile("guides_fake_nngrrt22.txt", ["TC" + gr], 6, "3")  # n = 28
    gfile("guides_fake_nngrrt23.txt", ["GTC" + gr], 6, "3")  # n = 29
    print("genomes written to %s" % out)


if __name__ == "__main__":
    main()
