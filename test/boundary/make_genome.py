#!/usr/bin/env python3
"""Synthetic genome for the CRISPRitz boundary test (test/run_boundary_test.sh).

One contig per case, so each planted site's distance to a boundary is exact. Every
site is the test guide (or a 1-4-mismatch or 1-bulge variant of it) next to a PAM
context that serves several PAMs at once. Written in target-strand orientation
(5' -> 3'), a 3' site is

    [distal side] guide [P[:t]] [PAM side]        P = TGGAGTCA

and a 5' site is

    [PAM side] [P5[-t:]] guide [distal side]      P5 = GCAGTTTA

so the first t context bases decide which PAMs it holds (3': NG from t = 2, NGG/NRG
from 3, NNGRRT from 6, an 8-wide PAM-less flank from 8; 5': TTTV from t = 4). Then:

  A  distal boundary: d = 0..4 random bases, then a contig end or an N run, beyond the
     guide's far end (the pad side; IDX-13), for perfect, 1- to 4-mismatch,
     DNA-bulge and RNA-bulge sites;
  B  PAM-side boundary: t = 0..8 context bases (t = 9: all 8 and one random base), then
     a contig end or an N run (IDX-15 at t = L, IDX-4 for t < K), perfect and
     3-mismatch sites;
  C  both sides tight: d = 0..2 and t in {2, 3, 4, 8}, contig ends on both sides;
  D  an N inside the guide region (never a site), and interior controls;
  E  for a search guide longer or shorter than the index's (docs/INDEX_FORMAT.md 5.1):
     the bases beyond the guide's far end are EXT (read outward), so the test's longer
     guides -- the guide plus EXT[:x] at its far end -- match there: d = 0..4 bases of
     EXT, then a contig end or an N run, perfect, 2-mismatch and 1-bulge sites, and an
     interior control; and the guide with its far x = 1, 2 bases cut off by a contig end
     or an N run, which a shorter guide matches but an index built for this guide's
     length does not hold.

Each case is written on the + strand as above and, reverse-complemented, on the -
strand, which also swaps contig start and end. A few IUPAC codes are sprinkled into
the random flanks so the variant-genome alphabet is exercised too.

usage: make_genome.py <out_dir> <3|5> [seed]
writes <out_dir>/genome/*.fa (one contig a file) and <out_dir>/cases.tsv
"""
import os
import random
import sys

GUIDE = "GACTTAGCCTAAGGTCTACG"
P3 = "TGGAGTCA"
P5 = "GCAGTTTA"
EXT = "TCAG"  # E cases: the bases beyond the guide's far end, read outward
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}
CHANGE = {"A": "C", "C": "A", "G": "T", "T": "G"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def main():
    out, side = sys.argv[1], sys.argv[2]
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 2026)

    def rnd(n, iupac=False):
        s = [rng.choice("ACGT") for _ in range(n)]
        if iupac and n > 6:
            for _ in range(n // 12):
                s[rng.randrange(n)] = rng.choice("RYSWKM")
        return "".join(s)

    # guide variants, with edits counted from the guide's far (distal) end: the 5' end
    # of the guide for a 3' PAM, its 3' end for a 5' PAM
    def at(k):
        return k if side == "3" else len(GUIDE) - 1 - k

    def variant(name):
        g = list(GUIDE)
        if name.startswith("m"):
            for k in (2, 7, 12, 17)[:int(name[1:])]:
                g[at(k)] = CHANGE[g[at(k)]]
        elif name == "dna1":  # one extra target base: a DNA bulge
            g.insert(at(3) + (0 if side == "3" else 1), "T")
        elif name == "rna1":  # one target base fewer: an RNA bulge
            del g[at(3)]
        return "".join(g)

    def ctx(t):  # the t PAM-side context bases next to the guide
        if side == "3":
            return P3[:t] if t <= 8 else P3 + rnd(t - 8)
        return P5[8 - t:] if t <= 8 else rnd(t - 8) + P5

    def site(name, d_side, t_side):
        """target-strand string: distal context + guide + PAM context (3'), or the
        mirror (5'); d_side / t_side are the strings beyond each end."""
        g = variant(name)
        if side == "3":
            return d_side + g + t_side
        return t_side + g + d_side

    cases = []

    def add(label, distal, guide_name, pamside):
        cases.append((label, site(guide_name, distal, pamside)))

    far = 25
    # A: distal boundary at d
    for btype in ("end", "N"):
        for d in range(0, 5):
            for name in ("m0", "m1", "m2", "m3", "m4", "dna1", "rna1"):
                inner = rnd(d)
                if side == "3":
                    distal = inner if btype == "end" else rnd(far, True) + "N" * 8 + inner
                    pamside = ctx(8) + rnd(far, True)
                else:
                    distal = inner if btype == "end" else inner + "N" * 8 + rnd(far, True)
                    pamside = rnd(far, True) + ctx(8)
                add(f"A_{btype}_d{d}_{name}", distal, name, pamside)
    # B: PAM-side boundary after t context bases
    for btype in ("end", "N"):
        for t in range(0, 10):
            for name in ("m0", "m3"):
                if side == "3":
                    pamside = ctx(t) + ("" if btype == "end" else "N" * 8 + rnd(far, True))
                    distal = rnd(far, True)
                else:
                    pamside = ("" if btype == "end" else rnd(far, True) + "N" * 8) + ctx(t)
                    distal = rnd(far, True)
                add(f"B_{btype}_t{t}_{name}", distal, name, pamside)
    # C: both ends tight
    for d in range(0, 3):
        for t in (2, 3, 4, 8):
            add(f"C_d{d}_t{t}", rnd(d), "m0", ctx(t))
    # D: an N in the guide region; interior controls
    g = list(GUIDE)
    g[10] = "N"
    cases.append(("D_guideN", rnd(far) + "".join(g) + ctx(8) + rnd(far) if side == "3"
                  else rnd(far) + ctx(8) + "".join(g) + rnd(far)))
    for name in ("m0", "m1", "m4", "dna1", "rna1"):
        add(f"D_interior_{name}", rnd(far, True), name, ctx(8) + rnd(far, True) if side == "3"
            else rnd(far, True) + ctx(8))

    # E: a known context beyond the far end (EXT), then a boundary; and cut guides
    def ext(d):  # the d bases beyond the far end, in genome (5' -> 3') order
        return EXT[:d][::-1] if side == "3" else EXT[:d]

    def pamside8():
        return ctx(8) + rnd(far, True) if side == "3" else rnd(far, True) + ctx(8)

    for btype in ("end", "N"):
        for d in range(0, 5):
            for name in ("m0", "m2", "dna1", "rna1"):
                if side == "3":
                    distal = ext(d) if btype == "end" else rnd(far, True) + "N" * 8 + ext(d)
                else:
                    distal = ext(d) if btype == "end" else ext(d) + "N" * 8 + rnd(far, True)
                add(f"E_{btype}_d{d}_{name}", distal, name, pamside8())
    add("E_interior_m0", rnd(far, True) + ext(4) if side == "3" else ext(4) + rnd(far, True),
        "m0", pamside8())
    for btype in ("end", "N"):
        for x in (1, 2):
            if side == "3":
                cut = GUIDE[x:]
                left = "" if btype == "end" else rnd(far, True) + "N" * 8
                cases.append((f"E_{btype}_cut{x}", left + cut + pamside8()))
            else:
                cut = GUIDE[:-x]
                right = "" if btype == "end" else "N" * 8 + rnd(far, True)
                cases.append((f"E_{btype}_cut{x}", pamside8() + cut + right))

    gdir = os.path.join(out, "genome")
    os.makedirs(gdir, exist_ok=True)
    with open(os.path.join(out, "cases.tsv"), "w") as tab:
        for k, (label, s) in enumerate(cases):
            for strand in "+-":
                name = f"c{k:03d}{'p' if strand == '+' else 'm'}_{label}"
                seq = s if strand == "+" else revcomp(s)
                with open(os.path.join(gdir, name + ".fa"), "w") as fh:
                    fh.write(f">{name}\n{seq}\n")
                tab.write(f"{name}\t{strand}\t{label}\t{len(seq)}\n")
    print(f"{len(cases) * 2} contigs written to {gdir}")


if __name__ == "__main__":
    main()
