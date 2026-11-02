#!/usr/bin/env python3
"""Independent brute-force oracle for the CRISPRitz boundary test.

Lists every row a CRISPRitz search must report for a genome, straight from the FASTA
files -- no index, no CRISPRitz code -- in searchTST's 10-column format:

    bulge type, crRNA, DNA, chromosome, position, cluster position, direction,
    mismatches, bulge size, total

What counts as a target follows docs/INDEX_FORMAT.md section 3 (index format 2):

* a site is a PAM occurrence on either strand whose guide region lies inside the
  contig and holds no N;
* its PAM: the search PAM's L bases must lie inside the contig and match (a genomic N
  matches no PAM letter) -- for a dedicated index (index PAM = search PAM) and for
  --prune-pam off a PAM-less index alike; without --prune-pam on a PAM-less (all-N)
  index of width K the K flank bases must lie inside the contig, whatever they hold;
* the pad, the mb bases beyond the guide's far end, may only be consumed (by DNA
  bulges) up to the first base that is off the contig or N.

How a guide aligns to a site follows searchTST's search (nearsearch): reading guide and
target from the PAM outward, each guide base is a match (their IUPAC bits meet; no
choice then), a mismatch, or an RNA bulge (guide base unpaired); a target base may be a
DNA bulge (target base unpaired); the alignment ends as soon as the last guide base is
placed. On the stored string's last base (the end of guide + pad) only a match or a
mismatch of the last guide base is possible, as there the tree node is a leaf.

usage: oracle.py --genome DIR --guides FILE --pam FILE --index-pam FILE --mb N
                 --mm N --bdna N --brna N [--prune]
Prints the rows, sorted, without a header.
"""
import argparse
import os
import sys

BITS = {"A": 1, "C": 2, "G": 4, "T": 8, "R": 5, "Y": 10, "S": 6, "W": 9, "K": 12,
        "M": 3, "B": 14, "D": 13, "H": 11, "V": 7, "N": 15}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}
OFF = None  # a position off the contig


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def read_pam_file(path):
    """(pattern, side, width): pattern = the PAM letters, side '3' or '5'."""
    with open(path) as fh:
        motif, limit = fh.readline().split()
    motif, limit = motif.upper(), int(limit)
    if limit < 0:
        return motif[:-limit], "5", -limit, len(motif)
    return motif[len(motif) - limit:], "3", limit, len(motif)


def read_contigs(genome_dir):
    for name in sorted(os.listdir(genome_dir)):
        if not name.endswith(".fa"):
            continue
        with open(os.path.join(genome_dir, name)) as fh:
            head = fh.readline().strip()[1:].split()[0]
            seq = "".join(line.strip() for line in fh).upper()
        yield head, seq


def pam_ok(bases, pattern):
    """bases (None = off the contig) match an IUPAC pattern; a genomic N never does."""
    for b, p in zip(bases, pattern):
        if b is OFF or b == "N" or not BITS[b] & BITS[p]:
            return False
    return True


def alignments(T, g, G, mb, pad_valid, mm, bdna, brna):
    """Every alignment searchTST reports between guide g and stored string T (both
    read from the PAM outward; T is guide region + pad, length G + mb). Yields
    (guide columns, target columns, mismatches, dna bulges, rna bulges)."""
    L = G + mb
    limit = G + pad_valid  # target bases an alignment may consume
    out = []

    def rec(i, j, d, bD, bR, gc, tc):
        if i == G:  # the last guide base is placed: report
            out.append((gc, tc, mm - d, bdna - bD, brna - bR))
            return
        if j >= limit:  # the next target base is not there
            # only an RNA bulge (which consumes no target base) can go on
            if j < L - 1 and bR > 0:
                rec(i + 1, j, d, bD, bR - 1, gc + g[i], tc + "-")
            return
        t = T[j]
        match = BITS[g[i]] & BITS[t]
        if j < L - 1:  # an inner tree node
            if match:
                rec(i + 1, j + 1, d, bD, bR, gc + g[i], tc + t)
            elif d > 0:
                rec(i + 1, j + 1, d - 1, bD, bR, gc + g[i], tc + t.lower())
            if bR > 0:
                rec(i + 1, j, d, bD, bR - 1, gc + g[i], tc + "-")
            if bD > 0:
                rec(i, j + 1, d, bD - 1, bR, gc + "-", tc + t)
        elif i == G - 1:  # the leaf: the last guide base on the last stored base
            if match:
                out.append((gc + g[i], tc + t, mm - d, bdna - bD, brna - bR))
            elif d > 0:
                out.append((gc + g[i], tc + t.lower(), mm - d + 1, bdna - bD, brna - bR))

    rec(0, 0, mm, bdna, brna, "", "")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--guides", required=True)
    ap.add_argument("--pam", required=True, help="the search PAM file")
    ap.add_argument("--index-pam", required=True, help="the PAM file the index was built with")
    ap.add_argument("--mb", type=int, required=True, help="the index's max bulges")
    ap.add_argument("--mm", type=int, required=True)
    ap.add_argument("--bdna", type=int, required=True)
    ap.add_argument("--brna", type=int, required=True)
    ap.add_argument("--prune", action="store_true")
    a = ap.parse_args()

    spat, side, L, slen = read_pam_file(a.pam)
    ipat, iside, K, _ = read_pam_file(a.index_pam)
    assert side == iside, "PAM sides differ"
    pamless = set(ipat) == {"N"}
    G = slen - L
    mb = a.mb
    if not a.prune:
        assert K == L, "without --prune-pam the search PAM must have the index's width"
    with open(a.guides) as fh:
        guides = [line.strip().upper() for line in fh if line.strip()]

    # the PAM bases a site needs, and those its row reports: under --prune-pam, and on
    # a dedicated index, the search PAM's L bases (matching); without --prune-pam on a
    # PAM-less index, all K bases of the flank (on the contig, any content)
    flank_only = pamless and not a.prune
    width = K if flank_only else L

    rows = set()
    for chrom, seq in read_contigs(a.genome):
        n = len(seq)
        for strand in "+-":
            Q = seq if strand == "+" else revcomp(seq)  # the target strand, 5' to 3'

            def at(x):
                return Q[x] if 0 <= x < n else OFF

            for g0 in range(0, n - G + 1):  # guide region Q[g0 : g0 + G]
                guide_region = Q[g0:g0 + G]
                if "N" in guide_region:
                    continue
                if side == "3":
                    pam = [at(x) for x in range(g0 + G, g0 + G + width)]
                    pad = [at(x) for x in range(g0 - 1, g0 - 1 - mb, -1)]
                    T = guide_region[::-1]
                else:
                    pam = [at(x) for x in range(g0 - width, g0)]
                    pad = [at(x) for x in range(g0 + G, g0 + G + mb)]
                    T = guide_region
                if flank_only:
                    if any(b is OFF for b in pam):
                        continue
                elif not pam_ok(pam, spat):
                    continue
                pad_valid = 0
                while pad_valid < mb and pad[pad_valid] not in (OFF, "N"):
                    pad_valid += 1
                T = T + "".join(pad[:pad_valid]) + "A" * (mb - pad_valid)
                pam_str = "".join(pam)
                for line in guides:
                    guide = line[:G] if side == "3" else line[L:L + G]
                    g = guide[::-1] if side == "3" else guide
                    for gc, tc, m, dna, rna in alignments(T, g, G, mb, pad_valid,
                                                         a.mm, a.bdna, a.brna):
                        consumed = G + dna - rna  # target bases of the alignment
                        if side == "3":
                            crrna = gc[::-1] + "N" * width
                            dna_s = tc[::-1] + pam_str
                            if strand == "+":
                                pos = g0 + G - consumed
                                cluster = g0
                            else:
                                pos = n - g0 - G - width
                                cluster = pos
                        else:
                            crrna = "N" * width + gc
                            dna_s = pam_str + tc
                            if strand == "+":
                                pos = g0 - width
                                cluster = pos
                            else:
                                pos = n - g0 - consumed
                                cluster = n - g0 - G
                        if dna and rna:
                            btype = "RNA,DNA"
                        elif rna:
                            btype = "RNA"
                        elif dna:
                            btype = "DNA"
                        else:
                            btype = "X"
                        rows.add("\t".join(map(str, (btype, crrna, dna_s, chrom, pos, cluster,
                                                     strand, m, dna + rna, m + dna + rna))))
    sys.stdout.write("".join(r + "\n" for r in sorted(rows)))


if __name__ == "__main__":
    main()
