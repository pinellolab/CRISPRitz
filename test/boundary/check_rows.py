#!/usr/bin/env python3
"""Every row a CRISPRitz search reports must be a stretch of its contig.

For each row of each targets file: the DNA column without its bulge gaps ('-'),
upper-cased, must be the contig's bases from Position on (reverse-complemented for
the - strand), and so lie inside the contig. A row made from bases the index does not
hold -- a position before the contig start, a base past its end or on an N run read as
the placeholder an index stores there -- fails. Independent of the oracle: it reads
only the FASTA and the targets file.

usage: check_rows.py <genome_dir> <targets.txt>...
Prints the number of rows checked; exits 1 naming the first bad rows.
"""
import os
import sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N", "Q": "Q"}


def main():
    genome_dir, files = sys.argv[1], sys.argv[2:]
    seqs = {}
    for name in os.listdir(genome_dir):
        if name.endswith(".fa"):
            with open(os.path.join(genome_dir, name)) as fh:
                head = fh.readline().strip()[1:].split()[0]
                seqs[head] = "".join(line.strip() for line in fh).upper()
    n = 0
    bad = []
    for path in files:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                col = line.rstrip("\n").split("\t")
                dna, chrom, pos, strand = col[2], col[3], int(col[4]), col[6]
                bases = dna.replace("-", "").upper()
                seq = seqs[chrom]
                n += 1
                if pos < 0 or pos + len(bases) > len(seq):
                    bad.append("%s: outside the contig (length %d): %s" % (path, len(seq), line.strip()))
                    continue
                ref = seq[pos:pos + len(bases)]
                if strand == "-":
                    ref = "".join(COMP[c] for c in reversed(ref))
                if ref != bases:
                    bad.append("%s: the contig holds %s there: %s" % (path, ref, line.strip()))
    if bad:
        sys.stdout.write("".join(b + "\n" for b in bad[:10]))
        sys.stdout.write("%d of %d rows are not their contig's bases\n" % (len(bad), n))
        sys.exit(1)
    print(n)


if __name__ == "__main__":
    main()
