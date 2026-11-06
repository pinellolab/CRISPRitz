#!/usr/bin/env python3
"""Multi-allelic enrichment: unit tests run against BOTH enricher ports.

A normalised VCF (`bcftools norm -m -any`, how every CRISPRme panel is built) splits
a multi-allelic site into one record per ALT. Up to CRISPRitz 2.9.1 each record
rebuilt the position's IUPAC code from the reference base and overwrote its SNP dict
entry, so only the LAST record's ALT reached the enriched genome and the dict. These
tests pin the merged behaviour on small synthetic inputs, with exact expected outputs:

  * every ALT of a split site is coded, in any record order, and the code equals the
    one the same site gets written as a single record listing every ALT;
  * a position whose reference and ALTs cover all four bases keeps the reference and
    its most frequent ALTs (N is a gap to the search, not "any base");
  * the dict lists every ALT, each entry with its own record's carriers, genotypes,
    rsID and AF, joined by '$' like a multi-ALT record's value;
  * non-PASS records, duplicate ALTs, '*', haploid / unphased / missing genotypes,
    non-adjacent records, an indel at the same position, a reference N or IUPAC base;
  * the indel log and fake-indel contigs do not move;
  * a position with one record behaves exactly as before.

Each test runs the real program as a subprocess: `enricher.py` under this interpreter
(needs pandas; ENRICHER_PY overrides the path) and the compiled `enricher`
(ENRICHER_BIN, else compiled here from enricher.cpp with $CXX / c++ / g++ and -lz). A run must write nothing to stderr: a
CRISPRme stage fails on any.

Usage: python3 test/test_enricher_multiallelic.py [-v]
"""

import gzip
import itertools
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENRICH_DIR = os.path.join(REPO, "sourceCode", "Python_Scripts", "Enrichment")
ENRICHER_PY = os.environ.get("ENRICHER_PY") or os.path.join(ENRICH_DIR, "enricher.py")
ENRICHER_CPP = os.path.join(ENRICH_DIR, "enricher.cpp")

CHROM = "chrT"
SAMPLES = ["S1", "S2", "S3", "S4"]

BITS = {"A": 1, "C": 2, "G": 4, "T": 8}
LETTER = {3: "M", 5: "R", 9: "W", 6: "S", 10: "Y", 12: "K",
          14: "B", 13: "D", 11: "H", 7: "V", 15: "N"}


def code_of(bases):
    """IUPAC letter of a set of A/C/G/T bases (one base -> that base)."""
    mask = 0
    for b in bases:
        mask |= BITS.get(b, 0)
    if mask in (1, 2, 4, 8):
        return [k for k, v in BITS.items() if v == mask][0]
    return LETTER[mask]


def make_reference(length=400, seed=11):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(length))


REF = make_reference()


def ref1(pos1):
    return REF[pos1 - 1]


def others(pos1):
    """The three bases that differ from the reference base at pos1, in A/C/G/T order."""
    return [b for b in "ACGT" if b != ref1(pos1)]


def rec(pos, rsid, ref, alt, info, gts, filt="PASS", fmt="GT"):
    return "\t".join([CHROM, str(pos), rsid, ref, alt, ".", filt, info, fmt] + list(gts))


def _build_cpp(workdir):
    cxx = os.environ.get("CXX") or shutil.which("c++") or shutil.which("g++")
    if not cxx:
        return None
    out = os.path.join(workdir, "enricher")
    r = subprocess.run([cxx, "-std=c++11", "-O2", ENRICHER_CPP, "-o", out, "-lz"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return out if r.returncode == 0 else None


_BUILD_DIR = tempfile.mkdtemp(prefix="enricher_build_")
CPP_BIN = os.environ.get("ENRICHER_BIN") or _build_cpp(_BUILD_DIR)


class Run(object):
    """One enricher run over ``records`` (a list of VCF data lines) in its own dir."""

    def __init__(self, cmd, records, mode="yes", fasta=REF, header=None):
        self.dir = tempfile.mkdtemp(prefix="enr_")
        os.mkdir(os.path.join(self.dir, "SNPs_genome"))
        fa = os.path.join(self.dir, CHROM + ".fa")
        with open(fa, "w") as f:
            f.write(">" + CHROM + "\n")
            for i in range(0, len(fasta), 60):
                f.write(fasta[i:i + 60] + "\n")
        lines = header or [
            "##fileformat=VCFv4.2",
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES),
        ]
        plain = os.path.join(self.dir, "variants.vcf")
        with open(plain, "w") as f:
            f.write("\n".join(lines + list(records)) + "\n")
        vcf = plain + ".gz"
        with open(plain, "rb") as src, gzip.open(vcf, "wb") as dst:
            dst.write(src.read())
        p = subprocess.run(cmd + [vcf, fa, "res", mode, plain], cwd=self.dir,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           universal_newlines=True)
        self.returncode, self.stdout, self.stderr = p.returncode, p.stdout, p.stderr
        self.mode = mode

    def path(self, *parts):
        return os.path.join(self.dir, *parts)

    @property
    def genome(self):
        with open(self.path("SNPs_genome", "res_enriched", CHROM + ".enriched.fa")) as f:
            f.readline()
            return f.read().rstrip("\n")

    def base(self, pos1):
        return self.genome[pos1 - 1]

    @property
    def dict(self):
        with open(self.path("SNPs_genome", "my_dict_" + CHROM + ".json")) as f:
            return json.load(f)

    def entries(self, pos1):
        return self.dict[CHROM + "," + str(pos1)].split("$")

    def raw(self, rel):
        with open(self.path(*rel.split("/")), "rb") as f:
            return f.read()

    def cleanup(self):
        shutil.rmtree(self.dir, ignore_errors=True)


class EnricherMultiallelicTests(object):
    """Mixed into one TestCase per port; ``cmd`` is the port's command prefix."""

    cmd = None

    def run_enricher(self, records, **kw):
        r = Run(self.cmd, records, **kw)
        self.addCleanup(r.cleanup)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(r.stderr, "", "a run must write nothing to stderr")
        return r

    # ---- the enriched genome -------------------------------------------------

    def test_split_two_alts_every_order(self):
        p = 50
        a1, a2 = others(p)[:2]
        r1 = rec(p, "rs1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"])
        r2 = rec(p, "rs2", ref1(p), a2, "AC=1;AF=0.2", ["0|0", "1|0", "0|0", "0|0"])
        want = code_of({ref1(p), a1, a2})
        for order in ([r1, r2], [r2, r1]):
            self.assertEqual(self.run_enricher(order).base(p), want)

    def test_split_code_equals_joined_code_every_order(self):
        """2 and 3 ALTs: every permutation of split records == the joined record."""
        for p, n in ((60, 2), (70, 3)):
            alts = others(p)[:n]
            afs = ["0.3", "0.2", "0.1"][:n]
            split = [rec(p, "rs%d" % k, ref1(p), a, "AC=1;AF=" + af,
                         ["0|1", "0|0", "0|0", "0|0"]) for k, (a, af) in enumerate(zip(alts, afs))]
            joined = self.run_enricher([rec(p, "rsJ", ref1(p), ",".join(alts),
                                            "AC=" + ",".join("1" * n) + ";AF=" + ",".join(afs),
                                            ["0|1", "0|2", "0|0", "0|0"])])
            for order in itertools.permutations(split):
                got = self.run_enricher(list(order))
                self.assertEqual(got.genome, joined.genome, "order %s" % [o.split("\t")[4] for o in order])

    def test_four_bases_keep_reference_and_most_frequent_alts(self):
        p = 80
        a, b, c = others(p)
        afs = {a: "0.05", b: "0.0001", c: "0.3"}  # b is the rarest -> left out
        split = [rec(p, "rs" + x, ref1(p), x, "AC=1;AF=" + afs[x], ["0|1", "0|0", "0|0", "0|0"])
                 for x in (a, b, c)]
        want = code_of({ref1(p), a, c})
        for order in itertools.permutations(split):
            r = self.run_enricher(list(order))
            self.assertEqual(r.base(p), want)
            self.assertNotEqual(r.base(p), "N")
            self.assertIn("NOTE: 1 position(s) of %s carry all four bases" % CHROM, r.stdout)
            # the dict still lists all three, in file order
            self.assertEqual([e.split(";")[1] for e in r.entries(p)],
                             [ref1(p) + "," + o.split("\t")[4] for o in order])

    def test_four_bases_tie_leaves_out_later_base(self):
        p = 90
        a, b, c = others(p)  # A/C/G/T order; all AF equal -> drop c (the latest)
        split = [rec(p, "rs" + x, ref1(p), x, "AC=1;AF=0.01", ["0|1", "0|0", "0|0", "0|0"])
                 for x in (c, a, b)]
        self.assertEqual(self.run_enricher(split).base(p), code_of({ref1(p), a, b}))
        # and with AF absent from the INFO field: every AF reads 0 -> same tie rule
        split = [rec(p, "rs" + x, ref1(p), x, "AC=1", ["0|1", "0|0", "0|0", "0|0"])
                 for x in (b, c, a)]
        self.assertEqual(self.run_enricher(split).base(p), code_of({ref1(p), a, b}))

    def test_four_bases_af_spellings(self):
        """Exponent, '.', missing AF: '.' and a missing value read as 0."""
        p = 100
        a, b, c = others(p)
        recs = [rec(p, "rsa", ref1(p), a, "AC=1;AF=2e-1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "rsb", ref1(p), b, "AC=1;AF=.", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "rsc", ref1(p), c, "AC=1;AF=1.5E-3", ["0|1", "0|0", "0|0", "0|0"])]
        self.assertEqual(self.run_enricher(recs).base(p), code_of({ref1(p), a, c}))

    def test_four_bases_ranked_by_af_fields(self):
        """No INFO/AF (a sites-only merge such as CRISPRme's mega: AF_<dataset>, AF_max):
        the ALTs rank by their highest AF_* value, not all 0 (which would leave the choice
        to base order); with AF and AF_* both present the higher counts."""
        p = 230
        a, b, c = others(p)
        # base order alone would leave out c; AF_max says b is the rarest
        recs = [rec(p, "m" + x, ref1(p), x, info, ["0|0", "0|0", "0|0", "0|0"])
                for x, info in ((a, "AF_one=0.01;AF_max=0.01"), (b, "AF_two=0.0002;AF_max=0.0002"),
                                (c, "AF_one=0.2;AF_two=0.3;AF_max=0.3"))]
        r = self.run_enricher(recs, header=None)
        self.assertEqual(r.base(p), code_of({ref1(p), a, c}))
        # AF 0.001 but AF_pop 0.4 for a: a is common somewhere, so c (0.002) is left out
        q = 235
        a, b, c = others(q)
        recs = [rec(q, "n" + x, ref1(q), x, info, ["0|1", "0|0", "0|0", "0|0"])
                for x, info in ((a, "AC=1;AF=0.001;AF_pop=0.4"), (b, "AC=1;AF=0.01"), (c, "AC=1;AF=0.002"))]
        self.assertEqual(self.run_enricher(recs).base(q), code_of({ref1(q), a, b}))

    def test_four_bases_joined_record_no_longer_writes_N(self):
        """A record listing three SNP ALTs used to write N, a gap to the search."""
        p = 110
        a, b, c = others(p)
        r = self.run_enricher([rec(p, "rsJ", ref1(p), ",".join([a, b, c]),
                                   "AC=1,1,1;AF=0.2,0.3,0.001", ["0|1", "0|2", "3|0", "0|0"])])
        self.assertEqual(r.base(p), code_of({ref1(p), a, b}))

    def test_four_bases_alt_on_two_records_counts_highest_af(self):
        p = 115
        a, b, c = others(p)
        recs = [rec(p, "r1", ref1(p), a, "AC=1;AF=0.001", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "r2", ref1(p), b, "AC=1;AF=0.01", ["0|0", "0|1", "0|0", "0|0"]),
                rec(p, "r3", ref1(p), c, "AC=1;AF=0.1", ["0|0", "0|0", "0|1", "0|0"]),
                rec(p, "r4", ref1(p), a, "AC=1;AF=0.5", ["0|0", "0|0", "0|0", "0|1"])]
        # a counts with 0.5, so b (0.01) is the rarest
        self.assertEqual(self.run_enricher(recs).base(p), code_of({ref1(p), a, c}))

    def test_reference_N_and_reference_iupac(self):
        fasta = list(REF)
        fasta[119] = "N"   # pos 120
        fasta[129] = "R"   # pos 130 (A/G in the FASTA itself)
        fasta = "".join(fasta)
        recs = [rec(120, "n1", "N", "A", "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(120, "n2", "N", "C", "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0|0"]),
                rec(130, "y1", "R", "C", "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(130, "y2", "R", "T", "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0|0"])]
        r = self.run_enricher(recs, fasta=fasta)
        # a reference N or IUPAC base adds no base of its own (the rule of 2.9.1)
        self.assertEqual(r.base(120), "M")
        self.assertEqual(r.base(130), "Y")
        # one record alone leaves a single base -> the FASTA letter stays, as before
        r1 = self.run_enricher(recs[:1] + recs[2:3], fasta=fasta)
        self.assertEqual(r1.base(120), "N")
        self.assertEqual(r1.base(130), "R")
        # four ALTs on a reference N: the rarest of the four is left out
        recs4 = [rec(120, "n" + x, "N", x, "AC=1;AF=" + af, ["0|1", "0|0", "0|0", "0|0"])
                 for x, af in (("A", "0.4"), ("C", "0.3"), ("G", "0.01"), ("T", "0.2"))]
        self.assertEqual(self.run_enricher(recs4, fasta=fasta).base(120), "H")

    def test_genome_is_the_same_without_the_dict(self):
        """Mode 'no' (genome only) codes every position exactly as mode 'yes'."""
        p, q = 140, 145
        a, b, c = others(p)
        recs = [rec(p, "r" + x, ref1(p), x, "AC=1;AF=0.0" + str(i + 1), ["0|1", "0|0", "0|0", "0|0"])
                for i, x in enumerate((c, b, a))]
        recs += [rec(q, "q1", ref1(q), others(q)[0], "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                 rec(q, "q2", ref1(q), others(q)[2], "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"])]
        yes = self.run_enricher(recs, mode="yes")
        no = self.run_enricher(recs, mode="no")
        self.assertEqual(yes.genome, no.genome)
        self.assertFalse(os.path.exists(no.path("SNPs_genome", "my_dict_" + CHROM + ".json")))

    # ---- the SNP dict ----------------------------------------------------------

    def test_dict_entries_keep_each_records_carriers_rsid_af(self):
        p = 150
        a1, a2 = others(p)[:2]
        recs = [rec(p, "rsA", ref1(p), a1, "AC=2;AF=0.25", ["0|1", "1|0", "0|0", "0|0"]),
                rec(p, "rsB;rsC", ref1(p), a2, "AC=1;AF=0.125", ["1|0", "0|0", "0|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.dict[CHROM + "," + str(p)],
                         "S1:0|1,S2:1|0;%s,%s;rsA;0.25$S1:1|0;%s,%s;rsB,rsC;0.125"
                         % (ref1(p), a1, ref1(p), a2))

    def test_dict_genotype_shapes(self):
        """Haploid, unphased, missing and GT:DP genotypes on split records."""
        p = 160
        a1, a2, a3 = others(p)
        recs = [rec(p, "h1", ref1(p), a1, "AC=1;AF=0.1", ["1", "0", ".", "0"]),
                rec(p, "h2", ref1(p), a2, "AC=1;AF=0.2", ["0/1", "./.", "1/.", "0/0"]),
                rec(p, "h3", ref1(p), a3, "AC=1;AF=0.3", [".|1:7", "0|0:3", "0|0:1", "1|1:9"],
                    fmt="GT:DP")]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), [
            "S1:1;%s,%s;h1;0.1" % (ref1(p), a1),
            "S1:0/1,S3:1/.;%s,%s;h2;0.2" % (ref1(p), a2),
            "S1:.|1,S4:1|1;%s,%s;h3;0.3" % (ref1(p), a3),
        ])

    def test_dict_passes_conflicting_phase_through(self):
        """Split records phased independently (as in HGDP + 1000 Genomes) can put two ALTs
        of one position on ONE haplotype slot: 0|1 under one ALT and 0|1 under the other
        (chr22:10661629, NA18869), 1|0 and 1|0, or more ALT copies than the slots hold
        (1|1 and 0|1), up to three ALTs on one slot. The enricher copies every record's
        genotypes into its entry verbatim -- it neither resolves the conflict nor drops a
        carrier -- and CRISPRme has to (docs/MULTIALLELIC_ENRICHMENT.md, section 6.3). S4
        is the consistent case, one ALT per slot."""
        p = 240
        a1, a2 = others(p)[:2]
        recs = [rec(p, "c1", ref1(p), a1, "AC=5;AF=0.625", ["0|1", "1|0", "1|1", "1|0"]),
                rec(p, "c2", ref1(p), a2, "AC=4;AF=0.5", ["0|1", "1|0", "0|1", "0|1"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), [
            "S1:0|1,S2:1|0,S3:1|1,S4:1|0;%s,%s;c1;0.625" % (ref1(p), a1),
            "S1:0|1,S2:1|0,S3:0|1,S4:0|1;%s,%s;c2;0.5" % (ref1(p), a2)])
        self.assertEqual(r.base(p), code_of({ref1(p), a1, a2}))
        q = 245
        b1, b2, b3 = others(q)
        recs = [rec(q, "t" + x, ref1(q), x, "AC=1;AF=" + af, ["0|1", "0|0", "0|0", "0|0"])
                for x, af in ((b1, "0.3"), (b2, "0.001"), (b3, "0.2"))]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(q), ["S1:0|1;%s,%s;t%s;%s" % (ref1(q), x, x, af)
                                        for x, af in ((b1, "0.3"), (b2, "0.001"), (b3, "0.2"))])
        self.assertEqual(r.base(q), code_of({ref1(q), b1, b3}))  # b2 is the rarest

    def test_record_without_carriers_is_listed(self):
        p = 165
        a1, a2 = others(p)[:2]
        recs = [rec(p, "c1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "c2", ref1(p), a2, "AC=0;AF=0", ["0|0", "0|0", "0|0", "0|0"])]
        self.assertEqual(self.run_enricher(recs).entries(p),
                         ["S1:0|1;%s,%s;c1;0.1" % (ref1(p), a1), ";%s,%s;c2;0" % (ref1(p), a2)])

    def test_non_pass_record_contributes_nothing(self):
        p = 170
        a1, a2, a3 = others(p)
        recs = [rec(p, "p1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "p2", ref1(p), a2, "AC=1;AF=0.9", ["1|1", "1|1", "1|1", "1|1"], filt="LowQual"),
                rec(p, "p3", ref1(p), a3, "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.base(p), code_of({ref1(p), a1, a3}))
        self.assertEqual([e.split(";")[2] for e in r.entries(p)], ["p1", "p3"])
        self.assertNotIn("NOTE:", r.stdout)

    def test_same_alt_on_two_records_keeps_the_later(self):
        p = 180
        a1, a2 = others(p)[:2]
        recs = [rec(p, "d1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "d2", ref1(p), a2, "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0|0"]),
                rec(p, "d3", ref1(p), a1, "AC=1;AF=0.3", ["0|0", "0|0", "1|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), ["S3:1|0;%s,%s;d3;0.3" % (ref1(p), a1),
                                        "S2:0|1;%s,%s;d2;0.2" % (ref1(p), a2)])
        self.assertEqual(r.base(p), code_of({ref1(p), a1, a2}))

    def test_joined_and_split_records_at_one_position(self):
        p = 190
        a1, a2, a3 = others(p)
        recs = [rec(p, "j1", ref1(p), a1 + "," + a2, "AC=1,1;AF=0.1,0.2", ["0|1", "2|0", "0|0", "0|0"]),
                rec(p, "s1", ref1(p), a3, "AC=1;AF=0.05", ["0|0", "0|0", "1|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), ["S1:0|1;%s,%s;j1;0.1" % (ref1(p), a1),
                                        "S2:2|0;%s,%s;j1;0.2" % (ref1(p), a2),
                                        "S3:1|0;%s,%s;s1;0.05" % (ref1(p), a3)])
        self.assertEqual(r.base(p), code_of({ref1(p), a1, a2}))  # a3 is the rarest

    def test_star_record_after_a_snp(self):
        p = 200
        a1 = others(p)[0]
        recs = [rec(p, "t1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "t2", ref1(p), "*", "AC=1;AF=0.2", ["0|0", "1|0", "0|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.base(p), code_of({ref1(p), a1}))
        self.assertEqual(r.entries(p), ["S1:0|1;%s,%s;t1;0.1" % (ref1(p), a1),
                                        "S2:1|0;%s,*;t2;0.2" % ref1(p)])

    def test_non_adjacent_records_merge(self):
        p, q = 210, 215
        a1, a2 = others(p)[:2]
        recs = [rec(p, "n1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(q, "x1", ref1(q), others(q)[0], "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "n2", ref1(p), ref1(p) + "TT", "AC=1;AF=0.1", ["0|0", "0|0", "1|0", "0|0"]),
                rec(p, "n3", ref1(p), a2, "AC=1;AF=0.2", ["0|0", "1|0", "0|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.base(p), code_of({ref1(p), a1, a2}))
        self.assertEqual([e.split(";")[2] for e in r.entries(p)], ["n1", "n3"])

    # ---- interactions with the indel outputs -------------------------------------

    def test_indel_outputs_unmoved_by_split_snps(self):
        """An indel at the same POS and split SNPs inside its window: the log and the
        fake-indel contig are those of the indel alone (cut from the reference)."""
        p = 260
        a1, a2, a3 = others(p)
        q = 270  # inside the 25/27-bp window of the indel at p
        b1, b2 = others(q)[:2]
        indel = [rec(p, "ins", ref1(p), ref1(p) + "GA", "AC=1;AF=0.3", ["1|0", "0|1", "0|0", "0|0"]),
                 rec(p, "del", REF[p - 1:p + 2], ref1(p), "AC=1;AF=0.2", ["0|0", "0|0", "1|0", "0|0"])]
        snps = [rec(p, "s1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "s2", ref1(p), a2, "AC=1;AF=0.2", ["0|0", "1|0", "0|0", "0|0"]),
                rec(p, "s3", ref1(p), a3, "AC=1;AF=0.3", ["0|0", "0|0", "0|1", "0|0"]),
                rec(q, "w1", ref1(q), b1, "AC=1;AF=0.1", ["1|0", "0|0", "0|0", "0|0"]),
                rec(q, "w2", ref1(q), b2, "AC=1;AF=0.1", ["0|0", "0|1", "0|0", "0|0"])]
        alone = self.run_enricher(indel)
        mixed = self.run_enricher([snps[0], indel[0], snps[1], snps[3], indel[1], snps[2], snps[4]])
        log = "SNPs_genome/log%s.txt" % CHROM
        fake = "fake_variants.vcf_%s/fake%s.fa" % (CHROM, CHROM)
        self.assertEqual(mixed.raw(log), alone.raw(log))
        self.assertEqual(mixed.raw(fake), alone.raw(fake))
        self.assertEqual(alone.raw(log).decode().count("\n"), 3)  # header + 2 indels
        self.assertEqual(mixed.base(p), code_of({ref1(p), a2, a3}))  # a1 is the rarest
        self.assertEqual(mixed.base(q), code_of({ref1(q), b1, b2}))
        self.assertEqual([e.split(";")[2] for e in mixed.entries(p)], ["s1", "s2", "s3"])
        self.assertEqual([e.split(";")[2] for e in mixed.entries(q)], ["w1", "w2"])

    # ---- a position with one record behaves as before ----------------------------

    def test_single_record_positions_unchanged(self):
        p, q, s = 300, 310, 320
        a = others(p)[0]
        b1, b2 = others(q)[:2]
        recs = [rec(p, "b", ref1(p), a, "AC=1;AF=0.1", ["0|1", "1|1", "0|0", "0|0"]),
                rec(q, "m", ref1(q), b1 + "," + b2, "AC=1,1;AF=0.1,0.2", ["1|2", "0|2", "0|0", "0|0"]),
                rec(s, "i", ref1(s), ref1(s) + "C", "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.base(p), code_of({ref1(p), a}))
        self.assertEqual(r.base(q), code_of({ref1(q), b1, b2}))
        self.assertEqual(r.base(s), ref1(s))
        self.assertEqual(r.dict, {
            CHROM + ",%d" % p: "S1:0|1,S2:1|1;%s,%s;b;0.1" % (ref1(p), a),
            CHROM + ",%d" % q: "S1:1|2;%s,%s;m;0.1$S1:1|2,S2:0|2;%s,%s;m;0.2"
                               % (ref1(q), b1, ref1(q), b2),
        })
        self.assertNotIn("NOTE:", r.stdout)
        # every other base is the reference
        g = r.genome
        self.assertEqual([i + 1 for i in range(len(REF)) if g[i] != REF[i]], [p, q])

    def test_snp_outside_the_chromosome_fails(self):
        """A SNP past the last base or at POS 0 (a VCF of another assembly): exit 1 with
        a message, instead of an IndexError, a coded last base, or an out-of-bounds write."""
        for pos in (len(REF) + 1, 0):
            r = Run(self.cmd, [rec(pos, "x", "A", "C", "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"])])
            self.addCleanup(r.cleanup)
            self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
            self.assertIn("ERROR: %s:%d is outside the chromosome (%d bp)" % (CHROM, pos, len(REF)),
                          r.stderr)


class PythonPort(EnricherMultiallelicTests, unittest.TestCase):
    cmd = [sys.executable, ENRICHER_PY]


@unittest.skipIf(CPP_BIN is None, "no compiled enricher (set ENRICHER_BIN or install a C++ compiler + zlib)")
class CppPort(EnricherMultiallelicTests, unittest.TestCase):
    cmd = [CPP_BIN] if CPP_BIN else None


def tearDownModule():
    shutil.rmtree(_BUILD_DIR, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
