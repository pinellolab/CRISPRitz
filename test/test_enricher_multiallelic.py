#!/usr/bin/env python3
"""Multi-allelic enrichment: unit tests run against BOTH enricher ports.

A normalised VCF (`bcftools norm -m -any`, how every CRISPRme panel is built) splits
a multi-allelic site into one record per ALT. Up to CRISPRitz 2.9.1 each record
rebuilt the position's IUPAC code from the reference base and overwrote its SNP dict
entry, so only the LAST record's ALT reached the enriched genome and the dict. These
tests pin the merged behaviour on small synthetic inputs, with exact expected outputs:

  * every ALT of a split site is coded, in any record order, and the code equals the
    one the same site gets written as a single record listing every ALT;
  * a position whose reference and ALTs cover all four bases is coded Q, which the
    search reads as any base (N is a gap to the search, not "any base"), and a
    reference N stays N: the enriched genome holds exactly the reference's N;
  * the dict lists every ALT, each entry with its own record's carriers, genotypes,
    rsID and AF, joined by '$' like a multi-ALT record's value; the same REF>ALT on two
    records is one entry with the carriers of both;
  * each record's AF is read by its key, whatever the INFO field's order;
  * non-PASS records, '*', haploid / unphased / missing genotypes, non-adjacent
    records, an indel at the same position, a reference N or IUPAC base, a reference
    that already holds Q (refused);
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

    def test_four_bases_coded_Q(self):
        """A position whose reference base and ALTs cover all four bases is Q, in every
        record order and whatever the ALTs' AFs, and the dict lists all three ALTs."""
        p = 80
        a, b, c = others(p)
        afs = {a: "0.05", b: "0.0001", c: "0.3"}
        split = [rec(p, "rs" + x, ref1(p), x, "AC=1;AF=" + afs[x], ["0|1", "0|0", "0|0", "0|0"])
                 for x in (a, b, c)]
        for order in itertools.permutations(split):
            r = self.run_enricher(list(order))
            self.assertEqual(r.base(p), "Q")
            self.assertEqual(r.stdout.count("NOTE: 1 position(s) of %s carry all four bases" % CHROM), 1)
            self.assertIn("each is coded Q, which the search reads as any base", r.stdout)
            # the dict still lists all three, in file order
            self.assertEqual([e.split(";")[1] for e in r.entries(p)],
                             [ref1(p) + "," + o.split("\t")[4] for o in order])

    def test_four_bases_whatever_the_af_fields(self):
        """Equal AFs, no AF, '.' or an exponent, only per-dataset AF_<x> fields (a
        sites-only merge such as CRISPRme's mega): no ALT is left out, so none ranks."""
        for p, infos in ((90, ("AC=1;AF=0.01",) * 3), (95, ("AC=1",) * 3),
                         (100, ("AC=1;AF=2e-1", "AC=1;AF=.", "AC=1;AF=1.5E-3")),
                         (230, ("AF_one=0.01;AF_max=0.01", "AF_two=0.0002;AF_max=0.0002",
                                "AF_one=0.2;AF_two=0.3;AF_max=0.3"))):
            recs = [rec(p, "r%d%s" % (p, x), ref1(p), x, info, ["0|1", "0|0", "0|0", "0|0"])
                    for x, info in zip(others(p), infos)]
            self.assertEqual(self.run_enricher(recs).base(p), "Q", infos)

    def test_four_bases_joined_record_no_longer_writes_N(self):
        """A record listing three SNP ALTs used to write N, a gap to the search."""
        p = 110
        a, b, c = others(p)
        r = self.run_enricher([rec(p, "rsJ", ref1(p), ",".join([a, b, c]),
                                   "AC=1,1,1;AF=0.2,0.3,0.001", ["0|1", "0|2", "3|0", "0|0"])])
        self.assertEqual(r.base(p), "Q")

    def test_four_bases_over_a_repeated_alt(self):
        """Four records, one ALT twice: Q, and the repeated ALT is one dict entry."""
        p = 115
        a, b, c = others(p)
        recs = [rec(p, "r1", ref1(p), a, "AC=1;AF=0.001", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "r2", ref1(p), b, "AC=1;AF=0.01", ["0|0", "0|1", "0|0", "0|0"]),
                rec(p, "r3", ref1(p), c, "AC=1;AF=0.1", ["0|0", "0|0", "0|1", "0|0"]),
                rec(p, "r4", ref1(p), a, "AC=1;AF=0.5", ["0|0", "0|0", "0|0", "0|1"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.base(p), "Q")
        self.assertEqual(r.entries(p), ["S1:0|1,S4:0|1;%s,%s;r4;0.5" % (ref1(p), a),
                                        "S2:0|1;%s,%s;r2;0.01" % (ref1(p), b),
                                        "S3:0|1;%s,%s;r3;0.1" % (ref1(p), c)])

    def test_reference_holding_Q_is_refused(self):
        """Q is the four-allele symbol, read as any base by the search: a reference FASTA
        holding one (in either case) is refused, with or without a variant there."""
        for k, letter in ((60, "Q"), (250, "q")):
            fasta = REF[:k - 1] + letter + REF[k:]
            r = Run(self.cmd, [rec(50, "x", ref1(50), others(50)[0], "AC=1;AF=0.1",
                                   ["0|1", "0|0", "0|0", "0|0"])], fasta=fasta)
            self.addCleanup(r.cleanup)
            self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
            self.assertIn("ERROR: %s holds 'Q' at %d: the reference FASTA must hold A/C/G/T, N "
                          "and IUPAC codes only (Q is reserved for four-allele SNP positions)"
                          % (CHROM, k), r.stderr)
            self.assertFalse(os.path.exists(r.path("SNPs_genome", "res_enriched")))

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
        # a reference N is a gap and stays one, whatever the ALTs (merged without that rule,
        # these two split records give M; 2.9.1 coded each alone and kept the N); a
        # reference IUPAC base adds no base of its own (the rule of 2.9.1)
        self.assertEqual(r.base(120), "N")
        self.assertEqual(r.base(130), "Y")
        self.assertEqual([e.split(";")[2] for e in r.entries(120)], ["n1", "n2"])
        # one record alone leaves a single base -> the FASTA letter stays, as before
        r1 = self.run_enricher(recs[:1] + recs[2:3], fasta=fasta)
        self.assertEqual(r1.base(120), "N")
        self.assertEqual(r1.base(130), "R")
        # four ALTs on a reference N, split or joined: still N, not Q
        recs4 = [rec(120, "n" + x, "N", x, "AC=1;AF=" + af, ["0|1", "0|0", "0|0", "0|0"])
                 for x, af in (("A", "0.4"), ("C", "0.3"), ("G", "0.01"), ("T", "0.2"))]
        self.assertEqual(self.run_enricher(recs4, fasta=fasta).base(120), "N")
        joined = rec(120, "nJ", "N", "A,C,G,T", "AC=1,1,1,1;AF=0.4,0.3,0.01,0.2",
                     ["0|1", "2|3", "0|4", "0|0"])
        self.assertEqual(self.run_enricher([joined], fasta=fasta).base(120), "N")
        # four ALTs on a reference IUPAC base: they cover every base, so Q
        recs4 = [rec(130, "y" + x, "R", x, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]) for x in "ACGT"]
        self.assertEqual(self.run_enricher(recs4, fasta=fasta).base(130), "Q")

    def test_enriched_genome_keeps_the_references_N(self):
        """The N-count invariant: whatever the records, the enriched genome holds an N
        exactly where the reference does -- a variant never writes one (a record listing
        three SNP ALTs wrote N up to 2.9.1) and never codes one away (a record listing two
        or three ALTs at a reference N was coded up to 2.9.1, and merged split records
        would be)."""
        fasta = list(REF)
        for k in (20, 21, 22, 23, 24, 200, 330):
            fasta[k - 1] = "N"
        fasta = "".join(fasta)
        g1 = ["0|1", "0|0", "0|0", "0|0"]
        recs = []
        for p, alts in ((20, "A"), (21, "A,C"), (22, "A,C,G"), (23, "A,C,G,T"), (200, "G,T"),
                        (330, "N,*")):
            recs.append(rec(p, "j%d" % p, "N", alts, "AC=1;AF=0.1", g1))
        for x in "ACGT":
            recs.append(rec(24, "s" + x, "N", x, "AC=1;AF=0.1", g1))
        for p, n in ((60, 3), (61, 3), (62, 2), (63, 1)):
            for x in others(p)[:n]:
                recs.append(rec(p, "q%d%s" % (p, x), ref1(p), x, "AC=1;AF=0.1", g1))
        recs.append(rec(70, "jq", ref1(70), ",".join(others(70)), "AC=1,1,1;AF=0.1,0.1,0.1", g1))
        recs.append(rec(75, "nn", ref1(75), "N", "AC=1;AF=0.1", g1))
        for mode in ("yes", "no"):
            g = self.run_enricher(recs, mode=mode, fasta=fasta).genome
            self.assertEqual(len(g), len(fasta))
            self.assertEqual([i for i, c in enumerate(g) if c == "N"],
                             [i for i, c in enumerate(fasta) if c == "N"], mode)
            self.assertEqual([i + 1 for i, c in enumerate(g) if c == "Q"], [60, 61, 70], mode)

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
        self.assertEqual(r.base(q), "Q")

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

    def test_same_alt_on_two_records_merges_carriers(self):
        """The same REF>ALT on two records is one entry listing the carriers of both (up
        to 2.9.1 the later record's entry replaced the earlier). A carrier of both phased
        alike gets, slot by slot, the ALT wherever either record puts it (S1 0|1 + 1|0 =
        1|1); otherwise, as unphased slots cannot be matched up, the genotype with more
        copies of it (S4 1 + 1/1 = 1/1; S3 0/1 + 1|. ties: the earlier). The rsID and AF
        are the more frequent record's."""
        p = 180
        a1, a2 = others(p)[:2]
        recs = [rec(p, "d1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0/1", "1"]),
                rec(p, "d2", ref1(p), a2, "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0"]),
                rec(p, "d3", ref1(p), a1, "AC=1;AF=0.3", ["1|0", "0|0", "1|.", "1/1"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), ["S1:1|1,S3:0/1,S4:1/1;%s,%s;d3;0.3" % (ref1(p), a1),
                                        "S2:0|1;%s,%s;d2;0.2" % (ref1(p), a2)])
        self.assertEqual(r.base(p), code_of({ref1(p), a1, a2}))
        # a less frequent second record keeps the first's rsID and AF; a tie too
        for af2, want in (("0.05", ("d1", "0.1")), ("0.1", ("d1", "0.1")), ("1e-1", ("d1", "0.1")),
                          (".", ("d1", "0.1"))):
            recs = [rec(p, "d1", ref1(p), a1, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                    rec(p, "d3", ref1(p), a1, "AC=1;AF=" + af2, ["0|0", "1|0", "0|0", "0|0"])]
            self.assertEqual(self.run_enricher(recs).entries(p),
                             ["S1:0|1,S2:1|0;%s,%s;%s;%s" % ((ref1(p), a1) + want)], af2)

    def test_same_alt_on_a_record_listing_several(self):
        """The same ALT on a split record and as the second ALT of a record listing
        several: the merged entry's genotypes are in split form (1 = that ALT, 0 = the
        reference or another ALT), so a 1|2 carrier of the joined record reads 0|1 there,
        while the joined record's other ALT keeps its genotypes verbatim."""
        p = 185
        a1, a2 = others(p)[:2]
        recs = [rec(p, "s1", ref1(p), a2, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "j1", ref1(p), a1 + "," + a2, "AC=2,2;AF=0.2,0.3", ["1|2", "2|0", "0|1", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), [
            "S1:0|1,S2:1|0;%s,%s;j1;0.3" % (ref1(p), a2),
            "S1:1|2,S3:0|1;%s,%s;j1;0.2" % (ref1(p), a1)])
        # the other way round: the joined record first, then the split record
        r = self.run_enricher(recs[::-1])
        self.assertEqual(r.entries(p), [
            "S1:1|2,S3:0|1;%s,%s;j1;0.2" % (ref1(p), a1),
            "S1:0|1,S2:1|0;%s,%s;j1;0.3" % (ref1(p), a2)])

    def test_af_read_per_record_by_key(self):
        """Each record's INFO/AF, wherever the field stands in its INFO: up to 2.9.1 the AF
        field's place came from the first record, so a record with its INFO fields in
        another order or number got no AF in the dict and the indel log (HPRC chrY lost
        91,950 of 102,665 indel AFs). A key that only starts with AF is not AF."""
        recs = [rec(50, "a", ref1(50), others(50)[0], "AC=1;AN=8", ["0|1", "0|0", "0|0", "0|0"]),
                rec(60, "b", ref1(60), others(60)[0], "AF=0.25;AC=1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(70, "c", ref1(70), others(70)[0], "AC=1;AN=8;AF_afr=0.5;AF=0.125", ["0|1", "0|0", "0|0", "0|0"]),
                rec(80, "d", ref1(80), ref1(80) + "T", "DP=9;AF=0.375;AC=3", ["0|1", "0|0", "0|0", "0|0"]),
                rec(90, "e", ref1(90), others(90)[0] + "," + others(90)[1], "AN=8;AC=1,2;AF=0.1,0.2",
                    ["0|1", "0|2", "0|0", "0|0"]),
                rec(100, "f", ref1(100), others(100)[0], "AF_max=0.4;AC=1", ["0|1", "0|0", "0|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual([r.entries(p)[0].split(";")[3] for p in (50, 60, 70, 90, 100)],
                         ["", "0.25", "0.125", "0.1", ""])
        self.assertEqual(r.entries(90)[1].split(";")[3], "0.2")
        log = r.raw("SNPs_genome/log%s.txt" % CHROM).decode().splitlines()
        self.assertEqual([row.split("\t")[3] for row in log[1:]], ["0.375"])

    def test_joined_and_split_records_at_one_position(self):
        p = 190
        a1, a2, a3 = others(p)
        recs = [rec(p, "j1", ref1(p), a1 + "," + a2, "AC=1,1;AF=0.1,0.2", ["0|1", "2|0", "0|0", "0|0"]),
                rec(p, "s1", ref1(p), a3, "AC=1;AF=0.05", ["0|0", "0|0", "1|0", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), ["S1:0|1;%s,%s;j1;0.1" % (ref1(p), a1),
                                        "S2:2|0;%s,%s;j1;0.2" % (ref1(p), a2),
                                        "S3:1|0;%s,%s;s1;0.05" % (ref1(p), a3)])
        self.assertEqual(r.base(p), "Q")

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
        self.assertEqual(mixed.base(p), "Q")
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
