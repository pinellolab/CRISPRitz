#!/usr/bin/env python3
"""Multi-allelic enrichment: unit tests run against BOTH enricher ports.

A normalised VCF (`bcftools norm -m -any`, how every CRISPRme panel is built) splits
a multi-allelic site into one record per ALT. Up to CRISPRitz 2.9.1 each record
rebuilt the position's IUPAC code from the reference base and overwrote its SNP dict
entry, so only the LAST record's ALT reached the enriched genome and the dict. These
tests pin the merged behaviour on small synthetic inputs, with exact expected outputs:

  * every ALT of a split site is coded, in any record order, and the code equals the
    one the same site gets written as a single record listing every ALT;
  * a position whose reference base (A, C, G or T) and ALTs cover all four bases is
    coded Q, which the search reads as any base (N is a gap to the search, not "any base"), and a
    reference N stays N: the enriched genome holds exactly the reference's N;
  * the dict lists every ALT, each entry with its own record's carriers, genotypes,
    rsID and AF, joined by '$' like a multi-ALT record's value; the same REF>ALT on two
    records is one entry with the carriers of both; every entry's genotypes are in split
    form (1 its own ALT), a record listing several ALTs' too;
  * each record's AF is read by its key, whatever the INFO field's order;
  * non-PASS records, '*', haploid / unphased / missing genotypes, non-adjacent
    records, an indel at the same position, a reference N, IUPAC base or byte that is
    neither (kept, as an N), a reference that already holds Q (refused);
  * the indel log and fake-indel contigs do not move; each contig is 30 reference bases,
    the REF with the ALT in its place, and 30 reference bases (the log's CHR field and
    FAKEPOS say where);
  * a position with one record of one ALT behaves exactly as before;
  * REF and ALTs in lower or mixed case are read upper-cased, as VCF 4.x reads them (a
    symbolic ALT as written): a lower-case SNP ALT is coded, a record differing from
    another only by case is one allele with it, an indel's window and log row are in
    upper case, and every output is that of the VCF written in upper case;
  * the dictionary examples README.md, docs/MULTIALLELIC_ENRICHMENT.md and CHANGELOG.md
    give are what the port writes.

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
import re
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


def vcf_header(samples=SAMPLES):
    return [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]


def recased(line, case):
    """VCF data line ``line`` with its REF and every ALT passed through ``case`` (a
    str -> str function), a symbolic ALT (holding '<') as written."""
    f = line.split("\t")
    f[3] = case(f[3])
    f[4] = ",".join(a if "<" in a else case(a) for a in f[4].split(","))
    return "\t".join(f)


def alternate(s):
    """s with every other letter in lower case, from the first: 'ACGT' -> 'aCgT'."""
    return "".join(c.lower() if i % 2 == 0 else c for i, c in enumerate(s))


# The documents whose dictionary examples test_the_documented_dictionary_examples runs.
DOC_FILES = ("README.md", os.path.join("docs", "MULTIALLELIC_ENRICHMENT.md"), "CHANGELOG.md")
_DOC_GT = r"[0-9.]+(?:[|/][0-9.]+)*"
# "a `1|2` carrier of `C>A,T` reads `1|0` under A and `0|1` under T"
_DOC_SENTENCE = re.compile(r"a `(%s)` carrier of `([ACGT])>([ACGT](?:,[ACGT])+)` reads `(%s)` "
                           r"under ([ACGT]) and `(%s)` under ([ACGT])\b" % (_DOC_GT, _DOC_GT, _DOC_GT))


def doc_sentence_examples(text):
    """Every 'a `<GT>` carrier of `<REF>><ALT>,<ALT>...` reads `<GT>` under <ALT> and `<GT>`
    under <ALT>' of a Markdown text (its line breaks read as spaces), as tuples (GT, REF,
    ALTs, GT under the first ALT named, that ALT, GT under the second, that ALT)."""
    return _DOC_SENTENCE.findall(" ".join(text.split()))


def doc_worked_example(text):
    """The worked example of docs/MULTIALLELIC_ENRICHMENT.md section 2.3: the fenced block
    that starts with a VCF header line (#CHROM ...; fields split on blanks), and the fenced
    block after it, which lists the "<chr>,<pos>": "<value>" pairs those records give in
    my_dict. Returns (rows of fields, header first; {key: value}), or None."""
    blocks = re.findall(r"^```[a-z]*\n(.*?)^```", text, flags=re.S | re.M)
    for i, block in enumerate(blocks[:-1]):
        if block.startswith("#CHROM"):
            rows = [line.split() for line in block.strip().splitlines()]
            pairs = json.loads("{" + ",".join(blocks[i + 1].strip().splitlines()) + "}")
            return rows, pairs
    return None


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
        lines = header or vcf_header()
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

    # the four outputs of a run in mode yes (mode no writes the enriched genome alone)
    OUTPUTS = ("SNPs_genome/res_enriched/%s.enriched.fa" % CHROM,
               "SNPs_genome/my_dict_%s.json" % CHROM,
               "SNPs_genome/log%s.txt" % CHROM,
               "fake_variants.vcf_%s/fake%s.fa" % (CHROM, CHROM))

    def assert_same_outputs(self, got, want, what):
        """Every output ``want`` wrote, byte for byte in ``got``, and stdout but its
        DONE line (which holds the run's seconds)."""
        for rel in self.OUTPUTS[:1 if want.mode == "no" else 4]:
            self.assertEqual(got.raw(rel), want.raw(rel), "%s: %s" % (what, rel))
        done = re.compile(r"^DONE IN .*$", re.M)
        self.assertEqual(done.sub("", got.stdout), done.sub("", want.stdout), what)

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
        """A position whose reference base (A, C, G or T) and ALTs cover all four bases is
        Q, in every record order and whatever the ALTs' AFs, and the dict lists all three
        ALTs."""
        p = 80
        a, b, c = others(p)
        afs = {a: "0.05", b: "0.0001", c: "0.3"}
        split = [rec(p, "rs" + x, ref1(p), x, "AC=1;AF=" + afs[x], ["0|1", "0|0", "0|0", "0|0"])
                 for x in (a, b, c)]
        for order in itertools.permutations(split):
            r = self.run_enricher(list(order))
            self.assertEqual(r.base(p), "Q")
            self.assertEqual(r.stdout.count("NOTE: 1 position(s) of %s are coded Q, the four-allele "
                                            "symbol, which the search reads as any base" % CHROM), 1)
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
            self.assertIn("ERROR: %s holds 'Q' at %d: a reference FASTA must not hold Q, which is "
                          "reserved for four-allele SNP positions" % (CHROM, k), r.stderr)
            self.assertFalse(os.path.exists(r.path("SNPs_genome", "res_enriched")))

    def test_reference_N_and_reference_iupac(self):
        fasta = list(REF)
        fasta[119] = "N"   # pos 120
        fasta[129] = "R"   # pos 130 (A/G in the FASTA itself)
        fasta[139] = "X"   # pos 140: no base, IUPAC code or N (RepeatMasker -x masks with X)
        fasta[149] = "-"   # pos 150
        fasta = "".join(fasta)
        recs = [rec(120, "n1", "N", "A", "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(120, "n2", "N", "C", "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0|0"]),
                rec(130, "y1", "R", "C", "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(130, "y2", "R", "T", "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0|0"]),
                rec(140, "x1", ref1(140), others(140)[0], "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(140, "x2", ref1(140), others(140)[1], "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0|0"]),
                rec(150, "d1", ref1(150), ",".join(others(150)[:2]), "AC=1,1;AF=0.1,0.2",
                    ["0|1", "2|0", "0|0", "0|0"])]
        r = self.run_enricher(recs, fasta=fasta)
        # a reference N is a gap and stays one, whatever the ALTs (merged without that rule,
        # these two split records give M; 2.9.1 coded each alone and kept the N); a
        # reference IUPAC base adds no base of its own (the rule of 2.9.1)
        self.assertEqual(r.base(120), "N")
        self.assertEqual(r.base(130), "Y")
        self.assertEqual([e.split(";")[2] for e in r.entries(120)], ["n1", "n2"])
        # so is a byte that is no base or IUPAC code, which buildTST and searchBruteForce
        # read as N: merged, the split records at 140 would code it, and up to 2.9.1 the
        # record listing two ALTs at 150 did
        self.assertEqual(r.base(140), "X")
        self.assertEqual(r.base(150), "-")
        self.assertEqual([e.split(";")[2] for e in r.entries(140)], ["x1", "x2"])
        self.assertEqual([e.split(";")[2] for e in r.entries(150)], ["d1", "d1"])
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
        # one record listing all four bases, which 2.9.1 coded N at both: Q over the
        # reference IUPAC base, and the byte that is no base or IUPAC code stays
        gt4 = ["0|1", "2|3", "0|4", "0|0"]
        joined4 = [rec(130, "yJ", "R", "A,C,G,T", "AC=1,1,1,1;AF=0.1,0.1,0.1,0.1", gt4),
                   rec(140, "xJ", ref1(140), "A,C,G,T", "AC=1,1,1,1;AF=0.1,0.1,0.1,0.1", gt4)]
        r4 = self.run_enricher(joined4, fasta=fasta)
        self.assertEqual([r4.base(130), r4.base(140)], ["Q", "X"])

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
        several: every entry's genotypes are in split form (1 = that ALT, 0 = the
        reference or another ALT), the merged one's and the joined record's other ALT's
        too, so a 1|2 carrier of the joined record reads 0|1 under its second ALT and 1|0
        under its first."""
        p = 185
        a1, a2 = others(p)[:2]
        recs = [rec(p, "s1", ref1(p), a2, "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                rec(p, "j1", ref1(p), a1 + "," + a2, "AC=2,2;AF=0.2,0.3", ["1|2", "2|0", "0|1", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), [
            "S1:0|1,S2:1|0;%s,%s;j1;0.3" % (ref1(p), a2),
            "S1:1|0,S3:0|1;%s,%s;j1;0.2" % (ref1(p), a1)])
        # the other way round: the joined record first, then the split record
        r = self.run_enricher(recs[::-1])
        self.assertEqual(r.entries(p), [
            "S1:1|0,S3:0|1;%s,%s;j1;0.2" % (ref1(p), a1),
            "S1:0|1,S2:1|0;%s,%s;j1;0.3" % (ref1(p), a2)])

    def test_a_record_listing_several_beside_a_split_record_numbers_its_own_alts(self):
        """C>G (a split record, S1 1|0) and then C>A,CT (S2 1|2: A on the first slot, the
        insertion on the second): A's entry holds S2's genotype in split form, 1|0, so
        its slot is told whatever else the position lists. Kept raw (1|2), CRISPRme's
        readers tied A to allele 2 against the split record's 1 and read A on the
        insertion's slot (the CRISPRme review's world X, site C); an unphased 1/2 stays
        unphased (0/1 under its second ALT), a haploid 2 is 1 under that ALT."""
        p = 187
        a1, a2 = others(p)[:2]
        recs = [rec(p, "s1", ref1(p), a2, "AC=1;AF=0.1", ["1|0", "0|0", "0|0", "0|0"]),
                rec(p, "j1", ref1(p), a1 + "," + ref1(p) + "T", "AC=1,1;AF=0.2,0.3",
                    ["0|0", "1|2", "2/1", "0|0"])]
        r = self.run_enricher(recs)
        self.assertEqual(r.entries(p), [
            "S1:1|0;%s,%s;s1;0.1" % (ref1(p), a2),
            "S2:1|0,S3:0/1;%s,%s;j1;0.2" % (ref1(p), a1)])
        q = 189
        b1, b2, b3 = others(q)
        recs = [rec(q, "m", ref1(q), b1 + "," + b2 + "," + b3, "AC=1,1,1;AF=0.1,0.2,0.3",
                    ["1|2", "3|1", "2", "0|3"])]
        self.assertEqual(self.run_enricher(recs).entries(q), [
            "S1:1|0,S2:0|1;%s,%s;m;0.1" % (ref1(q), b1),
            "S1:0|1,S3:1;%s,%s;m;0.2" % (ref1(q), b2),
            "S2:1|0,S4:0|1;%s,%s;m;0.3" % (ref1(q), b3)])

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
        # S2's 2|0 in split form under the joined record's second ALT (2.9.1: raw 2|0)
        self.assertEqual(r.entries(p), ["S1:0|1;%s,%s;j1;0.1" % (ref1(p), a1),
                                        "S2:1|0;%s,%s;j1;0.2" % (ref1(p), a2),
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
        q = 270  # inside the window of the indel at p (30 bases on each side)
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

    def test_indel_window_is_30_bases_on_each_side(self):
        """Each indel allele's contig is 30 reference bases, the REF (from POS on) with the
        ALT in its place, then 30 reference bases (fewer where the chromosome ends), one N
        between contigs: the INDEL_LEFT and INDEL_RIGHT crispritz.py states. The log's CHR
        field gives the window's real start and end (POS - 31 and POS - 1 + len(REF) + 30,
        0-based, the end unclipped), its FAKEPOS the contig on the fake chromosome, and
        refseq the reference from the window's start, as long as the contig. An indel whose
        window would start before the chromosome is skipped with a warning on stdout (POS
        30), while one at POS 31 starts at base 0. Up to 2.9.1 the window was 25 bases
        before the REF and 27 after it, and crispritz.py stated none."""
        with open(os.path.join(REPO, "crispritz.py")) as f:
            stated = dict(re.findall(r"^(INDEL_LEFT|INDEL_RIGHT) = (\d+)$", f.read(), re.M))
        self.assertEqual(stated, {"INDEL_LEFT": "30", "INDEL_RIGHT": "30"})
        end = len(REF)
        g = ["0|1", "0|0", "0|0", "0|0"]
        mnv_alt = others(220)[0] + others(221)[0]
        recs = [  # (POS, REF, ALT) as each record lists them
            (30, ref1(30), ref1(30) + "T"),             # skipped: its window starts at -1
            (31, ref1(31), ref1(31) + "GG"),            # its window starts at base 0
            (100, ref1(100), ref1(100) + "GA"),         # an insertion
            (160, REF[159:162], ref1(160)),             # a 3-base deletion
            (220, REF[219:221], mnv_alt),               # an MNV: REF and ALT differ at once
            (end - 10, REF[end - 11:end - 8], ref1(end - 10)),  # 8 bases left after its REF
        ]
        r = self.run_enricher([rec(p, "i%d" % p, ref, alt, "AC=1;AF=0.1", g)
                               for p, ref, alt in recs])
        self.assertIn("WARNING: skipping indel at %s:30 -- fewer than 31 bases from the "
                      "chromosome start" % CHROM, r.stdout)
        log = r.raw("SNPs_genome/log%s.txt" % CHROM).decode().splitlines()
        self.assertEqual(log[0], "CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq")
        rows = [x.split("\t") for x in log[1:]]
        flat, want_rows = "", []
        for k, (p, ref, alt) in enumerate(recs[1:]):
            start = p - 1 - 30
            contig = REF[start:p - 1] + alt + REF[p - 1 + len(ref):p - 1 + len(ref) + 30]
            want_rows.append([
                "%s_%d-%d_%d" % (CHROM, start, p - 1 + len(ref) + 30, k + 1), "S1", "i%d" % p,
                "0.1", "%s_%d_%s_%s" % (CHROM, p, ref, alt),
                "%d,%d" % (len(flat), len(flat) + len(contig)), REF[start:start + len(contig)]])
            flat += contig + "N"
        self.assertEqual(rows, want_rows)
        # the left flank, read back from the log as a reader derives it: POS - 1 - start
        self.assertEqual({int(x[4].split("_")[1]) - 1 - int(x[0].rsplit("_", 1)[0].rsplit("_", 1)[1]
                                                            .split("-")[0]) for x in rows}, {30})
        fake = r.raw("fake_variants.vcf_%s/fake%s.fa" % (CHROM, CHROM)).decode()
        self.assertEqual(fake, ">fake%s\n" % CHROM
                         + "".join(c + "\nN\n" for c in flat.split("N")[:-1]))
        # the right flank of the last contig is cut by the chromosome's end: 8 bases
        last = flat.split("N")[-2]
        self.assertEqual(len(last), 30 + 1 + 8)

    # ---- REF and ALT in lower or mixed case ----------------------------------------
    # VCF 4.x reads the bases of REF and ALT whatever their case, and so does the enricher
    # since 2.10.0: it upper-cases REF and every ALT as it reads a record (a symbolic ALT,
    # which holds '<', as written; '*' and '.' hold no letter). Up to 2.9.1 a lower-case
    # SNP ALT was not coded, a record differing from another only by case was another
    # allele, and an indel's window and log row held the VCF's case.

    def test_an_snv_alt_in_lower_case_is_coded(self):
        """A SNP ALT in lower case (C>g) is coded as its upper case, in both modes, and
        its dict entry names the record in upper case, as one whose REF is in lower case
        (c>G); three ALTs in lower case over an A, C, G or T base are Q. Up to 2.9.1 the
        genome kept the reference base where the ALT was in lower case, so the search
        never read the allele, while the dict listed its carriers (CRISPRme's world WLOW:
        a perfect target on a carrier's haplotype, never found)."""
        p, q, s = 40, 45, 55
        a, b = others(p)[0], others(q)[1]
        recs = [rec(p, "l1", ref1(p), a.lower(), "AC=3;AF=0.375", ["0|1", "1|1", "0|0", "0|0"]),
                rec(q, "l2", ref1(q).lower(), b, "AC=1;AF=0.125", ["0|0", "1|0", "0|0", "0|0"])]
        recs += [rec(s, "l3" + x, ref1(s).lower(), x.lower(), "AC=1;AF=0.1",
                     ["0|1", "0|0", "0|0", "0|0"]) for x in others(s)]
        r = self.run_enricher(recs)
        self.assertEqual(r.base(p), code_of({ref1(p), a}))
        self.assertEqual(r.base(q), code_of({ref1(q), b}))
        self.assertEqual(r.base(s), "Q")
        self.assertEqual(r.entries(p), ["S1:0|1,S2:1|1;%s,%s;l1;0.375" % (ref1(p), a)])
        self.assertEqual(r.entries(q), ["S2:1|0;%s,%s;l2;0.125" % (ref1(q), b)])
        self.assertEqual([e.split(";")[1] for e in r.entries(s)],
                         ["%s,%s" % (ref1(s), x) for x in others(s)])
        self.assertEqual(r.stdout.count("NOTE: 1 position(s) of %s are coded Q" % CHROM), 1)
        self.assertEqual(self.run_enricher(recs, mode="no").genome, r.genome)
        self.assert_same_outputs(r, self.run_enricher([recased(x, str.upper) for x in recs]),
                                 "the records in upper case")

    def test_indels_in_lower_or_mixed_case(self):
        """An indel or MNV whose REF or ALT is in lower or mixed case (bcftools norm -f
        writes the case of a soft-masked reference into an indel it realigns) is the
        record in upper case: its window holds the ALT upper-cased, its log row names it
        so, and the log and the fake-indel FASTA are byte for byte those of the records in
        upper case. Up to 2.9.1 the window held the ALT as written (buildTST reads a
        FASTA upper-cased, so its index was the same) and the log the VCF's text."""
        g = ["0|1", "1|0", "0|0", "0|0"]
        mnv = others(250)[0] + others(251)[0]
        recs = [  # (POS, REF, ALT) as the records write them
            (100, alternate(REF[99:103]), REF[99].lower()),        # a deletion, mixed REF
            (150, ref1(150).lower(), (ref1(150) + "ACG").lower()),  # an insertion in lower case
            (200, ref1(200), ref1(200) + "gTa"),                    # an insertion, mixed ALT
            (250, REF[249:251].lower(), alternate(mnv)),            # an MNV
            (300, REF[299:310].lower(), REF[299].lower()),          # a deletion of 10 bases
        ]
        lines = [rec(p, "i%d" % p, ref, alt, "AC=2;AF=0.2", g) for p, ref, alt in recs]
        r = self.run_enricher(lines)
        log = r.raw("SNPs_genome/log%s.txt" % CHROM).decode().splitlines()
        self.assertEqual([x.split("\t")[4] for x in log[1:]],
                         ["%s_%d_%s_%s" % (CHROM, p, ref.upper(), alt.upper()) for p, ref, alt in recs])
        fake = r.raw("fake_variants.vcf_%s/fake%s.fa" % (CHROM, CHROM)).decode()
        self.assertEqual(fake, ">fake%s\n" % CHROM + "".join(
            REF[p - 31:p - 1] + alt.upper() + REF[p - 1 + len(ref):p - 1 + len(ref) + 30] + "\nN\n"
            for p, ref, alt in recs))
        self.assert_same_outputs(r, self.run_enricher([recased(x, str.upper) for x in lines]),
                                 "the records in upper case")

    def test_lower_case_in_a_record_listing_several_alts(self):
        """c>a,G,ctt: each ALT is read in upper case, so A and G are coded, both dict
        entries name the record in upper case (genotypes in split form), and the
        insertion's window and log row hold CTT. Up to 2.9.1 the lower-case a was not
        coded and the entries kept the VCF's case."""
        p = 190
        a1, a2 = others(p)[:2]
        line = rec(p, "j", ref1(p).lower(), ",".join([a1.lower(), a2, (ref1(p) + "TT").lower()]),
                   "AC=1,2,1;AF=0.1,0.2,0.3", ["1|2", "0|3", "2|0", "0|0"])
        r = self.run_enricher([line])
        self.assertEqual(r.base(p), code_of({ref1(p), a1, a2}))
        self.assertEqual(r.entries(p), ["S1:1|0;%s,%s;j;0.1" % (ref1(p), a1),
                                        "S1:0|1,S3:1|0;%s,%s;j;0.2" % (ref1(p), a2)])
        log = r.raw("SNPs_genome/log%s.txt" % CHROM).decode().splitlines()
        self.assertEqual([x.split("\t")[1:5] for x in log[1:]],
                         [["S2", "j", "0.3", "%s_%d_%s_%sTT" % (CHROM, p, ref1(p), ref1(p))]])
        self.assert_same_outputs(r, self.run_enricher([recased(line, str.upper)]),
                                 "the record in upper case")

    def test_a_record_differing_only_by_case_is_the_same_allele(self):
        """The same REF>ALT on two records is one allele (CRISPRme decision B2), whatever
        their case: C>G, c>g and C>g at one position are one dict entry with the carriers
        of all three, genotypes merged and the rsID and AF of the most frequent, exactly
        as three records written alike; two indel records differing by case get a window
        and a log row each, in upper case, exactly as two written alike. Up to 2.9.1 each
        spelling was an allele of its own (an entry each, the lower-case ALTs not coded)."""
        p, q = 180, 260
        a = others(p)[0]
        lines = [rec(p, "d1", ref1(p), a, "AC=1;AF=0.1", ["0|1", "0|0", "0/1", "1"]),
                 rec(p, "d2", ref1(p).lower(), a.lower(), "AC=1;AF=0.3", ["1|0", "0|0", "1|.", "1/1"]),
                 rec(p, "d3", ref1(p), a.lower(), "AC=1;AF=0.2", ["0|0", "0|1", "0|0", "0"]),
                 rec(q, "e1", REF[q - 1:q + 1], ref1(q), "AC=1;AF=0.1", ["0|1", "0|0", "0|0", "0|0"]),
                 rec(q, "e2", REF[q - 1:q + 1].lower(), ref1(q).lower(), "AC=1;AF=0.2",
                     ["0|0", "1|0", "0|0", "0|0"])]
        r = self.run_enricher(lines)
        self.assertEqual(r.entries(p), ["S1:1|1,S2:0|1,S3:0/1,S4:1/1;%s,%s;d2;0.3" % (ref1(p), a)])
        self.assertEqual(r.base(p), code_of({ref1(p), a}))
        log = r.raw("SNPs_genome/log%s.txt" % CHROM).decode().splitlines()
        self.assertEqual([x.split("\t")[1:5] for x in log[1:]],
                         [[s, e, af, "%s_%d_%s_%s" % (CHROM, q, REF[q - 1:q + 1], ref1(q))]
                          for s, e, af in (("S1", "e1", "0.1"), ("S2", "e2", "0.2"))])
        self.assert_same_outputs(r, self.run_enricher([recased(x, str.upper) for x in lines]),
                                 "the records written alike")

    def test_every_output_is_that_of_the_records_in_upper_case(self):
        """A VCF of SNPs (split, joined, with '*', 'N' or '.' among the ALTs, at a
        four-base site, at a reference N, outside FILTER PASS), indels and MNVs (joined
        with SNPs too) and symbolic ALTs, each REF and ALT in a random mix of cases (a
        symbolic ALT as written): the enriched genome, the dict, the indel log, the
        fake-indel FASTA and stdout are byte for byte those of the same VCF in upper case,
        in both modes."""
        fasta = REF[:69] + "N" + REF[70:]
        g1, g2 = ["0|1", "0|0", "0|0", "0|0"], ["1|2", "0|1", "2/2", "0|0"]
        o = others
        base = [
            rec(40, "s1", ref1(40), o(40)[0], "AC=1;AF=0.1", g1),
            rec(40, "s2", ref1(40), o(40)[1], "AC=1;AF=0.2", ["0|0", "1|0", "0|0", "0|0"]),
            rec(50, "j1", ref1(50), o(50)[0] + "," + o(50)[1], "AC=1,1;AF=0.1,0.2", g2),
            rec(60, "j2", ref1(60), o(60)[2] + ",*", "AC=1,1;AF=0.1,0.2", g2),
            rec(65, "j3", ref1(65), "N," + o(65)[0], "AC=1,1;AF=0.1,0.2", g2),
            rec(66, "j4", ref1(66), ".", "AC=0;AF=0", g1),
            rec(70, "n1", "N", "A,C", "AC=1,1;AF=0.1,0.2", g2),
            rec(80, "q1", ref1(80), ",".join(o(80)), "AC=1,1,1;AF=0.1,0.2,0.3", ["1|2", "3|0", "0|0", "0|0"]),
            rec(90, "f1", ref1(90), o(90)[0], "AC=1;AF=0.1", g1, filt="LowQual"),
            rec(100, "y1", ref1(100), o(100)[0] + ",<DEL>", "AC=1,1;AF=0.1,0.2", g2),
            rec(105, "y2", ref1(105), "<del>," + o(105)[1], "AC=1,1;AF=0.1,0.2", g2),
            rec(120, "i1", ref1(120), ref1(120) + "GAT", "AC=1;AF=0.1", g1),
            rec(130, "d1", REF[129:136], ref1(130), "AC=1;AF=0.1", ["0|0", "1|1", "0|0", "0|0"]),
            rec(140, "m1", REF[139:142], o(140)[0] + o(141)[0] + o(142)[0], "AC=1;AF=0.1", g1),
            rec(150, "x1", ref1(150), o(150)[0] + "," + ref1(150) + "CC", "AC=1,1;AF=0.1,0.2", g2),
            rec(160, "x2", REF[159:162], REF[159] + "," + REF[159:162] + "A", "AC=1,1;AF=0.1,0.2", g2),
            rec(170, "t1", ref1(170), o(170)[0], "AC=1;AF=0.1", g1),
            rec(170, "t2", ref1(170), ref1(170) + "TTT", "AC=1;AF=0.1", ["0|0", "1|0", "0|0", "0|0"]),
        ]
        rnd = random.Random(7)
        cased = [recased(x, lambda s: "".join(c.lower() if rnd.random() < 0.5 else c for c in s))
                 for x in base]
        self.assertGreater(sum(x != y for x, y in zip(cased, base)), 12)
        for mode in ("yes", "no"):
            self.assert_same_outputs(self.run_enricher(cased, mode=mode, fasta=fasta),
                                     self.run_enricher(base, mode=mode, fasta=fasta),
                                     "mode %s, the VCF in upper case" % mode)

    # ---- the dictionary the documents describe -----------------------------------

    def test_the_documented_dictionary_examples(self):
        """What the documents say the dictionary holds is what this port writes: the worked
        example of docs/MULTIALLELIC_ENRICHMENT.md section 2.3 (VCF records and the my_dict
        values they give, compared whole), and each "a `<GT>` carrier of `<REF>><ALTs>`
        reads `<GT>` under <ALT> and `<GT>` under <ALT>" of README.md, that document and
        CHANGELOG.md, each of which must give one. Up to 1a9812d the README said every
        entry kept its record's genotypes verbatim and section 2.3 that a key holding one
        record was byte-identical to 2.9.1, both false for a record listing several ALTs
        since the split form, and neither gave an example a test could run."""
        problems = []
        runs = {}
        for name in DOC_FILES:
            with open(os.path.join(REPO, name)) as f:
                text = f.read()
            examples = doc_sentence_examples(text)
            if not examples:
                problems.append("%s gives no \"a `<GT>` carrier of `<REF>><ALTs>` reads `<GT>` "
                                "under <ALT> and `<GT>` under <ALT>\"" % name)
            for gt, ref, alts, g1, a1, g2, a2 in examples:
                key = (gt, ref, alts)
                if key not in runs:
                    p = 120
                    r = self.run_enricher(
                        [rec(p, "ex", ref, alts, "AF=" + ",".join(["0.1"] * len(alts.split(","))),
                             [gt, "0|0", "0|0", "0|0"])],
                        fasta=REF[:p - 1] + ref + REF[p:])
                    runs[key] = {}
                    for e in r.entries(p):  # "<sample:GT,...>;<REF>,<ALT>;<rsID>;<AF>"
                        carriers, ref_alt = e.split(";")[:2]
                        runs[key][ref_alt] = dict(c.rsplit(":", 1) for c in carriers.split(",") if c)
                for want, alt in ((g1, a1), (g2, a2)):
                    got = runs[key].get(ref + "," + alt, {}).get("S1")
                    if got != want:
                        problems.append("%s: a `%s` carrier of `%s>%s` reads %r under %s, not `%s`"
                                        % (name, gt, ref, alts, got, alt, want))
            if name.endswith("MULTIALLELIC_ENRICHMENT.md"):
                worked = doc_worked_example(text)
                if worked is None:
                    problems.append("%s gives no worked example (a fenced VCF block, then the "
                                    "my_dict values it gives)" % name)
                    continue
                rows, want = worked
                fasta = list(REF)
                records = []
                for fields in rows[1:]:
                    fasta[int(fields[1]) - 1] = fields[3]
                    records.append("\t".join([CHROM] + fields[1:]))
                r = self.run_enricher(records, fasta="".join(fasta), header=vcf_header(rows[0][9:]))
                got = dict((rows[1][0] + k[len(CHROM):], v) for k, v in r.dict.items())
                if got != want:
                    problems.append("%s: the worked example gives %s, not %s"
                                    % (name, json.dumps(got, sort_keys=True), json.dumps(want, sort_keys=True)))
        self.assertFalse(problems, "\n" + "\n".join(problems))

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
            # a record listing several ALTs: each entry's genotypes in split form (2.9.1
            # wrote the raw S1:1|2 under both)
            CHROM + ",%d" % q: "S1:1|0;%s,%s;m;0.1$S1:0|1,S2:0|1;%s,%s;m;0.2"
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
