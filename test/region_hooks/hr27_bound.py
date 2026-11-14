#!/usr/bin/env python3
"""The rule of docs/INDEX_FORMAT.md section 13.5 for a scan of the dense fake-indel windows
alone (review findings H-R2-7 and HK-R1-1), checked on planted genomes.

usage: hr27_bound.py genome <out_dir>
       hr27_bound.py check <out_dir> <tag> <work_dir> <G> <L> <mb> <bRNA> <Gs>...
       hr27_bound.py check-pam <out_dir> <tag> <work_dir> <G> <L> <mb>

A scan of the dense windows alone misses the lines a skip BED of them drops that cover the
ALT of a window that is not dense. The rule: with c = 1 where the interval of the dense
window a core reaches holds the separator N between it and the window it borders, and
c = 2 where it does not, none is missed while LEFT >= m - c and RIGHT >= m - c (for the
window before and the window after), m = n = Gs + L on the pad route (case 1) and
m = max(Gs, G + mb) + L on a PAM-less index searched without --prune-pam (case 2).

genome writes fake-indel genomes laid out as add-variants writes one: windows of LEFT
reference bases, an ALT (here 3 bases) and RIGHT reference bases, one N between two
windows, each window that is not dense holding one planted site between two dense ones
(random bases); <tag>/fakechrA.fa, <tag>.windows.tsv (contig, start, end, left, right,
dense, site), and the dense windows as a BED in the four forms a consumer may write them:
<tag>.win.bed [start, end); <tag>.trail.bed [start, end + 1), with the separator after
the window, as the haplotype-scan design's INDELS BED; <tag>.lead.bed [start - 1, end);
<tag>.both.bed [start - 1, end + 1).
  case 1: hr25 (LEFT/RIGHT 25/27, CRISPRitz up to 2.9.1) and hr30 (30/30, 2.10.0); an NNGRRT
     site of the guide region GR, PAM-right ('R<o>', + strand, its guide region starting o
     bases into the window, the o bases before it A) or PAM-left ('L<o>', - strand, its
     guide region ending o bases before the window's end, the o bases after it T), o = 0..5;
     guides_hr<Gs>.txt, 'A' * (Gs - 20) + GR, for Gs = 21..26.
  case 2: st28, st29, st30 (LEFT = RIGHT = 28, 29, 30); GR next to the window's edge on its
     PAM side, PAM-right ('R<j>', + strand, its guide region ending j bases before the
     window's end) or PAM-left ('L<j>', - strand, starting j bases into the window), j =
     0..7, with GR[:3] on its far side (so that 3 DNA bulges can consume the pad with no
     mismatch); guides_st.txt, GR. Searched on a PAM-less K=8 index without --prune-pam,
     whose flank is not tested: a line's PAM then runs across the separator.

check (case 1) reads <work_dir>/<tag>.<Gs>.U.targets.txt (the search without a BED) and
<tag>.<Gs>.<form>.targets.txt (with the skip BED of that form). For each Gs and form the
sites that lose a line lying wholly in their window and covering its ALT must be exactly
the rule's, site by site: a site at offset o loses lines iff o <= Gs - G - c and
o + G + L > LEFT (PAM-right; c of the window before) or RIGHT (PAM-left; the window after)
-- at the worst offset iff n >= LEFT + 1 + c, resp. RIGHT + 1 + c -- provided it has lines
at all: Gs - G - min(o, mb) <= bRNA for a search with no DNA bulge (the guide's extra bases
that find no pad base, the pad ending at the separator or after mb bases, are RNA bulges).
Every site that has lines must have them in U, no other site, and no line may run across a
separator (the index is a dedicated one: every PAM letter is tested, and N matches none).

check-pam (case 2) reads <tag>.U.targets.txt and <tag>.<form>.targets.txt: some lines must
run across a separator, and the sides that lose a line covering the ALT of a window that is
not dense must be the rule's: PAM-right lines (into the window after) iff
RIGHT < G + mb + L - c, PAM-left ones (into the window before) iff LEFT < G + mb + L - c.

Each check prints one line per search set and exits 1 on a difference.
"""
import bisect
import os
import random
import sys

GR = "GATCCTAGGCATTACGCATG"  # the planted guide region (20 nt)
PAM = "AAGAAT"  # an NNGRRT
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
FORMS = ("win", "trail", "lead", "both")
ALT = 3


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def write_genome(out, tag, left, right, kinds, plant, rng):
    """windows of left + ALT + right bases, one per kind ('-' dense, else a site that plant
    writes into the window), one N between two"""
    width = left + ALT + right
    seq, rows, pos = [], [], 0
    for i, kind in enumerate(kinds):
        w = [rng.choice("ACGT") for _ in range(width)]
        if kind != "-":
            plant(w, kind)
        if i:
            seq.append("N")
            pos += 1
        rows.append(("fakechrA", pos, pos + width, left, right, 1 if kind == "-" else 0, kind))
        seq.append("".join(w))
        pos += width
    contig = "".join(seq)
    d = os.path.join(out, tag)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "fakechrA.fa"), "w") as fh:
        fh.write(">fakechrA\n")
        for i in range(0, len(contig), 60):
            fh.write(contig[i:i + 60] + "\n")
    with open(os.path.join(out, tag + ".windows.tsv"), "w") as fh:
        for r in rows:
            fh.write("%s\t%d\t%d\t%d\t%d\t%d\t%s\n" % r)
    for form in FORMS:
        with open(os.path.join(out, "%s.%s.bed" % (tag, form)), "w") as fh:
            for k, (c, s, e, _, _, dense, _) in enumerate(rows):
                if not dense:
                    continue
                lead = 1 if form in ("lead", "both") and k > 0 else 0  # the separator before it
                trail = 1 if form in ("trail", "both") and k < len(rows) - 1 else 0  # the one after
                fh.write("%s\t%d\t%d\n" % (c, s - lead, e + trail))


def genome(out):
    rng = random.Random(1327)

    def plant_route(w, kind):  # case 1: the guide region o bases from the edge its core reaches
        o = int(kind[1:])
        if kind[0] == "R":
            w[0:o + 26] = list("A" * o + GR + PAM)
        else:
            s = revcomp("A" * o + GR + PAM)
            w[len(w) - len(s):] = list(s)

    def plant_pam(w, kind):  # case 2: the guide region j bases from the edge on its PAM side
        j = int(kind[1:])
        if kind[0] == "R":
            w[len(w) - j - 23:len(w) - j] = list(GR[:3] + GR)
        else:
            w[j:j + 23] = list(revcomp(GR[:3] + GR))

    def kinds(n):
        k = ["-"]
        for side in "RL":
            for o in range(n):
                k += ["%s%d" % (side, o), "-"]  # every site between two dense windows
        return k

    for left, right, tag in ((25, 27, "hr25"), (30, 30, "hr30")):
        write_genome(out, tag, left, right, kinds(6), plant_route, rng)
    for w in (28, 29, 30):
        write_genome(out, "st%d" % w, w, w, kinds(8), plant_pam, rng)
    for gs in range(21, 27):
        with open(os.path.join(out, "guides_hr%d.txt" % gs), "w") as fh:
            fh.write("A" * (gs - 20) + GR + "N" * len(PAM) + "\n")
    with open(os.path.join(out, "guides_st.txt"), "w") as fh:
        fh.write(GR + "N" * 8 + "\n")
    print("hr25, hr30, st28, st29, st30 written to %s" % out)


def read_lines(path):
    with open(path) as fh:
        return [l for l in fh if not l.startswith("#")]


def windows(out, tag):
    rows = [l.split() for l in open(os.path.join(out, tag + ".windows.tsv"))]
    return rows, [int(r[1]) for r in rows], int(rows[0][3]), int(rows[0][4])


def line_span(line):
    f = line.rstrip("\n").split("\t")
    return int(f[4]), len(f[2]) - f[2].count("-")


def check(out, tag, work, G, L, mb, brna, gss):
    rows, starts, left, right = windows(out, tag)
    sites = [r[6] for r in rows if r[6] != "-"]

    def site_of(line):
        """the site of the window that holds the whole line, and whether it covers its ALT"""
        p, span = line_span(line)
        k = bisect.bisect_right(starts, p) - 1
        if k < 0:
            return None, False
        s, e = int(rows[k][1]), int(rows[k][2])
        if p + span > e:
            return None, False
        return rows[k][6], p < e - right and p + span > s + left

    bad = []
    for gs in gss:
        U = read_lines(os.path.join(work, "%s.%d.U.targets.txt" % (tag, gs)))
        has_lines = {s for s in sites if gs - G - min(int(s[1:]), mb) <= brna}
        seen, across = set(), 0
        for line in U:
            site, _ = site_of(line)
            if site is None:
                across += 1
            elif site != "-":
                seen.add(site)
        if across:
            bad.append("%s Gs %d: %d lines run across a separator" % (tag, gs, across))
        if seen != has_lines:
            bad.append("%s Gs %d: the sites with lines are %s, the rule's %s"
                       % (tag, gs, sorted(seen), sorted(has_lines)))
        out_line = []
        for form in FORMS:
            S = set(read_lines(os.path.join(work, "%s.%d.%s.targets.txt" % (tag, gs, form))))
            missed = set()
            for line in U:
                if line in S:
                    continue
                site, covers = site_of(line)
                if site not in (None, "-") and covers:
                    missed.add(site)
            want = set()
            for s in has_lines:
                o = int(s[1:])
                if s[0] == "R":  # its core reaches the window before
                    c = 1 if form in ("trail", "both") else 2
                    if o <= gs - G - c and o + G + L > left:
                        want.add(s)
                else:  # the window after
                    c = 1 if form in ("lead", "both") else 2
                    if o <= gs - G - c and o + G + L > right:
                        want.add(s)
            if missed != want:
                bad.append("%s Gs %d %s: missed %s, the rule says %s" % (tag, gs, form, sorted(missed), sorted(want)))
            out_line.append("%s %s" % (form, ",".join(sorted(missed)) or "none"))
        print("%s n=%d (Gs %d): %d lines; sites missed: %s" % (tag, gs + L, gs, len(U), "; ".join(out_line)))
    for b in bad:
        print("  FAIL: " + b)
    return 1 if bad else 0


def check_pam(out, tag, work, G, L, mb):
    rows, starts, left, right = windows(out, tag)
    m = G + mb + L
    U = read_lines(os.path.join(work, "%s.U.targets.txt" % tag))
    bad, across, out_line = [], 0, []
    for line in U:
        p, span = line_span(line)
        if bisect.bisect_right(starts, p + span - 1) != bisect.bisect_right(starts, p):
            across += 1
    if not across:
        bad.append("%s: no line runs across a separator" % tag)
    for form in FORMS:
        S = set(read_lines(os.path.join(work, "%s.%s.targets.txt" % (tag, form))))
        missed = set()
        for line in U:
            if line in S:
                continue
            p, span = line_span(line)
            for r in rows[max(0, bisect.bisect_right(starts, p) - 1):bisect.bisect_right(starts, p + span - 1)]:
                s, e = int(r[1]), int(r[2])
                if r[5] == "0" and p < e - right and p + span > s + left:  # the ALT of a window not dense
                    missed.add(r[6][0])
        want = set()
        if right < m - (1 if form in ("lead", "both") else 2):  # PAM-right lines: the window after
            want.add("R")
        if left < m - (1 if form in ("trail", "both") else 2):  # PAM-left lines: the window before
            want.add("L")
        if missed != want:
            bad.append("%s %s: lines missed on side(s) %s, the rule says %s"
                       % (tag, form, "".join(sorted(missed)) or "none", "".join(sorted(want)) or "none"))
        out_line.append("%s %s" % (form, "".join(sorted(missed)) or "none"))
    print("%s (%d/%d, G + mb + L = %d): %d lines, %d across a separator; sides missed: %s"
          % (tag, left, right, m, len(U), across, "; ".join(out_line)))
    for b in bad:
        print("  FAIL: " + b)
    return 1 if bad else 0


def main():
    a = sys.argv
    if len(a) == 3 and a[1] == "genome":
        genome(a[2])
        return 0
    if len(a) >= 10 and a[1] == "check":
        G, L, mb, brna = (int(x) for x in a[5:9])
        return check(a[2], a[3], a[4], G, L, mb, brna, [int(x) for x in a[9:]])
    if len(a) == 8 and a[1] == "check-pam":
        return check_pam(a[2], a[3], a[4], int(a[5]), int(a[6]), int(a[7]))
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())
