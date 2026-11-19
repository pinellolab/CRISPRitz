#!/usr/bin/env python3
"""Checks of the region hooks (docs/INDEX_FORMAT.md section 13) on one search.

usage: check_hooks.py --unskipped U --skip S --skip-out S.stdout --keep K --keep-out K.stdout
         --bed BED --pam-side 3|5 --G G --K K --mb MB --gs GS --L L [--prune]
         [--skip-noprune S0] [--dump D --genome DIR] [--leaves N] [--ordered]
         [--min-old-lost N] [--min-old-dup N] [--windows W.tsv --hr27 any|none]
         [--min-dropped N] [--min-kept N] [--label TEXT]

U, S, K, S0 are the targets files of the search without a BED, with the skip BED, with
the same BED as a keep BED, and with the skip BED and CRISPRITZ_SKIP_PRUNE=0. Everything
is computed from the lines and the FASTA, sharing nothing with the engine:

* THE core of a line: n = len(crRNA) - d (d gaps in crRNA, r in DNA); a PAM-right line
  (a 3' PAM on +, a 5' PAM on -) has the core [p + d - r, p + d - r + n), a PAM-left one
  [p, p + n). S must be exactly the lines of U whose core meets no interval, and K those
  whose core meets one (in U's order with --ordered: single-threaded searches), S0 = S;
* the leaf of a line (its window start ws): PAM-right ws = p + d - r - mb - off, PAM-left
  ws = p - pp (off = G - Gs, pp = K - L under --prune with a 3' PAM, else 0);
* the SKIPPED/KEPT LEAVES lines: S's and K's add up to --leaves, and are each other's
  mirror; the leaf dump holds exactly SKIPPED leaves, the LEAF DUMP count, each one's
  core meeting an interval, every leaf of a dropped line of U among them, and every
  letter that of the FASTA (upper case; a byte that is no base, IUPAC code or Q read as
  N) but '-' exactly where docs/INDEX_FORMAT.md section 13.3 says: a pad base at or past
  pad_valid (the pad bases from the guide outward before the first one off the contig or
  N), a flank base off the contig and the N next to it at the flank's distal end;
* --min-old-lost / --min-old-dup: at least that many lines of U on which the earlier
  prototype's leaf rule ([ws, + G + K + pp) PAM-left, [ws + mb, + G + K + pp) PAM-right)
  and THE core disagree: lines it would have dropped though their core meets no interval
  (lost), and lines it would have kept though it does (duplicated). The engine must follow
  THE core on all of them (the first check);
* --windows/--hr27: on a fake-indel genome, the dropped lines whose alignment lies wholly
  in a window that is not dense and covers its ALT (docs/INDEX_FORMAT.md section 13.5):
  'any' requires at least one, 'none' none.
Prints one line of counts; exits 1 naming each failure.
"""
import argparse
import bisect
import os
import sys

IUPAC = set("ACGTRYSWKMBDHVNQ")


def read_lines(path):
    with open(path) as fh:
        return [l for l in fh if not l.startswith("#")]


def read_bed(path):
    iv = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c, s, e = line.rstrip("\r\n").split("\t")[:3]
            iv.setdefault(c, []).append((int(s), int(e)))
    merged = {}
    for c, v in iv.items():
        v.sort()
        m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[c] = ([a for a, b in m], [b for a, b in m])
    return merged


def meets(bed, c, s, e):
    if c not in bed:
        return False
    starts, ends = bed[c]
    k = bisect.bisect_right(starts, s) - 1
    if k >= 0 and ends[k] > s:
        return True
    return k + 1 < len(starts) and starts[k + 1] < e


def read_genome(d):
    g = {}
    for name in sorted(os.listdir(d)):
        if not name.endswith(".fa"):
            continue
        with open(os.path.join(d, name)) as fh:
            head = fh.readline().strip()[1:].split()[0]
            seq = "".join(l.strip() for l in fh).upper()
        g[head] = "".join(ch if ch in IUPAC else "N" for ch in seq)
    return g


def counts(stdout_path):
    s = k = dump = None
    with open(stdout_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "SKIPPED LEAVES" and len(f) == 4 and f[2] == "KEPT LEAVES":
                if s is not None:
                    return None, None, "two SKIPPED LEAVES lines"
                s, k = int(f[1]), int(f[3])
            elif f[0] == "LEAF DUMP" and len(f) == 2:
                dump = int(f[1])
    return (s, k), dump, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unskipped", required=True)
    ap.add_argument("--skip", required=True)
    ap.add_argument("--skip-out", required=True)
    ap.add_argument("--keep", required=True)
    ap.add_argument("--keep-out", required=True)
    ap.add_argument("--skip-noprune")
    ap.add_argument("--bed", required=True)
    ap.add_argument("--pam-side", required=True, choices=("3", "5"))
    ap.add_argument("--G", type=int, required=True)
    ap.add_argument("--K", type=int, required=True)
    ap.add_argument("--mb", type=int, required=True)
    ap.add_argument("--gs", type=int, required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--prune", action="store_true")
    ap.add_argument("--dump")
    ap.add_argument("--genome")
    ap.add_argument("--leaves", type=int)
    ap.add_argument("--ordered", action="store_true")
    ap.add_argument("--min-old-lost", type=int, default=0)
    ap.add_argument("--min-old-dup", type=int, default=0)
    ap.add_argument("--min-dropped", type=int, default=1)
    ap.add_argument("--min-kept", type=int, default=1)
    ap.add_argument("--windows")
    ap.add_argument("--hr27", choices=("any", "none"))
    ap.add_argument("--dup-ok", action="store_true",
                    help="a leaf may be dumped twice, in identical lines (a 2.9.1 PAM-less index holds index 0 twice)")
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    side3 = a.pam_side == "3"
    n_want = a.gs + a.L
    off = a.G - a.gs
    pp = a.K - a.L if (a.prune and side3) else 0
    bed = read_bed(a.bed)
    bad = []

    U, S, Kp = read_lines(a.unskipped), read_lines(a.skip), read_lines(a.keep)

    def parse(line):
        f = line.rstrip("\n").split("\t")
        cr, dna, c, p, strand = f[1], f[2], f[3], int(f[4]), f[6]
        d, r = cr.count("-"), dna.count("-")
        n = len(cr) - d
        right = (strand == "+") == side3
        cs = p + d - r if right else p
        ws = p + d - r - a.mb - off if right else p - pp
        return c, cs, n, ws, "R" if right else "L", p, len(dna) - r

    inside, outside = [], []
    leaves_of_dropped = set()
    old_lost = old_dup = 0
    bad_n = 0
    old_len = a.G + a.K + pp
    for line in U:
        c, cs, n, ws, sd, p, span = parse(line)
        if n != n_want:
            bad_n += 1
        m = meets(bed, c, cs, cs + n)
        (inside if m else outside).append(line)
        if m:
            leaves_of_dropped.add((c, ws, sd))
        o0 = ws if sd == "L" else ws + a.mb
        om = meets(bed, c, o0, o0 + old_len)
        if om and not m:
            old_lost += 1
        elif m and not om:
            old_dup += 1
    if bad_n:
        bad.append("%d lines of U have a core of other than Gs + L = %d bases" % (bad_n, n_want))
    if len(set(U)) != len(U):
        bad.append("U holds %d duplicate lines" % (len(U) - len(set(U))))
    if a.ordered:
        if S != outside:
            bad.append("S is not U's lines whose core meets no interval, in U's order (%d vs %d lines; "
                       "%d of S not among them, %d of them not in S)"
                       % (len(S), len(outside), len(set(S) - set(outside)), len(set(outside) - set(S))))
        if Kp != inside:
            bad.append("K is not U's lines whose core meets an interval, in U's order (%d vs %d lines)"
                       % (len(Kp), len(inside)))
    else:
        if sorted(S) != sorted(outside):
            bad.append("S is not the lines of U whose core meets no interval (%d vs %d; %d lost, %d extra)"
                       % (len(S), len(outside), len(set(outside) - set(S)), len(set(S) - set(outside))))
        if sorted(Kp) != sorted(inside):
            bad.append("K is not the lines of U whose core meets an interval (%d vs %d)" % (len(Kp), len(inside)))
    if len(inside) < a.min_dropped:
        bad.append("only %d lines dropped by the BED (want >= %d): the case tests nothing" % (len(inside), a.min_dropped))
    if len(outside) < a.min_kept:
        bad.append("only %d lines kept (want >= %d)" % (len(outside), a.min_kept))
    if old_lost < a.min_old_lost:
        bad.append("%d lines on which the prototype's rule loses a line (want >= %d)" % (old_lost, a.min_old_lost))
    if old_dup < a.min_old_dup:
        bad.append("%d lines on which the prototype's rule duplicates a line (want >= %d)" % (old_dup, a.min_old_dup))
    if a.skip_noprune is not None:
        with open(a.skip_noprune) as f0, open(a.skip) as f1:
            if f0.read() != f1.read():
                bad.append("CRISPRITZ_SKIP_PRUNE=0 gives another output than the pruned skip")

    sk, dump_count, err = counts(a.skip_out)
    if err:
        bad.append("the skip search printed " + err)
    kk, _, err = counts(a.keep_out)
    if err:
        bad.append("the keep search printed " + err)
    if sk is None or sk[0] is None:
        bad.append("the skip search printed no SKIPPED LEAVES line")
    if kk is None or kk[0] is None:
        bad.append("the keep search printed no SKIPPED LEAVES line")
    if sk and kk and sk[0] is not None and kk[0] is not None:
        if sk[0] != kk[1] or sk[1] != kk[0]:
            bad.append("skip's SKIPPED/KEPT %d/%d are not keep's KEPT/SKIPPED %d/%d" % (sk[0], sk[1], kk[1], kk[0]))
        if a.leaves is not None and sk[0] + sk[1] != a.leaves:
            bad.append("SKIPPED + KEPT = %d, the index holds %d leaves" % (sk[0] + sk[1], a.leaves))

    n_dump = n_letters = n_dash = 0
    if a.dump:
        G = read_genome(a.genome)
        width = a.K + a.G + a.mb
        seen = {}
        conflicts = wrong = bad_core = bad_len = 0
        dumped = set()
        dlines_seen = set()
        examples = []
        with open(a.dump) as fh:
            dlines = fh.read().splitlines()
        if dlines != sorted(dlines, key=lambda x: (x.split("\t")[0], int(x.split("\t")[1]), x.split("\t")[3], x.split("\t")[2])):
            bad.append("the leaf dump is not sorted by contig, window start, side")
        for line in dlines:
            c, ws, w, sd = line.split("\t")
            ws = int(ws)
            n_dump += 1
            if (c, ws, sd) in dumped and not (a.dup_ok and line in dlines_seen):
                bad.append("the leaf dump lists %s %d %s twice" % (c, ws, sd))
            dumped.add((c, ws, sd))
            dlines_seen.add(line)
            if len(w) != width:
                bad_len += 1
                continue
            cs = ws + pp if sd == "L" else ws + a.mb + off
            if not meets(bed, c, cs, cs + n_want):
                bad_core += 1
            seq = G[c]
            Lc = len(seq)
            if sd == "L":
                flank = list(range(ws, ws + a.K))  # distal first
                guide0 = ws + a.K
                pad = [guide0 + a.G + j for j in range(a.mb)]  # from the guide outward
            else:
                guide0 = ws + a.mb
                flank = list(range(guide0 + a.G + a.K - 1, guide0 + a.G - 1, -1))  # distal first
                pad = [guide0 - 1 - j for j in range(a.mb)]
            dash = set()
            valid = 0
            for q in pad:
                if 0 <= q < Lc and seq[q] != "N":
                    valid += 1
                else:
                    break
            dash.update(pad[valid:])
            if any(not 0 <= q < Lc for q in flank):  # FLANK_OFF: the distal run off the contig or N
                for q in flank:
                    if not 0 <= q < Lc or seq[q] == "N":
                        dash.add(q)
                    else:
                        break
            for i, ch in enumerate(w):
                q = ws + i
                want = "-" if (q in dash or not 0 <= q < Lc) else seq[q]
                if ch != want:
                    wrong += 1
                    if len(examples) < 3:
                        examples.append("%s:%d dumped %r, want %r (%s)" % (c, q, ch, want, line))
                if ch == "-":
                    n_dash += 1
                    continue
                n_letters += 1
                old = seen.setdefault((c, q), ch)
                if old != ch:
                    conflicts += 1
        if bad_len:
            bad.append("%d dumped windows are not K + G + mb = %d letters" % (bad_len, width))
        if bad_core:
            bad.append("%d dumped leaves have a core that meets no interval" % bad_core)
        if wrong:
            bad.append("%d dumped letters differ from the FASTA and its '-' rule: %s" % (wrong, "; ".join(examples)))
        if conflicts:
            bad.append("%d letters differ between overlapping dumped windows" % conflicts)
        missing = leaves_of_dropped - dumped
        if missing:
            bad.append("%d leaves of lines the skip dropped are not in the dump, e.g. %s"
                       % (len(missing), sorted(missing)[:3]))
        if sk and sk[0] is not None and n_dump != sk[0]:
            bad.append("the dump holds %d leaves, SKIPPED LEAVES says %d" % (n_dump, sk[0]))
        if dump_count != n_dump:
            bad.append("LEAF DUMP says %s, the dump holds %d" % (dump_count, n_dump))

    hr27 = 0
    if a.windows:
        wins = {}
        for line in open(a.windows):
            c, s, e, left, right, dense = line.split()
            wins.setdefault(c, []).append((int(s), int(e), int(left), int(right), dense == "1"))
        wstart = {c: [w[0] for w in v] for c, v in wins.items()}
        for line in set(U) - set(S):
            c, cs, n, ws, sd, p, span = parse(line)
            k = bisect.bisect_right(wstart[c], p) - 1
            if k < 0:
                continue
            s, e, left, right, dense = wins[c][k]
            if p + span <= e and not dense and p < e - right and p + span > s + left:
                hr27 += 1
        if a.hr27 == "any" and hr27 == 0:
            bad.append("no dropped line lies in a window that is not dense and covers its ALT (want some)")
        if a.hr27 == "none" and hr27:
            bad.append("%d dropped lines lie in a window that is not dense and cover its ALT" % hr27)

    print("%s: U %d = S %d + K %d; old rule: %d lost, %d duplicated; SKIPPED/KEPT %s/%s; dump %d leaves, "
          "%d letters, %d '-'; H-R2-7 lines %d; %s"
          % (a.label, len(U), len(S), len(Kp), old_lost, old_dup, sk[0] if sk else None,
             sk[1] if sk else None, n_dump, n_letters, n_dash, hr27, "OK" if not bad else "FAIL"))
    for b in bad:
        print("  FAIL: " + b)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
