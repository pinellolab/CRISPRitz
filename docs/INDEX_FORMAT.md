# CRISPRitz index format and target windows (format 2, CRISPRitz 2.9.2)

This document specifies what a CRISPRitz TST index holds, how it is laid out on
disk, and how the search reads it. It covers the boundary/window fixes planned for
CRISPRitz 2.9.2 (audit ids IDX-13, IDX-14, IDX-15 and IDX-4) and the index-format
version stamp (B5) that those fixes need, because two of them change what an index
holds.

Status: design, implemented on branch `feat/engine-boundaries`. Section 10 holds the
measurements; section 11 the tests.

| id | audit | what goes wrong in 2.9.1 | fix | needs a rebuild |
|----|-------|--------------------------|-----|-----------------|
| B1 | IDX-14 | a leaf with index 0 is built as a minus-strand window (3' PAM) but read as plus strand, so a reverse-strand PAM at a contig's first base is reported on the + strand at base `max_bulges` | search side: read index 0 as PAM-left for a 3' PAM; drop the index-0 leaf the builder made from a plus-strand window | no (works on format-1 indexes) |
| B2 | IDX-15 | a PAM that ends on a contig's last base is never found, by any index or by the brute-force scan | scan bound `nt <= len - pam_limit` on both PAM sides, in `buildTST` and `searchBruteForce` | for indexes yes; brute force no |
| B3 | IDX-13 | a site whose guide lies within `max_bulges` nt of a contig start/end or an N run is dropped, even a perfect match: the whole pad+guide+PAM window must fit and be free of N | emit a site when its guide+PAM core fits; store the pad positions that are off the contig or on an N as non-consumable | yes (format 2) |
| B4 | IDX-4 | `--prune-pam` with flank `K` wider than the PAM `L` loses + strand PAMs with 1..K-L bases to the contig end and - strand PAMs 0..K-L-1 from the start, and reports a - strand PAM exactly K-L from the start on the + strand | store a PAM-less window whose distal flank runs off the contig (flank bases there stored as N); store the window start as a signed coordinate with the strand in its own bit | yes (format 2) |
| B5 | - | an index does not say what it holds | a header with a magic number, a format version and the build parameters; old indexes keep loading and searching as with 2.9.1 | - |

## 1. Terms

* `G` guide length (20 for SpCas9), `K` stored PAM/flank width (the build PAM
  file's declared length), `L` the search PAM's declared length (`L = K` except
  under `--prune-pam`), `mb` the index's `max_bulges` (`-bMax`).
* A **site** is one PAM occurrence on one strand. Its **window** is the stretch of
  the contig the index stores for it: the guide region (`G` bases), the PAM/flank
  (`K` bases) on its PAM side, and a **pad** of `mb` bases on the far (distal) side
  of the guide, which only DNA-bulge alignments consume.
* Seen on the + strand, a window is one of two shapes:

  | window | layout on the + strand | used for |
  |--------|------------------------|----------|
  | PAM-left  | `[flank K][guide G][pad mb]` | 3' PAM on the - strand; 5' PAM on the + strand |
  | PAM-right | `[pad mb][guide G][flank K]` | 3' PAM on the + strand; 5' PAM on the - strand |

  `ws`, the **window start**, is the + strand coordinate of the window's leftmost
  base (0-based).
* The **core** of a site is its guide region plus the PAM bases the search
  constrains. The pad is not part of it.
* The TST stores, per site, a string of `G + mb` characters read from the PAM
  outward (guide first, then pad), and per leaf the `K` flank bases as 4-bit
  nibbles.
* Reported coordinates (unchanged by this work): for a PAM-left window the PAM's
  first base, `ws + (K - L)` for a 3' PAM, `ws` for a 5' one; for a PAM-right
  window the leftmost aligned base, `ws + mb - dna + rna` (`dna`/`rna` = bulges
  used). Cluster Position: the same as Position for PAM-left windows, the guide
  start `ws + mb` for PAM-right windows.

## 2. Format 1 (CRISPRitz up to 2.9.1)

```
int32  numLeaves                 (>= 1)
int32  len_guide                 (G)
numLeaves x {
    int32  guideIndex            > 0: PAM-right, ws = guideIndex
                                 < 0: PAM-left,  ws = -guideIndex
                                 = 0: ambiguous (see below)
    ceil(K/2) bytes              flank, two nibbles a byte (A=1 C=2 G=4 T=8, IUPAC
                                 sums, N=0); K comes from the file name
    uint8  '0'                   no next leaf with the same string
         | '_' int32 next        -(i+1) of the next leaf with the same string
}
int32  numNodes
tree   pre-order, 4-bit symbols: base (1..14), 0 = no lo/hi kid, 15 = end of
       string, followed by int32 -(i+1) of the first leaf of its chain
```

The file name is `<PAM>_<contig>_<part>.bin`; `<PAM>` (the build PAM region) is the
only record of `K` and of the PAM.

Defects of format 1 and its builder:

* **Index 0 is ambiguous (B1, IDX-14).** The builder pushes `-nt` for a PAM-left
  window and `ws` for a PAM-right one, and classifies index 0 as PAM-left for a 3'
  PAM (`mainParallel.cpp`, `if (pamIndices[i] > 0)`) and as PAM-right for a 5' PAM
  (`if (pamIndices[i] < 0)`). The search reads only `< 0` as PAM-left
  (`searchOnTST.cpp:718`, `:907`). So for a 3' PAM the leaf at index 0 is written
  as a PAM-left window and read as a PAM-right one: a - strand PAM at base 0 is
  reported as `+` at base `mb` (or `mb - dna + rna`), and with `--prune-pam` a -
  strand PAM exactly `K-L` from the start lands there too. A window of the other
  shape that starts at 0 -- for a 3' PAM a PAM-right one (a + strand site whose guide
  starts at exactly `mb`), for a 5' PAM a PAM-left one (a + strand PAM at base 0) --
  is built as the shape index 0 is read as, from the wrong bases: its real site is
  lost, and the stored artefact can surface as a match whose PAM is not there.
* **The scan stops one short (B2, IDX-15).** `pamGenerator.cpp:344` and
  `pam_searching.cpp:570` scan `nt < len - pam_limit`, so a PAM on the last
  `pam_limit` bases is never found. The 5' branches bound on `- pam_len`
  (`pamGenerator.cpp:387`, `pam_searching.cpp:612`), which also loses every
  reverse-strand 5' PAM within `G` bases of the end.
* **The whole window must fit and be free of N (B3, IDX-13).** A PAM-right window
  is kept only if `ws >= 0` and a PAM-left one only if `ws + G + K + mb <= len`
  (`pamGenerator.cpp:371`, `:378`, `:414`, `:421`), and `guideRegionHasN`
  (`mainParallel.cpp:68`) discards a window with an N anywhere outside its PAM, pad
  included. A 0-bulge alignment needs none of the pad, so sites within `mb` nt of a
  contig end or an N run are lost, and an index with a larger `mb` loses more.
* **The flank must fit (B4, IDX-4).** For a PAM-less index the flank is `K` wide but
  a `--prune-pam` search constrains only its `L` proximal bases; the window still
  has to fit, so PAMs within `K-L` nt of the end (+ strand) or the start (- strand)
  are lost. A PAM-left window starting before base 0 cannot even be written: the
  sign of `guideIndex` is the strand.

## 3. What a format-2 index holds

Let `n` be the contig length. Positions outside `[0, n)` are **off the contig**.

### 3.1 Which sites are stored

* **Dedicated PAM** (any build PAM with at least one non-N letter): every position
  `nt` in `[0, n-K]` (B2: inclusive) where the `K` bases match the PAM (forward, or
  reverse complement), as in 2.9.1 (a genomic N matches no PAM letter, not even
  `N`), **and** whose guide region lies inside the contig. The pad may run off the
  contig or into an N.
* **PAM-less index** (build PAM all `N`, the index `--prune-pam` serves): every
  window whose guide region lies inside the contig and whose most proximal flank
  base (the one next to the guide) is on the contig. The rest of the flank may run
  off the contig; those bases are stored as `N` (nibble 0) and the leaf is flagged
  `FLANK_OFF`. A window whose flank is entirely off the contig is not stored: no
  PAM of length >= 1 can match it.
* In both cases a window whose **guide region** holds an N is discarded, as in
  2.9.1: the tree has no symbol for N and the core is unusable.

### 3.2 The pad

`pad_valid` (0..`mb`) is the number of pad positions, counted from the guide
outward, before the first one that is off the contig or an N. An alignment that
consumes `dna - rna > pad_valid` pad bases needs a base that is not there and is not
reported. Pad positions from `pad_valid` on are stored in the tree as `A`: no
reported alignment can depend on them (an alignment reaching them is rejected by the
`pad_valid` test at leaf emission), so any symbol serves, and the 14 nibble values
of the tree's alphabet stay as they are. For every site of 2.9.1 (window inside the
contig, no N outside the PAM), `pad_valid = mb` and the stored string is unchanged.

### 3.3 The flank

As in 2.9.1 the flank holds the genome's bases, with a genomic N as nibble 0. A
flank base off the contig is also nibble 0. Under `--prune-pam` only the `L`
positions the search PAM covers are tested, and nibble 0 matches no PAM letter, so
a site whose real PAM runs off the contig is rejected exactly as a dedicated index
never stores it. Without `--prune-pam`, a PAM-less index's flank is not tested at
all (2.9.1 semantics: a genomic N in the flank does not drop the site), but a leaf
flagged `FLANK_OFF` is skipped, so a non-pruning search returns the windows 2.9.1
could hold, plus the B2/B3 sites, and never a coordinate before base 0.

## 4. Format 2 on disk

All integers are 32-bit, in the byte order of the host that wrote the file (little
endian on every platform CRISPRitz is built for), exactly as format 1.

```
header
  off  size  field
  0    4     magic        bytes 43 5A 49 FE ("CZI" 0xFE); as an int32 it is
                          negative (-28,747,197), so a format-1 reader refuses it
  4    4     version      2
  8    4     header_len   28 + K: offset of numLeaves
  12   4     flags        bit 0: PAM at the 5' end; bit 1: PAM-less (all-N) build
                          PAM; other bits 0 (reserved, must be 0)
  16   4     guide_len    G
  20   4     pam_len      K
  24   4     max_bulges   mb (0..31)
  28   K     pam          the build PAM region, ASCII upper case ("NGG",
                          "NNNNNNNN", "TTTV")
body (from header_len)
  int32  numLeaves        (>= 1)
  int32  len_guide        G again (format 1 keeps it here)
  numLeaves x {
      int32  ws           window start, signed (may be < 0 for a window whose pad
                          or flank runs off the contig start)
      ceil(K/2) bytes     flank nibbles, as format 1
      uint8  lflags       bit 7 HAS_NEXT, bit 6 PAM_LEFT, bit 5 FLANK_OFF,
                          bits 0-4 pad_valid (0..31)
      int32  next         present iff HAS_NEXT: -(i+1) of the next leaf with the
                          same string
  }
  int32  numNodes
  tree                    unchanged
trailer
  int32  magic            the same 4 bytes; a file not ending in them is cut short
```

Notes:

* The per-leaf record is the same size as in format 1 (the `lflags` byte takes the
  place of the `'0'`/`'_'` marker), so a format-2 index differs in size from a
  format-1 index of the same genome only by the sites B2/B3/B4 add and 32 + K bytes
  of header and trailer per `.bin` (section 10 has the measurement).
* The tree serialization, the sort and the leaf chains are unchanged. A contig with
  no boundary site therefore has the same tree as in 2.9.1 and every search over it
  returns the same rows in the same order.
* `mb` is capped at 31 by the 5-bit `pad_valid`; `buildTST` refuses a larger
  `-bMax` with an error. (The search's guide buffers are 100 characters, so such
  values never worked in practice.)
* `buildTST` checks every write and the final close, and exits 1 naming the file
  if one fails (a full disk used to leave a cut-short `.bin` behind an exit 0).

## 5. How the search reads an index

`loadTST` reads the first 4 bytes of each `.bin`. If they are the magic it parses
the header and reads the body as format 2; otherwise they are `numLeaves` and the
file is format 1. Each `.bin` is decoded on its own, so a directory mixing formats
still searches correctly (though `buildTST` never makes one).

The search refuses, with an error on stderr and exit 1:

* a format version it does not know (> 2), or reserved header bits set;
* a header `pam_len` that differs from the file name's PAM prefix, or a PAM side
  (3'/5') that differs from the search PAM's;
* a format-2 file that does not end in the trailer (cut short).

If `max_bulges` on the command line differs from a format-2 header's, the header
wins for decoding coordinates and a note goes to stdout (stderr is fatal in
CRISPRme). Before searching, `searchTST` prints the formats it found on stdout, e.g.
`INDEX FORMAT: 2 x24` or `INDEX FORMAT: 1 (CRISPRitz <= 2.9.1) x501`.

Per leaf the search derives `PAM_LEFT`, `ws`, `pad_valid` and `FLANK_OFF`:

| | format 1 | format 2 |
|-|----------|----------|
| `PAM_LEFT` | `guideIndex < 0`, or `guideIndex == 0` and the PAM is 3' (B1) | `lflags` bit 6 |
| `ws` | `abs(guideIndex)` | `ws` |
| `pad_valid` | `mb` (2.9.1 kept only whole windows) | `lflags` bits 0-4 |
| `FLANK_OFF` | 0 | `lflags` bit 5 |
| skip | an index-0 leaf whose flank does not match the build PAM (B1) | - |

At emission a leaf is skipped when `--prune-pam` is on and its flank fails the
search PAM (unchanged), when the index is PAM-less, `--prune-pam` is off and the
leaf is `FLANK_OFF`, or when the alignment used `dna - rna > pad_valid` pad bases.
Coordinates are then computed as in section 1, with `mb` from the header for format
2 and from the command line for format 1 (as 2.9.1 does).

The B1 skip: for a 3' PAM, index 0 holds either a real PAM-left window at 0 (its
flank matches the build PAM, which the scan required) or a PAM-right window at 0
that the 2.9.1 builder built as PAM-left from the wrong bases. Such a leaf's flank
is the complement of the contig's first `K` bases and matches the build PAM only if
those bases do hold a reverse PAM, in which case it is the identical real window
(the scan found it too) and is kept. The mirror holds for a 5' PAM. For a PAM-less
build PAM every window is real and nothing is skipped. So on a format-1 index B1
changes only the rows of index-0 leaves: real ones get their strand and coordinate
right, artefacts are dropped. A format-1 index never holds the lost sites; only a
rebuild adds them.

## 6. Compatibility

| | format-1 index (built by <= 2.9.1) | format-2 index (built by 2.9.2) |
|-|------------------------------------|----------------------------------|
| CRISPRitz <= 2.9.0 | as before | searchTST aborts: `std::length_error` (`vector::_M_default_append`) on stderr, exit 134, empty targets file (tested with the 2.9.0 binary of the CRISPRme 2.8.0 image; its `crispritz.py` does not check exit statuses, CRISPRme fails the stage on the stderr line) |
| CRISPRitz 2.9.1 | as before | refused: `'...' is not a valid CRISPRitz index (.bin) file (numLeaves=-28747197)`, exit 1 (tested) |
| CRISPRitz 2.9.2 searchTST | byte-identical to 2.9.1 except index-0 rows (B1) | full B1-B4 semantics |
| CRISPRitz 2.9.2 searchBruteForce | - (no index) | - ; B2 applies at once |

Byte-identical means the same `.targets.txt` bytes, row order included; tested on
the repo's fixtures and on indexes built by the 2.9.1 binary (section 11).
Users' existing indexes keep working unchanged. They keep missing the sites of
B2-B4 until rebuilt; B1 corrects their index-0 rows at once.

A format-2 index of a genome holds every site the format-1 index of the same
genome and parameters holds, with the same string, flank and coordinates, plus the
boundary sites (measured on hg38 chr22 and chr1, section 10). A 2.9.2 search on it
returns every row a 2.9.2 search on the format-1 index returns, plus rows at sites
whose guide lies within `mb` nt of a contig end or an N run (B3, including a + strand
3' site whose guide starts at exactly `mb`), whose PAM lies within `K - L` nt of a
contig end under `--prune-pam` (B4), or whose PAM ends on a contig's last base (B2).

## 7. Reading the stamp without loading the tree

Any reader (CRISPRme, scripts) can classify an index from the first bytes of one
`.bin`. Reference implementation (also `crispritz.py index-info <index_dir>`, which
prints this for every `.bin` as JSON):

```python
import struct

CZI_MAGIC = b"CZI\xfe"

def read_index_header(path):
    """Format and build parameters of a CRISPRitz .bin, without reading the tree."""
    with open(path, "rb") as fh:
        head = fh.read(28)
        if head[:4] != CZI_MAGIC:              # format 1: numLeaves, len_guide
            leaves, guide = struct.unpack("=ii", head[:8])
            return {"format": 1, "num_leaves": leaves, "guide_len": guide}
        version, header_len, flags, guide, pam_len, mb = struct.unpack("=6i", head[4:28])
        pam = fh.read(pam_len).decode("ascii")
        fh.seek(header_len)
        leaves = struct.unpack("=i", fh.read(4))[0]
    return {"format": version, "pam_at_start": bool(flags & 1),
            "pamless": bool(flags & 2), "guide_len": guide, "pam_len": pam_len,
            "max_bulges": mb, "pam": pam, "num_leaves": leaves}
```

## 8. CRISPRme changes needed

CRISPRme (crisprme-plus) parses `.bin` files itself and mirrors `buildTST`'s rules
in Python, so it must learn format 2 before it can use CRISPRitz 2.9.2 indexes. As
of `a5379e0`:

1. **`PostProcess/crispritz_bins.py:bin_problem`** reads the first int as
   `numLeaves` and walks the leaf records by their `'0'`/`'_'` marker. It must branch
   on the magic: for format 2 skip `header_len`, take `numLeaves` there, read each
   record as `int32 ws`, the flank bytes, `lflags` and an `int32 next` iff bit 7, and
   check the trailer. (A format-2 file cut short fails the trailer test on its own.)
2. **`_leaves`** (used by `missing_bins` to see whether a contig's last part is full)
   must read `numLeaves` at `header_len` for format 2.
3. **`holds_target` / `count_targets`** reproduce 2.9.1's site rule (window inside
   the contig, no N outside the PAM, index 0 as minus strand, the short scan bound).
   For an index built by 2.9.2 they must use section 3's rule instead, or a contig
   whose only sites lie at a boundary is reported as "holds a target but has no
   .bin" (B3/B4 add such sites) and a full-part count is off. Choose the rule by the
   engine version that built the index (`crispritz.py version` >= 2.9.2) or, when a
   `.bin` exists, by its header.
4. **Engine/index pairing.** An index built by 2.9.2 cannot be searched by
   CRISPRitz <= 2.9.1. Where CRISPRme resolves or downloads an index, read the format
   of one `.bin` (section 7) and refuse a format-2 index when the pinned engine is
   older, naming both; record `index_format` in the index's `index_meta.json` and
   on the published archive's card.
5. **Boundary rows downstream.** Format 2 returns sites whose pad runs off a contig
   or into an N, and, with `--prune-pam`, sites whose distal flank does. Any step
   that cuts a fixed window around a target (scoring context, variant lookup at
   `pos + offset`, indel FAKE-contig mapping) must clamp at the contig ends and not
   assume `mb` (or `K`) bases exist beyond the target. The target itself (guide +
   PAM as reported) always lies inside the contig and holds no N in its guide.
6. **Equivalence claims.** With B2-B4, a `--prune-pam` search of a format-2 PAM-less
   index of flank `K` returns the same rows as a dedicated format-2 index for any
   3' PAM `L <= K`, contig ends included, and an index with a larger `mb` returns a
   superset of a smaller one's rows. These claims do not hold for format-1 indexes
   (the audit's counts for hg38 stand for them).

## 9. searchBruteForce (B2 only)

The no-bulge scan has no index; it gets the scan bound of B2 on both PAM sides (a
3' PAM ending on the last base; every reverse-strand 5' PAM whose window ends at the
contig end) and returns early on a contig shorter than one target, where the
unsigned bound used to wrap. Its other rules are already those of section 3.1 for
`mb = 0`: the window must lie inside the contig and an N in the guide region drops
it (`currentmissmatch += 100`).

## 10. Cost (measured)

`crispritz.py index-genome <name> <dir with one chromosome> <PAM file> -bMax 2 -th 8`, hg38
chromosome FASTA (chr22 md5 319cba9f..., chr1 md5 9465e0f0...), CRISPRitz 2.9.1 (the
shipped binaries, md5-identical to a compile of tag v2.9.1) against this branch, each
timed by a wrapper that records wall time, user CPU and the child's peak RSS. arm64:
Docker on Apple Silicon (24 CPUs, g++ 14.2 -O3), two rounds alternating the versions.

| index | engine | wall s (round 1 / 2) | user CPU s | peak RSS GB | .bin bytes | leaves |
|-------|--------|----------------------|------------|-------------|------------|--------|
| chr22 NGG, `-bMax 2` | 2.9.1 | 6.1 / 5.6 | 11.1 / 11.1 | 1.45 | 129,640,072 | 5,295,027 |
| | 2.9.2 | 6.1 / 5.8 | 11.1 / 11.4 | 1.46 / 1.47 | 129,639,385 | 5,295,032 |
| chr22 PAM-less K=8 | 2.9.1 | 108.4 / 94.9 | 277.4 / 268.0 | 18.51 | 1,864,232,876 | 78,317,538 |
| | 2.9.2 | 107.0 / 136.2 | 274.5 / 300.1 | 18.88 | 1,864,235,403 | 78,317,730 |
| chr1 NGG, `-bMax 2` | 2.9.1 | 31.9 / 38.5 | 74.9 / 74.4 | 5.53 / 5.54 | 579,345,896 | 24,957,801 |
| | 2.9.2 | 28.8 / 32.2 | 71.4 / 70.9 | 5.63 | 579,348,060 | 24,957,823 |

amd64: ml008 (255 cores shared with other users, load average 250-500 during the runs;
g++ 14.2 -O3, the rc281 image's toolchain), the 2.9.2 build after the 2.9.1 one:

| index | engine | wall s | user CPU s | peak RSS GB | .bin bytes | leaves |
|-------|--------|--------|------------|-------------|------------|--------|
| chr1 PAM-less K=8 | 2.9.1 | 1,022.7 | 2,890.0 | 101.01 | 10,468,103,808 | 460,955,944 |
| | 2.9.2 | 1,093.6 | 2,819.2 | 102.86 | 10,468,113,334 | 460,956,484 |

* Build time is unchanged within the run-to-run noise of a shared host.
* Peak memory grows by 4 bytes a site during the build (a site is now a window start
  plus its shape, 8 bytes, instead of a signed int): +0.1 GB on chr1 NGG, +0.37 GB
  on chr22 K=8, +1.85 GB on chr1 K=8 (about 2% each). searchTST's in-memory leaf
  stays at 32 bytes.
* Index size changes by bytes: a leaf record keeps its size, each `.bin` gains 32 + K
  bytes of header and trailer, and the added sites are few. The sizes are reproducible
  (identical in both rounds).
* The added leaves are the boundary sites: chr22 NGG +5, chr1 NGG +22, chr22 K=8
  +192, chr1 K=8 +540. Read back from the NGG `.bin` files of chr22 and chr1, every
  format-1 leaf (window start, shape, stored PAM) is present in the format-2 index of
  the same chromosome, all with the whole pad; each added leaf has `pad_valid < 2`,
  and the chr1 ones match an independent scan of hg38 for perfect-match sites 0-1 nt
  from an N run (e.g. chr1 `-` 16799140 d=0, `-` 121782408 d=1, `+` 586008 d=0).
* Search time is unchanged: see the benchmark in section 11.

## 11. Tests

1. **Boundary suite** (`test/run_boundary_test.sh`, `test/boundary/`): synthetic
   contigs with planted perfect, 1-4-mismatch and 1-bulge sites at every distance
   `d = 0..4` from a contig start, a contig end and an N run, on both strands, for a
   3' PAM (NGG, NRG, NG, NNGRRT) and a 5' PAM (TTTV), with the PAM 0..9 bases from a
   contig end or an N run; dedicated indexes at `-bMax` 0-3; PAM-less K=3, K=8 and
   5' K=4 indexes with and without `--prune-pam`; searches at 0-3 DNA/RNA bulges.
   Every indexed search, run through `crispritz.py`, must equal
   `test/boundary/oracle.py` -- a brute-force aligner over the FASTA that enumerates
   the sites of section 3 and aligns with the engine's alignment rules, sharing no
   code with it -- as a set of rows (141 searches, 5,889,164 rows), and every
   0-bulge search must equal `crispritz.py search` without bulges, i.e.
   `searchBruteForce` (32). The 2.9.1 binaries fail 142 of those comparisons; each
   new mechanism, disabled on its own (pad limit, `FLANK_OFF` skip, brute-force
   bound, guide-only N rule), fails 10 to 52; so does the oracle without its model of
   the engine's leaf node.
2. **Format-1 compatibility** (`test/run_format1_test.sh`): indexes built by the 2.9.1
   `buildTST` (committed fixtures, byte-reproducible on arm64 and amd64): rows away
   from index 0 byte-identical to the 2.9.1 `searchTST`'s, in order; index-0 rows equal
   to the oracle's.
3. **Index format** (`test/run_index_format_test.sh`): header and trailer, `index-info`,
   a directory mixing both formats, `max_bulges` from the header, and the refusals.
4. **Existing tests** (`run_golden_test.sh`, `run_minimal_test.sh`,
   `run_pam_filter_test.sh`, `run_flank_n_test.sh`, `run_short_contig_test.sh`) pass;
   the golden outputs are byte-identical (their genome holds no boundary site).
   `run_short_contig_test.sh` now also requires a contig exactly one target long to
   be indexed.
5. **AddressSanitizer**: all three binaries, over every test above, on arm64 and amd64.
6. **Platforms**: arm64 (Docker on Apple Silicon: g++ 14.2, and `make -f
   Makefile_conda` in `gcc:12`) and amd64 (ml008, g++ 14.2); every test passes on each.
7. **hg38**: see section 10 for build cost and leaf-level comparison, and:
   * the shipped format-1 `NGG_3_hg38` (all 455 contigs, 501 trees, not rebuilt),
     searched with the guides of the 11 contigs that open with a reverse NGG
     (`-mm 0 -bDNA 1 -bRNA 1`, 32 threads): 63,095 rows from both engines; the 11
     contigs' rows move from the + strand at 2-4 to the minus strand at 0 under 2.9.2
     (same rows otherwise); every other contig's rows are identical (in order for the
     132 single-tree contigs; the 6 multi-tree chromosomes interleave their trees in
     thread-completion order, the same rows).
   * 246 NGG contig-end sites on 171 alt/random/Un contigs, each searched with its own
     perfect guide: a dedicated format-1 index finds 235 with 2.9.1 (11 come out on the
     + strand) and 246 with 2.9.2; a format-1 K=8 index with `--prune-pam` finds 0; a
     format-2 K=8 index with `--prune-pam` finds 246, and its rows equal the format-2
     dedicated index's (368,586 rows). 2.9.2 on the format-1 dedicated index changes the
     rows of 12 contigs only: the 11 move to the minus strand at base 0, and 6 rows of
     the builder's index-0 artefact on chrUn_KI270539v1 (reported on the + strand with
     the PAM `GGC`) are dropped; the other 159 contigs are byte-identical.
   * chr22 and chr1, 50 guides (45 whose perfect match lies 0-2 nt from an N run on
     chr1/chr22, and TRAC, EMX1, VEGFA site 3, HEK4, FANCF), `-mm 4 -bDNA 1 -bRNA 1`:
     - the format-1 index gives byte-identical output under 2.9.1 and 2.9.2: chr22
       5,232,494 rows, chr1 14,358,326 rows (each of chr1's 5 trees byte-identical,
       searched alone; with 16 threads the trees' rows interleave in completion order);
     - the format-2 index returns every row of the format-1 index and 1,827 (chr22) and
       5,512 (chr1) more, every one at a site whose pad runs into an N run: the 5 and 22
       added sites;
     - the 0-bulge search of the format-2 index equals the brute force (chr22 32,150
       rows, chr1 81,522); the 2.9.1 index misses 10 and 24 of those rows, including 5
       and 22 perfect matches;
     - `--prune-pam` NGG off the format-2 K=8 index equals the format-2 dedicated NGG
       index on chr22 (5,234,321 rows).
   * search time, chr22 NGG, the 50 guides, three alternating rounds (wall s): 2.9.1 on
     the format-1 index 81.5 / 87.4 / 80.7, 2.9.2 on it 70.3 / 82.9 / 78.5 (output
     byte-identical each round), 2.9.2 on the format-2 index 79.6 / 81.3 / 80.3.

## 12. Not changed

* A 5' PAM shorter than the stored flank is still refused under `--prune-pam`.
* A guide-region N still discards the site, even for an alignment with RNA bulges
  that would not reach the N.
* Coordinates are still 32-bit ints (contigs up to 2^31 - 1 - (mb + K) bases).
* searchTST still tries no RNA bulge on the last base of a stored string (a leaf node
  of the tree), so an alignment that would leave the last guide base unpaired after
  consuming all but one of the `G + mb` stored bases is not reported, and whether such
  an alignment is reported depends on the index's `-bMax`. This holds for interior
  sites as much as boundary ones, is part of the golden outputs, and is modelled by
  the oracle (without that model it predicts 352-1,786 extra rows per bulged search).
* `main.cpp` (the old single-threaded builder, not built by any recipe) is left as
  is; it does not compile against the current `pamGenerator.cpp` and never did since
  the Aho-Corasick code was commented out.
