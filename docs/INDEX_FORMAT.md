# CRISPRitz index format and target windows (format 2, CRISPRitz 2.10.0)

This document specifies what a CRISPRitz TST index holds, how it is laid out on
disk, and how the search reads it. It covers the boundary/window fixes of CRISPRitz
2.10.0 (audit ids IDX-13, IDX-14, IDX-15 and IDX-4), the index-format version stamp
(B5) that those fixes need, because two of them change what an index holds, and Q,
the fifth symbol of a SNP position whose reference base and ALTs cover all four bases
(section 4.1).

Status: implemented (branches `feat/engine-boundaries` and `feat/multiallelic-enrich`,
merged into CRISPRitz 2.10.0, which adds Q). Section 10 holds the measurements;
section 11 the tests.

| id | audit | what goes wrong in 2.9.1 | fix | needs a rebuild |
|----|-------|--------------------------|-----|-----------------|
| B1 | IDX-14 | a leaf with index 0 is built as a minus-strand window (3' PAM) but read as plus strand, so a reverse-strand PAM at a contig's first base is reported on the + strand at base `max_bulges` | search side: read index 0 as PAM-left for a 3' PAM; drop the index-0 leaf the builder made from a plus-strand window | no (works on format-1 indexes) |
| B2 | IDX-15 | a PAM that ends on a contig's last base is never found, by any index or by the brute-force scan | scan bound `nt <= len - pam_limit` on both PAM sides, in `buildTST` and `searchBruteForce` | for indexes yes; brute force no |
| B3 | IDX-13 | a site whose guide lies within `max_bulges` nt of a contig start/end or an N run is dropped, even a perfect match: the whole pad+guide+PAM window must fit and be free of N | emit a site when its guide+PAM core fits; store the pad positions that are off the contig or on an N as non-consumable | yes (format 2) |
| B4 | IDX-4 | `--prune-pam` with flank `K` wider than the PAM `L` loses + strand PAMs with 1..K-L bases to the contig end and - strand PAMs 0..K-L-1 from the start, and reports a - strand PAM exactly K-L from the start on the + strand | store a PAM-less window whose distal flank runs off the contig (flank bases there stored as N); store the window start as a signed coordinate with the strand in its own bit | yes (format 2) |
| B5 | - | an index does not say what it holds | a header with a magic number, a format version and the build parameters; old indexes keep loading and searching as with 2.9.1 | - |

## 1. Terms

* `G` the guide length the index was built for (the N's of the build PAM file
  before a 3' PAM or after a 5' one: 20 for SpCas9), `K` stored PAM/flank width (the
  build PAM file's declared length), `L` the search PAM's declared length (`L = K`
  except under `--prune-pam`), `mb` the index's `max_bulges` (`-bMax`).
* `Gs` the search guide length (the search PAM file's), and `off = G - Gs` (2.9.1's
  `offset_guide_len`). The search does not require `Gs = G`; section 5.1 says what a
  search guide of another length finds. Elsewhere this document takes `Gs = G`
  (`off = 0`) unless it says otherwise.
* A **site** is one PAM occurrence on one strand. Its **window** is the stretch of
  the contig the index stores for it: the guide region (`G` bases), the PAM/flank
  (`K` bases) on its PAM side, and a **pad** of `mb` bases on the far (distal) side
  of the guide, which only alignments that consume more than `G` bases reach: DNA
  bulges, or a search guide longer than `G`.
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
  window the leftmost aligned base, `ws + mb + off - dna + rna` (`dna`/`rna` =
  bulges used). Cluster Position: the same as Position for PAM-left windows, the
  start of a `Gs`-base guide, `ws + mb + off`, for PAM-right windows.
* An alignment **consumes** `Gs + dna - rna` stored bases, read from the PAM outward:
  the `G` of the guide region first, then `Gs + dna - rna - G = dna - rna - off` pad
  bases when that is positive.

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

**N is a gap.** A genomic `N`, like a position off the contig, matches nothing: no
guide base, no PAM letter (not even the search PAM's `N`), and it ends the pad. It is
never "any base" and the enricher never writes one for a variant (it keeps the
reference's `N` exactly, even where a VCF lists ALTs at one). The base that matches
everything is `Q` (section 4.1), the symbol of a position where the reference base and
the SNP ALTs cover A, C, G and T: it counts as a base everywhere above (in the guide
region, the pad and the flank), and a site, its pad and its PAM may hold it.

### 3.2 The pad

`pad_valid` (0..`mb`) is the number of pad positions, counted from the guide
outward, before the first one that is off the contig or an N. An alignment consumes
`Gs + dna - rna` stored bases (section 1), of which `Gs + dna - rna - G` (that is
`dna - rna - off`; `dna - rna` when the search guide has the index's length) are pad
bases. If that is more than `pad_valid`, the alignment needs a base that is not there
and is not reported. A search guide longer than the index's (`Gs > G`) reaches into
the pad without any bulge, so at a site whose pad runs off the contig or into an N
after fewer than `Gs - G` bases only alignments with enough RNA bulges are found; a
shorter one (`Gs < G`) needs `G - Gs` DNA bulges more than RNA bulges to reach the
pad. Pad positions from `pad_valid` on are stored
in the tree as `A`: no reported alignment can depend on them (an alignment reaching
them is rejected by the `pad_valid` test at leaf emission), so any symbol serves, and
the 14 nibble values of the tree's alphabet stay as they are. For every site of 2.9.1
(window inside the contig, no N outside the PAM), `pad_valid = mb` and the stored
string is unchanged.

So every base a reported alignment consumes is a base of the contig and not an N:
the guide region's `G` by the site rule, the pad's by this test. Every reported
target (its guide part and its PAM) lies inside the contig, at any search guide
length (the boundary suite checks every row it gets against the contig's bases,
section 11).

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
                          PAM; bit 2: FOUR_ALLELE, the file holds Q (section 4.1);
                          other bits 0 (reserved, must be 0: a reader refuses a
                          file with a bit it does not know)
  16   4     guide_len    G (the body's len_guide must equal it)
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

### 4.1 Q, the four-allele symbol

`add-variants` writes `Q` where a position's reference base (A, C, G or T) and its SNP
ALTs cover all four bases (docs/MULTIALLELIC_ENRICHMENT.md, section 2.2). Its IUPAC
letter would be `N`, which section 3 reads as a gap, so up to the four-base rule of the
multi-allelic work such a position removed every target through it from the index.
`buildTST`, `searchTST` and `searchBruteForce` read `Q` as the mask `1111`: it matches
every guide base and every PAM letter, so it is never a mismatch.

Every 4-bit slot of the tree already uses 15 of its 16 values: a root, lo or hi slot
holds `0` (no kid) or one of the 14 IUPAC letters (1..14, the letter's A|C|G|T mask),
and an eq slot holds one of the 14 letters or `15`, the end of a string (followed by the
int32 leaf index). So `Q` takes the value its slot never uses:

| where | `Q` is written as |
|-------|-------------------|
| a root, lo or hi slot of the tree | `15` (`1111`, its mask) |
| an eq slot of the tree | `0` |
| a flank nibble (a leaf's stored PAM/flank) | `15` |

Only a file whose header has the `FOUR_ALLELE` flag (bit 2) holds these values; the
builder sets it on a `.bin` when one of its leaves stores `Q` in its guide region, its
pad or its flank, and a `.bin` without `Q` is written byte for byte as without the fifth
symbol. A reader that knows bit 2 decodes the tree slot by slot (in such a file a `1111`
high nibble is not always an end of string, and its byte's low nibble is not always
padding); a reader that does not refuses the file by section 5's rule on reserved bits
(the format-2 reader of `feat/engine-boundaries`, before `Q`, prints `has a malformed
index header (version 2, flags 4, ...)`), and CRISPRitz 2.9.1 and earlier refuse every
format-2 file. `buildTST` prints on stdout how many targets hold `Q` and how many files it
flagged (`FOUR-ALLELE (Q): ...`), `searchTST` how many of the files it reads are flagged
(`INDEX FORMAT: 2 x42 (40 holding Q, the four-allele symbol)`), and `crispritz.py
index-info` reports `four_allele` for each file (section 7).

A reference FASTA must not hold `Q`: `add-variants` refuses one, in a chromosome it
enriches and in one it copies (`ERROR: <chr> holds 'Q' at <pos>: ...`). `buildTST` cannot
tell a reference from an enriched genome and indexes a `Q` wherever it finds one.

## 5. How the search reads an index

`loadTST` reads the first 4 bytes of each `.bin`. If they are the magic it parses
the header and reads the body as format 2; otherwise they are `numLeaves` and the
file is format 1. Each `.bin` is decoded on its own, so a directory mixing formats
still searches correctly (though `buildTST` never makes one).

The search refuses, with an error on stderr and exit 1, before it searches (the leaf
count: as it loads the file):

* a file it cannot open, or too short to hold the magic (fewer than 4 bytes);
* a format version it does not know (> 2);
* a malformed header: the magic followed by a version other than 2, a flag bit other
  than bits 0-2 set, `pam_len` outside 1..1000, `header_len` other than
  `28 + pam_len`, `guide_len` below 1, or `max_bulges` outside 0..31;
* a format-2 file cut short: inside its header, or not ending in the trailer;
* a header `guide_len` that differs from the body's `len_guide` (malformed);
* a header PAM other than the file name's PAM prefix;
* a leaf count below 1, in either format (2.9.1 refuses that too);
* a leaf count the file cannot hold: a leaf record takes at least `4 + ceil(K/2) + 1`
  bytes, so a file shorter than `header_len + 12 + numLeaves x (4 + ceil(K/2) + 1)`
  (format 1: `header_len` 0, `K` the length of the PAM it is named for) is damaged; the
  search used to allocate that many leaves and fail in `std::bad_alloc`, or read past
  the file's end;
* a directory whose `.bin` files are named for more than one PAM (the prefix of
  `<PAM>_<contig>_<part>.bin`): the search reads every file of a directory as an index
  of the PAM that names the first it lists, so such a directory is not one index. A
  macOS `._` (AppleDouble) file copied along with an index is one; whatever order the
  directory lists them in, the search refuses it before reading a file (until the fix,
  only when the `._` file came first: listed after the index's own files, it was read as
  a tree of 118,883,584 leaves and the search failed after 22.5 s);
* and, depending on the search: a search PAM that does not fit the stored PAM width
  (the file name's PAM prefix; for format 2 also `pam_len`): it must be as wide, or
  under `--prune-pam` a narrower 3' PAM, N-padded to it (section 12); or a PAM side
  (3'/5') other than the search PAM's.

`crispritz.py index-info` makes every one of these checks that does not depend on the
search (section 7). Loading a file, the search also refuses a node count below 1, and
reports on stderr a flank nibble that is no IUPAC code (both as 2.9.1 did); those lie in
or past the leaves, which `index-info` does not read.

If `max_bulges` on the command line differs from a format-2 header's, the header
wins for decoding coordinates and a note goes to stdout (stderr is fatal in
CRISPRme). Before searching, `searchTST` prints the formats it found on stdout, e.g.
`INDEX FORMAT: 2 x24`, `INDEX FORMAT: 1 (CRISPRitz <= 2.9.1) x501`, or, with files
flagged `FOUR_ALLELE`, `INDEX FORMAT: 2 x24 (3 holding Q, the four-allele symbol)`.

The guide length the index was built for, `G` (the header's `guide_len`; format 1:
the body's `len_guide`), is read for each `.bin` and kept per search thread with its
`off = G - Gs`, as is `mb`. A search guide length `Gs` other than `G` is not refused
(2.9.1 searched it, and CRISPRme sends such searches, section 8): it is searched with
`off` in the coordinates (section 1) and in the pad test below, and `searchTST` says on
stdout, once per `G`, what that changes, e.g.

```
NOTE: the index was built for 20-nt guides and the search guides are 21 nt; it holds
22 bases a site, so an alignment is found only if its DNA bulges minus RNA bulges are
at most 1 (an index built for 21-nt guides with max_bulges 2 allows 2). See
docs/INDEX_FORMAT.md, section 5.1.
```

Section 5.1 gives the rules and how the rows differ from those of an index built for
`Gs`.

Per leaf the search derives `PAM_LEFT`, `ws`, `pad_valid` and `FLANK_OFF`:

| | format 1 | format 2 |
|-|----------|----------|
| `PAM_LEFT` | `guideIndex < 0`, or `guideIndex == 0` and the PAM is 3' (B1) | `lflags` bit 6 |
| `ws` | `abs(guideIndex)` | `ws` |
| `pad_valid` | no limit (255: 2.9.1 kept only whole windows, and an alignment cannot consume more than the stored string) | `lflags` bits 0-4 |
| `FLANK_OFF` | 0 | `lflags` bit 5 |
| skip | an index-0 leaf whose flank does not match the build PAM (B1) | - |

At emission a leaf is skipped when `--prune-pam` is on and its flank fails the
search PAM (unchanged), when the index is PAM-less, `--prune-pam` is off and the
leaf is `FLANK_OFF`, or when the alignment used more pad bases than `pad_valid`:
`dna - rna - off > pad_valid` (section 3.2). Coordinates are then computed as in
section 1, with `mb` from the header for format 2 and from the command line for
format 1 (as 2.9.1 does), and `off` from the `.bin` the leaf came from.

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

### 5.1 A search guide of another length than the index's

`searchTST` takes the guide length from the search PAM file (`Gs`) and the index's from
the `.bin` (`G`), and searches when they differ, as 2.9.1 did. CRISPRme does this: its
`21bp-NNGRRT-SaCas9.txt` (21-nt guides) is served under `--prune-pam` by a PAM-less
index built for 20-nt guides, and `23bp-TTTV-Cas12a.txt` (23-nt) by a 20-nt TTTV or
NNNN index (section 8). The rules, for both formats:

1. **Sites are the index's**: a `G`-base guide region inside the contig with no N,
   and the PAM of section 3.1. When `Gs < G`, a target whose `Gs`-base guide fits but
   whose `G`-base region does not -- its far end within `G - Gs` nt of a contig end or
   an N run -- is not in the index and is not found.
2. **An alignment consumes `Gs + dna - rna` of the `G + mb` stored bases**, so it is
   found only if `dna - rna <= G + mb - Gs`: a longer guide has fewer DNA bulges to
   spend (with `Gs > G + mb`, none is found without an RNA bulge).
3. **The pad test** (section 3.2): the `dna - rna - off` pad bases it consumes must
   all exist. So a longer guide's extra bases are never paired with a base off the
   contig or on an N, and every reported row lies inside the contig.
4. **The leaf rule** of section 12 applies where the stored string ends, `G + mb`
   bases from the PAM.
5. **Coordinates** carry `off` (section 1).

The rows are therefore exact for what the index holds, and on a format-2 index a
0-bulge search with a longer guide (`G < Gs <= G + mb`) returns what
`searchBruteForce` returns (the boundary suite checks both). But they can differ from
the rows of an index built for `Gs`: not without bulges when `Gs > G`, but with
bulges, and for a shorter guide at contig ends and N runs. Measured with
`crispritz.py` on the boundary suite's genome (`test/boundary/make_genome.py`, 346
contigs a PAM side), rows only in one of the two:

| search guides | the index for `G` | the index for `Gs` | -mm/-bDNA/-bRNA | rows (`G` / `Gs`) | only `G`'s | only `Gs`'s |
|---|---|---|---|---|---|---|
| 21-nt NNGRRT | K=8 PAM-less, 20 nt, `-bMax 2`, `--prune-pam` | NNGRRT, 21 nt, `-bMax 2` | 4/0/0 | 290 / 290 | 0 | 0 |
| | | | 1/1/0 | 274 / 274 | 0 | 0 |
| | | | 1/0/1 | 324 / 272 | 52 | 0 |
| | | | 1/2/2 | 222,162 / 246,728 | 3,912 | 28,478 |
| 23-nt TTTV | TTTV, 20 nt, `-bMax 3` | TTTV, 23 nt, `-bMax 3` | 4/0/0 | 122 / 122 | 0 | 0 |
| | | | 1/1/0 | 16 / 76 | 0 | 60 |
| | | | 1/0/1 | 102 / 100 | 26 | 24 |
| | | | 1/2/2 | 67,280 / 91,196 | 1,920 | 25,836 |
| 19-nt NGG | NGG, 20 nt, `-bMax 2` | NGG, 19 nt, `-bMax 2` | 4/0/0 | 492 / 508 | 0 | 16 |
| | | | 1/1/0 | 412 / 416 | 0 | 4 |
| | | | 1/0/1 | 640 / 684 | 0 | 44 |
| | | | 1/2/2 | 432,068 / 448,670 | 690 | 17,292 |

Why they differ, row by row:

* only the `G` index's, `Gs > G`: RNA-bulged alignments that stop short of the guide
  region's far end at a site where it runs off the contig or into an N (the `Gs`
  index needs its whole `Gs`-base region and does not hold the site). All at boundary
  sites.
* only the `Gs` index's, `Gs > G`: alignments with more DNA bulges than rule 2
  allows (the 60 rows at 1/1/0: no DNA bulge fits 23-nt guides into 20 + 3 stored
  bases; 166 rows with 2 more DNA than RNA bulges at 1/2/2 for 21-nt guides), and
  the leaf rule: an alignment that leaves its last guide base unpaired after
  `G + mb - 1` bases is not tried where the `G` index's string ends, but is by the
  longer `Gs` index. These occur at interior sites too (1,650 of the 28,478; 1,482
  of the 25,836).
* only the `Gs` index's, `Gs < G`: rule 1's sites (16 perfect and mismatched
  targets at 4/0/0, the suite's guides cut by a contig end or an N run), all at
  boundary sites.
* only the `G` index's, `Gs < G`: the leaf rule the other way round (690 rows at
  1/2/2, 40 of them at interior sites).

What does hold for any `Gs`: a `--prune-pam` search of a PAM-less index returns the
same rows as a dedicated index **built for the same `G` and `mb`** (rows only in one:
0 at each of the four settings of the table's 21-nt NNGRRT lines, comparing the K=8
index with a 20-nt NNGRRT index; the boundary suite checks it for `Gs` = 19, 20, 21
and 23, section 11).

For the rows of an index built for the search's guide length, search one built for
it: `index-genome` with the search's PAM file (its N's are the guide length).

## 6. Compatibility

| | format-1 index (built by <= 2.9.1) | format-2 index (built by 2.10.0) |
|-|------------------------------------|----------------------------------|
| CRISPRitz <= 2.9.0 | as before | searchTST aborts: `std::length_error` (`vector::_M_default_append`) on stderr, exit 134, empty targets file (tested with the 2.9.0 binary of the CRISPRme 2.8.0 image; its `crispritz.py` does not check exit statuses, CRISPRme fails the stage on the stderr line) |
| CRISPRitz 2.9.1 | as before | refused: `'...' is not a valid CRISPRitz index (.bin) file (numLeaves=-28747197)`, exit 1 (tested) |
| CRISPRitz 2.10.0 searchTST | byte-identical to 2.9.1 except index-0 rows (B1) | full B1-B4 semantics, and `Q` (section 4.1) |
| CRISPRitz 2.10.0 searchBruteForce | - (no index) | - ; B2 applies at once |

Byte-identical means the same `.targets.txt` bytes, row order included; tested on
the repo's fixtures and on indexes built by the 2.9.1 binary (section 11).
Users' existing indexes keep working unchanged. They keep missing the sites of
B2-B4 until rebuilt; B1 corrects their index-0 rows at once.

A format-2 index of a genome holds every site the format-1 index of the same
genome and parameters holds, with the same string, flank and coordinates, plus the
boundary sites (measured on hg38 chr22 and chr1, section 10). A 2.10.0 search on it
returns every row a 2.10.0 search on the format-1 index returns, plus rows at sites
whose guide lies within `mb` nt of a contig end or an N run (B3, including a + strand
3' site whose guide starts at exactly `mb`), whose PAM lies within `K - L` nt of a
contig end under `--prune-pam` (B4), or whose PAM ends on a contig's last base (B2).

## 7. Reading the stamp without loading the tree

Any reader (CRISPRme, scripts) can classify an index from the first bytes of its `.bin`
files, and must refuse what `searchTST` refuses before it relies on them. The function
below is the reference implementation; it is `crispritz.py`'s own, and
`test/run_index_format_test.sh` checks that this copy returns what `crispritz.py
index-info` does for every file of its corpus (good files of both formats, and 20 damaged
in one way each). `crispritz.py index-info <index_dir>` prints its result for every `.bin`
of a directory as JSON (`four_allele`: the file holds `Q`, section 4.1), with
`writes_format`, the format this CRISPRitz builds (its `INDEX_FORMAT`, section 8), and
`pams`, the PAMs the files are named for, and exits 1 when the directory holds no `.bin`,
a `.bin` with an `error` entry, or `.bin` files named for more than one PAM, which
`searchTST` refuses too (section 5). A macOS `._` (AppleDouble) file copied along with an
index is one: its first int, read as a format-1 leaf count, is positive, but it is named
for `.`, and its 118,883,584 leaves cannot fit in its 163 bytes.

An `error` entry means `searchTST` refuses the file whatever the search: the function
makes every check of section 5 that does not depend on the search, in `searchTST`'s order,
reading the header, the two ints that open the body, the last 4 bytes and the file's
size. (A format-1 file of 4 to 7 bytes, which it calls too short, `searchTST` refuses when
it reads the node count.) A format-1 file has no end marker, so one cut short after its
first 8 bytes is not seen unless its leaf count cannot fit; nor is damage inside the
leaves or the tree of either format, which `searchTST` finds only as it loads them, if at
all (CRISPRme's `bin_problem` walks the leaves, section 8).

```python
CZI_MAGIC = b"CZI\xfe"


def read_index_header(path):
    """The index format and build parameters of the CRISPRitz .bin ``path``, read from
    its header, the two ints that open its body and its last 4 bytes, without loading
    the tree (docs/INDEX_FORMAT.md section 7). An "error" entry says why searchTST
    refuses the file: every check searchTST makes on a file's own bytes and name
    before it reads its leaves is made here, in the same order. A format-1 file has no
    end marker, so one cut short after its first 8 bytes is not seen here."""
    import os
    import struct

    def refused(info, why):
        info["error"] = why
        return info

    def damaged(info, header_len, width, size):
        # a leaf record takes at least 4 + ceil(K/2) + 1 bytes (window start, flank
        # nibbles, marker or flags byte), so a leaf count the file cannot hold is damage
        need = header_len + 12 + info["num_leaves"] * (4 + (width + 1) // 2 + 1)
        if need > size:
            return refused(info, "damaged: its %d leaves need at least %d bytes, the file holds %d"
                           % (info["num_leaves"], need, size))
        return info

    # the PAM the file is named for, <PAM>_<contig>_<part>.bin: searchTST reads every
    # .bin of a directory at the width of that PAM
    named = os.fsencode(os.path.basename(path)).split(b"_", 1)[0]
    try:
        fh = open(path, "rb")
    except OSError as exc:
        return refused({"format": None}, "cannot be opened (%s)" % exc.strerror)
    with fh:
        size = os.fstat(fh.fileno()).st_size
        head = fh.read(28)
        if head[:4] != CZI_MAGIC:
            # format 1 (CRISPRitz <= 2.9.1): no header, the file opens with numLeaves
            # and len_guide
            if len(head) < 8:
                return refused({"format": None}, "too short to be a CRISPRitz index (.bin) file")
            leaves, guide = struct.unpack("=ii", head[:8])
            info = {"format": 1, "num_leaves": leaves, "guide_len": guide}
            if leaves <= 0:
                return refused(info, "not a valid CRISPRitz index (.bin) file (numLeaves=%d)" % leaves)
            return damaged(info, 0, len(named), size)
        if len(head) < 28:
            return refused({"format": None}, "cut-short index header")
        version, header_len, flags, guide, pam_len, mb = struct.unpack("=6i", head[4:28])
        info = {
            "format": version if version >= 2 else None,
            "pam_at_start": bool(flags & 1),
            "pamless": bool(flags & 2),
            "four_allele": bool(flags & 4),
            "guide_len": guide,
            "pam_len": pam_len,
            "max_bulges": mb,
        }
        if version > 2:
            return refused(info, "index format %d, newer than this CRISPRitz reads (formats 1 and "
                                 "2): update CRISPRitz, or rebuild the index with this one" % version)
        # bit 0: PAM at the 5' end, bit 1: PAM-less, bit 2: holds Q, the four-allele
        # symbol; max_bulges fits pad_valid's 5 bits
        if (version != 2 or flags & ~7 or not 1 <= pam_len <= 1000 or header_len != 28 + pam_len
                or guide < 1 or not 0 <= mb <= 31):
            return refused(info, "malformed index header (version %d, flags %d, PAM length %d, "
                                 "header length %d, guide length %d, max bulges %d)"
                           % (version, flags, pam_len, header_len, guide, mb))
        pam = fh.read(pam_len)
        info["pam"] = pam.decode("ascii", "replace")
        if len(pam) < pam_len:
            return refused(info, "cut-short index header")
        fh.seek(-4, os.SEEK_END)
        info["complete"] = fh.read(4) == CZI_MAGIC
        if not info["complete"]:
            return refused(info, "cut short: it does not end in the end-of-index marker")
        fh.seek(header_len)
        body = fh.read(8)
    if len(body) < 8:
        return refused(info, "cut short: no leaf count and guide length after its header")
    leaves, body_guide = struct.unpack("=ii", body)
    info["num_leaves"] = leaves
    if body_guide != guide:
        return refused(info, "malformed index: its header gives guide length %d, its body %d"
                       % (guide, body_guide))
    # searchTST reads a format-2 file only if its header PAM is its name's prefix
    if pam != named:
        return refused(info, "holds the PAM/flank region '%s' but is named for '%s'"
                       % (info["pam"], os.fsdecode(named)))
    if leaves <= 0:
        return refused(info, "not a valid CRISPRitz index (.bin) file (numLeaves=%d)" % leaves)
    return damaged(info, header_len, pam_len, size)
```

## 8. CRISPRme changes needed

CRISPRme (crisprme-plus) parses `.bin` files itself and mirrors `buildTST`'s rules
in Python, so it must learn format 2 before it can use CRISPRitz 2.10.0 indexes. As
of `a5379e0`:

1. **`PostProcess/crispritz_bins.py:bin_problem`** reads the first int as
   `numLeaves` and walks the leaf records by their `'0'`/`'_'` marker, so it calls
   every format-2 `.bin` "not a tree (numLeaves=-28747197)". It must branch on the
   magic: for format 2 make section 7's checks, take `numLeaves` at `header_len`, read
   each record as `int32 ws`, the flank bytes, `lflags` and an `int32 next` iff bit 7,
   and check the trailer. (A format-2 file cut short fails the trailer test on its own.)
2. **`_leaves`** (used by `missing_bins` to see whether a contig's last part is full)
   must read `numLeaves` at `header_len` for format 2. As it is, it returns the magic
   (-28,747,197), never `TARGETS_A_TREE`, so on a format-2 index the full-part check
   never runs.
3. **`holds_target` / `count_targets`** implement 2.9.1's site rule (window inside
   the contig, no N outside the PAM, index 0 as minus strand, the short scan bound).
   A format-2 index holds more sites (section 3: B2-B4), so `missing_bins` must apply
   the rule of the format the index was built in:
   * **Per index, from the header of its `.bin` files** (section 7; `crispritz.py
     index-info`): format 1, the 2.9.1 rule; format 2, section 3's rule. Not from
     `crispritz.py version`: that is the installed engine, not the one that built an
     index already on disk, so with 2.10.0 installed it would apply section 3's rule to
     every 2.9.1-built index (and a build of the branch before the release bumped
     `VERSION` said 2.9.1). `index-genome` removes an index directory before it builds one,
     so a directory whose `.bin` files differ in format was assembled from two builds:
     treat it as incomplete and rebuild it.
   * **A directory with no `.bin`** has no header to read. When CRISPRme's own build
     has just made it (`build-index-only`; `pool_index_indels.py`'s
     `_record_no_targets`), use the format of the engine that ran: `INDEX_FORMAT` in
     the `crispritz.py` on PATH (2 from 2.10.0; read from the file as
     `generate_report._crispritz_version` reads `VERSION`; CRISPRitz <= 2.9.1 has none
     and writes format 1), and record it with the index (`index_meta.json`, item 4; the
     `NO_TARGETS` record). For one built earlier (`pool_search_indels.py`'s
     `_empty_index`), use that record, else the format of the index's other `.bin`
     files (its other sub-indexes), else section 3's rule: on a directory 2.9.1 built,
     its only error is a false "holds a target", which asks for a rebuild, where the
     2.9.1 rule on one 2.10.0 built would pass a lost `.bin` off as an empty result.
   * **What the wrong rule does.** Measured with `crispritz_bins` of `a5379e0` on a
     synthetic boundary genome (885 contigs; NGG `-bMax 2` indexes built by 2.9.1 and
     by this branch with `crispritz.py index-genome`), and with a section-3 rule
     written from section 3.1 put in place of `holds_target` (its site count equals
     the format-2 leaf count on all 885 contigs):
     - The 2.9.1 rule on a format-2 index raises no false report: format 2 holds every
       2.9.1 site (on every contig the format-2 leaf count is at least the format-1
       one, which `count_targets` gives exactly). But it misses lost `.bin` files.
       With the `.bin` deleted of the two contigs whose only site is a boundary site
       (both 23 nt, one target long: `exact_p`, a + strand site whose PAM ends on the
       last base, B2, and whose pad runs off the start, B3; `exact_m`, a - strand site
       whose pad runs off the end, B3), `missing_bins` returns nothing; with section
       3's rule it reports `exact_m`, the first it reaches. `count_targets` undercounts a
       format-2 contig by its boundary sites (467 sites on 361 of the 885 contigs),
       so a lost last part goes unseen when it would hold no more targets than that;
       and with `_leaves` as it is (item 2) the full-part check does not run at all.
     - Section 3's rule on a format-1 index raises false reports: on the complete
       2.9.1-built index it reports `exact_m` as "holds a target but has no .bin", as
       it would every contig whose only sites are boundary sites, so a good 2.9.1
       index fails the check. Its count exceeds the format-1 leaf count on 361 of the
       885 contigs, so a contig whose last part is full can be taken for one that
       lacks a part.
4. **Engine/index pairing.** An index built by 2.10.0 cannot be searched by
   CRISPRitz <= 2.9.1. Where CRISPRme resolves or downloads an index, read the formats
   of its `.bin` files (section 7: `index-info` exits 1 on one that `searchTST` would
   refuse) and refuse a format-2 index when the pinned engine is older, naming both;
   record `index_format` in the index's `index_meta.json` and on the published
   archive's card.
5. **Boundary rows downstream.** Format 2 returns sites whose pad runs off a contig
   or into an N, and, with `--prune-pam`, sites whose distal flank does. Any step
   that cuts a fixed window around a target (scoring context, variant lookup at
   `pos + offset`, indel FAKE-contig mapping) must clamp at the contig ends and not
   assume `mb` (or `K`) bases exist beyond the target. The target itself (guide +
   PAM as reported) always lies inside the contig and holds no N in its guide, at
   any search guide length (section 3.2; the boundary suite checks every row).
6. **Equivalence claims.** With B2-B4, a `--prune-pam` search of a format-2 PAM-less
   index of flank `K` returns the same rows as a dedicated format-2 index **built for
   the same guide length `G` and `mb`**, for any 3' PAM `L <= K` and any search guide
   length, contig ends included; and of two indexes built for the same `G`, the one
   with the larger `mb` returns a superset of the other's rows. Neither holds across
   guide lengths: the rows of an index built for 20-nt guides, searched with 21-nt
   guides, are not those of an index built for 21-nt guides once bulges are searched
   (section 5.1). Neither holds for format-1 indexes either (the audit's counts for
   hg38 stand for them).
7. **Guide-length pairing.** `_resolve_index_dir`
   (`PostProcess/submit_job_automated_new_multiple_vcfs.sh`) and `_pam_serves`
   (`crisprme.py`) pick an index by its PAM name and width and `-bMax` only; nothing
   compares guide lengths, and CRISPRitz searches a guide of any length (section 5.1).
   So `21bp-NNGRRT-SaCas9.txt` and `21bp-NNNRRT-SaCas9KKH.txt` (21-nt guides) are
   served under `--prune-pam` by a PAM-less index built for 20-nt guides, such as the
   K=8 general index, and `23bp-TTTV-Cas12a.txt` (23-nt) by a 20-nt TTTV or NNNN index
   (`23bp-TTTN-Cas12a.txt` by a 20-nt TTTN or NNNN one). The engine handles the
   difference (every row is a real
   target and, for a longer guide without bulges, the rows equal the brute force), but
   with bulges the rows are those of a 20-nt index, which section 5.1 shows differ
   from a 21- or 23-nt index's, at interior sites too, and a shorter guide misses
   targets at contig ends and N runs. CRISPRme must pair the guide length as well:
   read `guide_len` from one `.bin` (section 7; `crispritz.py index-info`) and prefer
   an index whose `guide_len` is the PAM file's guide length, building one for that
   length when none is available, or else use the other knowingly and say in the
   report which guide length the index was built for. `searchTST`'s stdout `NOTE: the
   index was built for <G>-nt guides ...` line (section 5) marks a run where it was
   not.
8. **The four-allele symbol Q** (section 4.1). A genome enriched by CRISPRitz 2.10.0
   holds `Q` where a SNP position's reference base and ALTs cover all four bases, and so
   do the targets a search reports there (in the DNA column, and in a PAM). Every reader
   of an enriched genome, a fake-indel genome overlaid from it, or a target must read `Q`
   as any base, `{A, C, G, T}`, and never as `N` (a gap, which matches nothing): the
   IUPAC tables of CRISPRme 2.8.1 have no `Q` (`pam_filter._IUPAC`/`_TARGET_IUPAC`,
   `overlay_indel_snps._IUPAC_AMBIG`,
   `new_simple_analysis.iupac_code`/`iupac_code_set`/`iupac_nucleotides`,
   `analisi_indels_NNN.iupac_code`, `twopass_find.IUPAC_TO_SET`, among the ~20 modules
   that decode IUPAC letters; docs/MULTIALLELIC_ENRICHMENT.md sections 6.3 and 6.4 say what
   the SNP path and the overlay lose without it), and `crispritz_bins`' site rule must
   count `Q` as a base, which matches every PAM letter. A header flag says which `.bin`
   files hold `Q` (`FOUR_ALLELE`, bit 2; `index-info`'s `four_allele`), and a format-2
   reader must accept that bit. The dictionaries list every ALT of such a position, each in
   its own entry (docs/MULTIALLELIC_ENRICHMENT.md, section 2.3).

## 9. searchBruteForce

The no-bulge scan has no index; it gets the scan bound of B2 on both PAM sides (a
3' PAM ending on the last base; every reverse-strand 5' PAM whose window ends at the
contig end) and returns early on a contig shorter than one target, where the
unsigned bound used to wrap. Its other rules are those of section 3.1 for `mb = 0`:
the window must lie inside the contig, and an N in the guide region drops it whatever
the mismatch budget (until 2.10.0 an N counted 101 mismatches, so a search allowing
more than 100 reported windows holding one), and it reads `Q` as any base.

## 10. Cost (measured)

Measured with the build of `feat/engine-boundaries` (`c286bf8`), before `Q`: 2.10.0 writes
the same bytes for a genome without `Q` (section 4.1), so the reference-genome numbers
stand for it.

`crispritz.py index-genome <name> <dir with one chromosome> <PAM file> -bMax 2 -th 8`, hg38
chromosome FASTA (chr22 md5 319cba9f..., chr1 md5 9465e0f0...), CRISPRitz 2.9.1 (the
shipped binaries, md5-identical to a compile of tag v2.9.1) against this branch, each
timed by a wrapper that records wall time, user CPU and the child's peak RSS. arm64:
Docker on Apple Silicon (24 CPUs, g++ 14.2 -O3), two rounds alternating the versions.

| index | engine | wall s (round 1 / 2) | user CPU s | peak RSS GB | .bin bytes | leaves |
|-------|--------|----------------------|------------|-------------|------------|--------|
| chr22 NGG, `-bMax 2` | 2.9.1 | 6.1 / 5.6 | 11.1 / 11.1 | 1.45 | 129,640,072 | 5,295,027 |
| | 2.10.0 | 6.1 / 5.8 | 11.1 / 11.4 | 1.46 / 1.47 | 129,639,385 | 5,295,032 |
| chr22 PAM-less K=8 | 2.9.1 | 108.4 / 94.9 | 277.4 / 268.0 | 18.51 | 1,864,232,876 | 78,317,538 |
| | 2.10.0 | 107.0 / 136.2 | 274.5 / 300.1 | 18.88 | 1,864,235,403 | 78,317,730 |
| chr1 NGG, `-bMax 2` | 2.9.1 | 31.9 / 38.5 | 74.9 / 74.4 | 5.53 / 5.54 | 579,345,896 | 24,957,801 |
| | 2.10.0 | 28.8 / 32.2 | 71.4 / 70.9 | 5.63 | 579,348,060 | 24,957,823 |

amd64: ml008 (255 cores shared with other users, load average 250-500 during the runs;
g++ 14.2 -O3, the rc281 image's toolchain), the 2.10.0 build after the 2.9.1 one:

| index | engine | wall s | user CPU s | peak RSS GB | .bin bytes | leaves |
|-------|--------|--------|------------|-------------|------------|--------|
| chr22 NGG, `-bMax 2` | 2.9.1 | 15.4 | 28.1 | 1.44 | 129,640,072 | 5,295,027 |
| | 2.10.0 | 16.1 | 27.5 | 1.46 | 129,639,385 | 5,295,032 |
| chr22 PAM-less K=8 | 2.9.1 | 209.0 | 483.9 | 18.51 | 1,864,232,876 | 78,317,538 |
| | 2.10.0 | 202.9 | 479.0 | 18.88 | 1,864,235,403 | 78,317,730 |
| chr1 NGG, `-bMax 2` | 2.9.1 | 78.3 | 145.2 | 5.54 | 579,345,896 | 24,957,801 |
| | 2.10.0 | 71.8 | 139.5 | 5.63 | 579,348,060 | 24,957,823 |
| chr1 PAM-less K=8 | 2.9.1 | 1,022.7 | 2,890.0 | 101.01 | 10,468,103,808 | 460,955,944 |
| | 2.10.0 | 1,093.6 | 2,819.2 | 102.86 | 10,468,113,334 | 460,956,484 |

The `.bin` files are byte-identical across the two architectures (chr22 NGG md5s
compared for both versions; the other sizes agree to the byte).

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
   5' K=4 indexes with and without `--prune-pam`; searches at 0-3 DNA/RNA bulges; and
   search guides longer and shorter than the index's (section 5.1: 21-nt NNGRRT, as
   CRISPRme's SaCas9 PAM file, off the K=8 index under `--prune-pam` and off 20- and
   21-nt NNGRRT indexes; 21- and 19-nt NGG; 23-nt TTTV, as CRISPRme's Cas12a PAM file,
   off 20-nt TTTV and 5' K=4 indexes; 19-nt TTTV), with sites planted for them (the
   bases beyond the guide's far end known, 0..4 of them before a contig end or an N
   run; and guides cut by one or two bases at a contig end or an N run).
   Every indexed search, run through `crispritz.py`, must equal
   `test/boundary/oracle.py` -- a brute-force aligner over the FASTA that enumerates
   the sites of section 3 and aligns with the engine's alignment rules, sharing no
   code with it -- as a set of rows (251 searches, 12,400,832 rows); every row must
   be its contig's own bases at its Position (`test/boundary/check_rows.py`, which
   reads only the FASTA); every 0-bulge search must equal `crispritz.py search`
   without bulges, i.e. `searchBruteForce` (58), except that a search guide shorter
   than the index's must find a subset of it, the brute force then being the oracle's
   dedicated index for the shorter guide (8 of the 58); and each `--prune-pam` search
   must give the same rows as the dedicated index built for the same guide length
   and `-bMax` (90 comparisons, at `Gs` = 19, 20, 21 and 23). The 2.9.1 `buildTST` and
   `searchTST` fail 530 of those comparisons (run with the fixed `searchBruteForce` as
   the 0-bulge reference: the 2.9.1 one segfaults on the 5' genome). This branch
   before the guide-length fix (`546333f`, which counted `dna - rna` pad bases for any
   guide length) fails 102, all in the 84 searches whose guide length differs from
   the index's: 54 differ from the oracle (rows paired with bases off the contig or
   on an N, at Position -1 among others, for a longer guide; DNA-bulged rows lost at
   boundaries for a shorter one), 40 report rows that are not their contig's bases,
   and 8 0-bulge searches with a longer guide differ from the brute force. On the
   first version of the matrix (141 searches), each other mechanism, disabled on its
   own (pad limit, `FLANK_OFF` skip, brute-force bound, guide-only N rule), failed 10
   to 52; so did the oracle without its model of the engine's leaf node.
2. **Format-1 compatibility** (`test/run_format1_test.sh`): indexes built by the 2.9.1
   `buildTST` (committed fixtures, byte-reproducible on arm64 and amd64): rows away
   from index 0 byte-identical to the 2.9.1 `searchTST`'s, in order; index-0 rows equal
   to the oracle's. Seven searches, two of them with 21- and 19-nt guides on the 20-nt
   index (313 and 427 rows byte-identical to 2.9.1's away from index 0).
3. **Index format** (`test/run_index_format_test.sh`): header and trailer, `index-info`,
   `INDEX_FORMAT` against `CZI_FORMAT_VERSION` and the files `buildTST` writes, a
   directory mixing both formats, `max_bulges` from the header, the refusals (a
   header guide length other than the body's among them), the stdout note for a
   21- and a 19-nt search guide on a 20-nt index (none for a 20-nt one), and
   `index-info` against `searchTST` on 20 files damaged in one way each (every check
   of section 5 that does not depend on the search, a reserved flag bit (bit 3),
   `max_bulges` 32 and -1, a header PAM other than the file name's, a leaf count of 0
   and one of 2,000,000 the file cannot hold in either format among them) and, run as
   a user other than root, an unreadable one: both refuse each, for the same reason,
   and the reader of section 7 returns what `index-info` does for the 20 and the good
   files (38). The `index-info` of `520a69a` accepts 12 of the first 18. `searchTST`
   and `index-info` also refuse a directory holding a macOS `._` file beside the
   index's `.bin`. The boundary suite runs `index-info` on each index it builds (exit 0,
   nothing on stderr, the header's parameters those of the build).
4. **The four-allele symbol** (`test/run_four_allele_test.sh`, `test/four_allele/`):
   46 synthetic contigs (each case on both strands) whose split and multi-ALT SNP
   records put `Q` in guide regions (where it turns a reference mismatch into a
   match, up to six in one guide), in PAMs that exist only through it (NGG's GG, one
   G or its N; TTTV's T or V), beyond a guide's far end, next to N runs and on contig
   ends, a run of 30, ALTs on reference `N` (which stay `N`), a 20-kb background with
   3% `Q`, IUPAC codes and N runs, and a guide with one `N` beside an NGG PAM, through
   `crispritz.py`: `add-variants` with
   both enricher ports (identical; `Q` exactly at the four-base positions, `N` exactly
   the reference's), `index-genome` for NGG (`-bMax` 0, 2), NRG, a PAM-less K=8 flank,
   5' TTTV and a 5' PAM-less K=4 flank (a `.bin` flagged `four_allele` exactly when
   one of its leaves stores `Q`), and 23 searches, with and without bulges and
   `--prune-pam`, each equal to the oracle (Q as any base; 9,328,182 rows, 9,128,468
   through `Q`), the 11 without bulges equal to `searchBruteForce`, no row holding an
   `N`, every row its contig's bases, the 38 planted perfect sites found, and
   `--prune-pam` equal to the dedicated index; and `searchBruteForce` at `-mm 120`
   reports no window with an N in its guide, though 6 NGG windows hold exactly one
   there (the planted guide, both strands): the brute force of 2.9.1 and of `524b594`,
   which counted an `N` as 101 mismatches, reports 12 rows through them. The pre-`Q`
   engine (`d75c603`) fails it:
   its `buildTST` writes to stderr on `Q`, its `searchBruteForce` reads `Q` as a gap
   (109 rows only the index finds), and its enricher writes no `Q`.
5. **Existing tests** (`run_golden_test.sh`, `run_minimal_test.sh`,
   `run_pam_filter_test.sh`, `run_flank_n_test.sh`, `run_short_contig_test.sh`) pass;
   the golden outputs are byte-identical (their genome holds no boundary site).
   `run_short_contig_test.sh` now also requires a contig exactly one target long to
   be indexed.
6. **AddressSanitizer**: all three binaries, over every test above, on arm64 and amd64
   (the short-contig test without its `ulimit -v` cap, under which ASan cannot reserve
   its shadow memory).
7. **Platforms**: arm64 (Docker on Apple Silicon: g++ 14.2, and `make -f
   Makefile_conda` in `gcc:12`) and amd64 (ml008, g++ 14.2); every test passes on each.
8. **hg38**: see section 10 for build cost and leaf-level comparison, and:
   * the shipped format-1 `NGG_3_hg38` (all 455 contigs, 501 trees, not rebuilt),
     searched with the guides of the 11 contigs that open with a reverse NGG
     (`-mm 0 -bDNA 1 -bRNA 1`, 32 threads): 63,095 rows from both engines; the 11
     contigs' rows move from the + strand at 2-4 to the minus strand at 0 under 2.10.0
     (same rows otherwise); every other contig's rows are identical (in order for the
     132 single-tree contigs; the 6 multi-tree chromosomes interleave their trees in
     thread-completion order, the same rows).
   * 246 NGG contig-end sites on 171 alt/random/Un contigs, each searched with its own
     perfect guide: a dedicated format-1 index finds 235 with 2.9.1 (11 come out on the
     + strand) and 246 with 2.10.0; a format-1 K=8 index with `--prune-pam` finds 0; a
     format-2 K=8 index with `--prune-pam` finds 246, and its rows equal the format-2
     dedicated index's (368,586 rows). 2.10.0 on the format-1 dedicated index changes the
     rows of 12 contigs only: the 11 move to the minus strand at base 0, and 6 rows of
     the builder's index-0 artefact on chrUn_KI270539v1 (reported on the + strand with
     the PAM `GGC`) are dropped; the other 159 contigs are byte-identical.
   * chr22 and chr1, 50 guides (45 whose perfect match lies 0-2 nt from an N run on
     chr1/chr22, and TRAC, EMX1, VEGFA site 3, HEK4, FANCF), `-mm 4 -bDNA 1 -bRNA 1`:
     - the format-1 index gives byte-identical output under 2.9.1 and 2.10.0: chr22
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
     the format-1 index 81.5 / 87.4 / 80.7, 2.10.0 on it 70.3 / 82.9 / 78.5 (output
     byte-identical each round), 2.10.0 on the format-2 index 79.6 / 81.3 / 80.3.
   * guide length (section 5.1), chr22, the format-2 NGG `-bMax 2` index (20-nt), the
     same 50 guides, `-mm 4 -bDNA 1 -bRNA 1`, amd64, this branch against its state
     before the guide-length fix (`546333f`; the two `buildTST`s build byte-identical
     `.bin` files):
     - 20-nt guides: byte-identical output (5,234,321 rows);
     - 21-nt guides (an `A` added at the far end): 3,264,009 rows; `546333f` reported
       1,102 more, and those are exactly its rows that are not chr22's bases at their
       Position (`check_rows.py`), e.g. a "perfect match" at 16,320,445 (-) whose 21st
       base is an N; at `-bDNA 0 -bRNA 0` the indexed search equals `searchBruteForce`
       (22,499 rows), where `546333f` had 6 rows more;
     - 19-nt guides (the far base dropped): 6,020,396 rows, 24 more than `546333f`,
       all DNA-bulged alignments at sites next to an N run that use no pad base; none
       lost;
     - every row of the three searches is chr22's bases at its Position.

## 12. Not changed

* A 5' PAM shorter than the stored flank is still refused under `--prune-pam`.
* A guide-region N still discards the site, even for an alignment with RNA bulges
  that would not reach the N.
* An index still holds the sites of its own guide length `G` only, so a search guide
  shorter than `G` does not find a target whose far end lies within `G - Gs` nt of a
  contig end or an N run, and one of another length gets, with bulges, rows that
  differ from those of an index built for its length (section 5.1). Making an index
  serve every guide length alike would mean storing every PAM whatever its guide
  region holds, with a per-leaf count of usable bases wider than `pad_valid`'s 5
  bits: a format change not made here.
* Coordinates are still 32-bit ints (contigs up to 2^31 - 1 - (mb + K) bases).
* searchTST still tries no RNA bulge on the last base of a stored string (a leaf node
  of the tree), so an alignment that would leave the last guide base unpaired after
  consuming all but one of the `G + mb` stored bases is not reported, and whether such
  an alignment is reported depends on the index's `-bMax` (and, for a search guide of
  another length, on the index's `G`: section 5.1). This holds for interior
  sites as much as boundary ones, is part of the golden outputs, and is modelled by
  the oracle (without that model it predicts 352-1,786 extra rows per bulged search).
* `main.cpp` (the old single-threaded builder, not built by any recipe) is left as
  is; it does not compile against the current `pamGenerator.cpp` and never did since
  the Aho-Corasick code was commented out.
