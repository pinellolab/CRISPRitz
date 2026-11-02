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
  used). Cluster Position: `ws` + the same pad for PAM-left windows, the guide start
  `ws + mb` for PAM-right windows.

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
  strand PAM exactly `K-L` from the start lands there too. A PAM-right window that
  starts at 0 (a 3' + strand site whose guide starts at exactly `mb`, or a 5' +
  strand PAM at base 0) is built as the other shape, from the wrong bases, and its
  real site is lost.
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
CRISPRme).

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
| CRISPRitz <= 2.9.0 | as before | fails loudly: `numLeaves` negative, `std::length_error` |
| CRISPRitz 2.9.1 | as before | fails loudly: `'...' is not a valid CRISPRitz index (.bin) file (numLeaves=-28747197)`, exit 1 |
| CRISPRitz 2.9.2 searchTST | byte-identical to 2.9.1 except index-0 rows (B1) | full B1-B4 semantics |
| CRISPRitz 2.9.2 searchBruteForce | - (no index) | - ; B2 applies at once |

Byte-identical means the same `.targets.txt` bytes, row order included; tested on
the repo's fixtures and on indexes built by the 2.9.1 binary (section 11).
Users' existing indexes keep working unchanged. They keep missing the sites of
B2-B4 until rebuilt; B1 corrects their index-0 rows at once.

A format-2 index of a genome holds every site the format-1 index of the same
genome and parameters holds, with the same string, flank and coordinates, plus the
boundary sites. A search on it returns a superset of the rows of a search on the
format-1 index, the extra rows all within `mb` (B3) or `K - L` (B4) nt of a contig
end or an N run, or on the last PAM bases of a contig (B2), and minus the index-0
artefacts (B1).

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

To be filled in from the hg38 chr1/chr22 runs (section 11.5): build time and
`.bin` bytes for 2.9.1 vs 2.9.2, dedicated NGG and PAM-less K=8.

## 11. Tests

1. **Boundary suite** (`test/run_boundary_test.sh`, `test/boundary/`): synthetic
   contigs with planted perfect and 1-4-mismatch sites at every distance
   `d = 0..mb+1` from a contig start, a contig end and an N run, on both strands,
   for a 3' PAM (NGG, NRG) and a 5' PAM (TTTV), with the PAM also planted on the
   contig's last bases and within `K` of either end; `-bMax` 0-3; searches at
   0-2 DNA/RNA bulges; `--prune-pam` from PAM-less indexes with `K = L` and
   `K > L` (NGG, NG and NRG off K=3 and K=8). Every indexed search, driven through
   `crispritz.py`, must equal an independent brute-force aligner written for the
   test (it enumerates the windows of section 3 straight from the FASTA and aligns
   with the engine's alignment rules) as a set of rows, and every 0-bulge search
   must also equal `crispritz.py search` without bulges (`searchBruteForce`).
2. **Format-1 compatibility**: indexes built by the 2.9.1 `buildTST`, searched by the
   2.9.1 and the 2.9.2 `searchTST`: byte-identical except index-0 rows, which must be
   exactly the B1 corrections.
3. **Existing tests** (`run_golden_test.sh`, `run_minimal_test.sh`,
   `run_pam_filter_test.sh`, `run_flank_n_test.sh`, `run_short_contig_test.sh`)
   unchanged and passing; the golden outputs byte-identical.
4. **AddressSanitizer** build of all three binaries over the boundary suite.
5. **Platforms and scale**: arm64 (Docker on Apple Silicon) and amd64 (ml008), and
   hg38 chr1 and chr22: build time and size, 2.9.1 vs 2.9.2.

## 12. Not changed

* A 5' PAM shorter than the stored flank is still refused under `--prune-pam`.
* A guide-region N still discards the site, even for an alignment with RNA bulges
  that would not reach the N.
* Coordinates are still 32-bit ints (contigs up to 2^31 - 1 - (mb + K) bases).
* `main.cpp` (the old single-threaded builder, not built by any recipe) is left as
  is; it does not compile against the current `pamGenerator.cpp` and never did since
  the Aho-Corasick code was commented out.
