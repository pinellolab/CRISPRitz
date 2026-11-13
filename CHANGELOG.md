# Changelog

All notable changes to CRISPRitz are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When cutting a release, move items out of `[Unreleased]` into a new dated
version section and update the link-reference footer. See `docs/RELEASING.md`
and the `release-crispritz` skill (`.claude/skills/release-crispritz/SKILL.md`).

CRISPRitz is a compiled (C++/OpenMP) Bioconda package that CRISPRme pins to an
exact version. Releasing therefore also drives a downstream repin in CRISPRme
(see `docs/RELEASING.md` for the crispritz -> Bioconda -> CRISPRme chain).

## [Unreleased]

## [2.10.0] - unreleased

Not yet tagged: the date is set when it is. Released together with rebuilt indexes and the
CRISPRme release that pins it (CRISPRme 2.10.0): the IDX-13, IDX-4 and IDX-15 fixes below
change what an index holds, and every allele of a split multi-allelic site, the
four-allele symbol `Q`, the merged carriers, the split form of a record listing several
ALTs and the wider fake-indel window change what an enrichment holds, so they take effect
only on genomes enriched and indexes built by this version (the brute-force scan gains
IDX-15 at once, and IDX-14 applies to existing indexes).

### Added
- **`Q`, the four-allele symbol.** A position whose SNP ALTs, with its reference base when
  that is A, C, G or T, cover all four bases is written `Q` in the enriched genome (but a
  reference `N`, or another byte that is no base or IUPAC code, stays: see Changed), and
  `buildTST`, `searchTST` and `searchBruteForce` read it as any base (`1111`): it matches
  every guide base and every PAM letter. `N` stays a gap that matches nothing
  ([docs/INDEX_FORMAT.md](docs/INDEX_FORMAT.md) section 4.1,
  [docs/MULTIALLELIC_ENRICHMENT.md](docs/MULTIALLELIC_ENRICHMENT.md) section 2.2). In the
  tree `Q` takes the value its slot never uses (`1111` in a root, lo or hi slot, `0000` in
  an eq slot, where `1111` ends a string) and in a flank `1111`; an index file holding it is
  flagged in its header (bit 2, `FOUR_ALLELE`), so that no engine that cannot read it
  searches it, and `Q` changes no byte of a file that does not hold it. `buildTST` says on
  stdout how many of its targets (the sites it stores, with `Q` in the guide region, pad or
  flank) and files hold `Q`, `searchTST` how many of the files it reads are flagged, and
  `index-info` reports `four_allele` per format-2 file. On HGDP + 1000 Genomes chr22, 2,448
  positions carry three SNP ALTs.
- **Index format 2** ([docs/INDEX_FORMAT.md](docs/INDEX_FORMAT.md)). Every `.bin` built by
  this version starts with a header -- the magic bytes `CZI\xfe`, the format version, and
  the build parameters (PAM side, all-N PAM, guide length, PAM length, max bulges, the PAM)
  -- and ends with the magic, so an index says what it holds and a cut-short file is
  recognised. The per-leaf `'0'`/`'_'` marker becomes a flags byte (strand shape, flank off
  the contig, bases of the bulge pad that exist), so leaf records keep their size; the tree
  is unchanged but for `Q` (above). `searchTST` reads formats 1 and 2 and refuses an unknown
  version, a cut-short format-2 file, or a PAM width (but a narrower 3' PAM under
  `--prune-pam`, N-padded to it) or, in format 2, a PAM side that differs from the search's,
  before searching. `buildTST` checks its writes (a full disk used to leave a cut-short
  `.bin` behind an exit 0) and refuses `-bMax` above 31. The header flags are bit 0 (PAM at
  the 5' end), bit 1 (all-N build PAM) and bit 2 (holds `Q`); a reader refuses any other.
- `crispritz.py index-info <index_dir>` prints each `.bin`'s format and build parameters as
  JSON, reading only its header, the start of its body and its last bytes (the reference
  reader for CRISPRme, docs/INDEX_FORMAT.md section 7), and exits 1, with an `error` entry
  saying why, on a `.bin` that `searchTST` refuses whatever the search for what those bytes
  and its size show (an unknown or malformed header, reserved flag bits, `max_bulges` above
  31, a cut-short format-2 file, a header guide length other than the body's, a header PAM
  other than the file name's, a leaf count below 1 or one the file cannot hold), and on a
  directory whose `.bin` files are named for more than one PAM (such as a macOS `._` file
  copied along with an index, which `searchTST` refuses too). Its `writes_format`, and
  `INDEX_FORMAT` in `crispritz.py`, give the index format this CRISPRitz builds, for a
  reader that must know the format of an index it has just built that holds no `.bin` to
  read (docs/INDEX_FORMAT.md section 8).
- `test/run_boundary_test.sh`: sites whose guide's far end lies 0-4 nt from a contig end or
  an `N` run, or whose PAM lies 0-9 nt from one, both strands, 3' and 5' PAMs, dedicated and
  PAM-less indexes, with and without `--prune-pam`, search guides of the index's length and
  longer and shorter ones (21-nt SaCas9 and 23-nt Cas12a guides off 20-nt indexes, as
  CRISPRme searches them), through `crispritz.py`; every search must equal
  `test/boundary/oracle.py`, an independent brute-force aligner over the FASTA, every row
  must be its contig's own bases at its position (`test/boundary/check_rows.py`), every
  0-bulge search must equal `searchBruteForce` (a shorter guide: a subset of it), and every
  `--prune-pam` search the dedicated index built for the same guide length and `-bMax` (251
  searches, 12,400,832 rows; 58 brute-force and 90 dedicated-index comparisons; the 2.9.1
  `buildTST`/`searchTST` fail 530).
- `test/run_format1_test.sh`: indexes built by the 2.9.1 `buildTST` (committed fixtures)
  search byte-identically to 2.9.1 away from index 0, with the index's guide length and
  with 21- and 19-nt guides, and the index-0 rows equal the oracle.
- `test/run_index_format_test.sh`: the format-2 header and trailer, `index-info`, a
  directory mixing both formats, `max_bulges` taken from the header, the reader's refusals
  (cut-short file, PAM side mismatch, unknown version, a header guide length other than the
  body's) and `buildTST`'s of `-bMax` above 31, the stdout note for a search guide longer or
  shorter than the index's, and `index-info` against `searchTST` on 20 files damaged in one
  way each: both refuse every one, for the same reason, and the reader in
  docs/INDEX_FORMAT.md section 7 returns what `index-info` does.
- `test/run_four_allele_test.sh`: `Q` end to end through `crispritz.py` -- `add-variants`
  with both enricher ports on 46 synthetic contigs whose split and multi-ALT records put `Q`
  in guides, in PAMs that exist only through it, in the bulge pad, beside `N` runs and on
  contig ends, then 3' and 5', dedicated and PAM-less indexes, and 23 searches with and
  without bulges and `--prune-pam`: every search equals the independent oracle, every
  0-bulge one `searchBruteForce` (every target the brute force finds, the index finds), no
  row holds an `N`, a `.bin` is flagged exactly when it holds `Q`, and the brute force at
  `-mm 120` reports no window with an `N` in its guide, where a guide with one `N` is planted
  beside a PAM (an `N` cost 101 mismatches: 2.9.1's `searchBruteForce` reports it). The
  engine before `Q` fails it. Last, a contig holding `X`, `.` and `-` in guide regions, pads
  and PAM-less flanks, beside `Q`: `buildTST` writes nothing to stderr, every search equals
  the oracle and, without bulges, `searchBruteForce`, and exactly the planted perfect sites
  are found.
- `test/test_enricher_multiallelic.py`: 24 unit tests run against both ports (48), with
  exact expected outputs (every ALT, `Q`, the N-count invariant, merged carriers, each
  record's AF, a reference holding `Q` or `q`); 44 of them fail on 2.9.1. In CI with
  `verify_enricher.sh`, which grows from 17 to 34 cases (split multi-allelic sites in every
  order, genotype shapes, SNP + indel at one position, four bases, two ALTs on one haplotype
  slot, AF by key, one allele on two records, ALTs on a reference `N`, `X` or `-`, REF and
  ALTs in lower or mixed case) and now also checks the outputs against what they must be:
  the enriched genome against an independent reconstruction (with the N-count invariant),
  the dict where one allele's records merge, AFs are read by key, records sit at a reference
  `N`, `X` or `-`, a record lists three SNP ALTs or records are in lower case (cases 24 and
  31-34), the indel log's AF read by key (case 31) and its rows against the VCF's indel
  alleles in upper case, and a reference holding `q` refused.
- `docs/MULTIALLELIC_ENRICHMENT.md`: the behaviour, its measurement, what CRISPRme must change
  to use it, the owner's decisions, and the published indexes and stores to rebuild.
- **Region hooks for CRISPRme's haplotype scan: `crispritz.py search --skip-bed <bed>`,
  `--keep-bed <bed>` and `--leaf-dump <file>`** (docs/INDEX_FORMAT.md section 13; `searchTST`
  takes them as its 13th-15th arguments or as `CRISPRITZ_SKIP_BED`, `CRISPRITZ_KEEP_BED`,
  `CRISPRITZ_LEAF_DUMP`). With a skip BED (contig, start, end) an indexed search skips every
  site whose core -- its search PAM and the search-guide bases next to it, `Gs + L` bases --
  meets an interval: it returns exactly the lines it returns without the BED whose core meets
  none, read from the line as `[Position + d - r, + n)` (a PAM-right line) or
  `[Position, + n)` (PAM-left), and with the BED as a keep BED exactly the others, so the two
  partition it, each in its order. The core is the search's, not the stored window's: on the
  K=8 index searched for NGG under `--prune-pam` it leaves out the 5 distal flank bases, and
  for a guide longer than the index's it takes the pad bases the longer guide reads (the
  prototype's window rule lost lines there, and wrote others twice). The dropped leaves are
  flagged as their tree is loaded and the search enters no subtree without another leaf, so
  skipping the dense windows costs a fraction of the search: on chr22's IGL locus with HGDP +
  1000 Genomes, one guide at 4 mismatches and 2+2 bulges, the INDELS index's search,
  204,453,499 lines in 1,440 s and 157 GB without the BED of its dense windows (the 25/27
  fake-indel windows the hooks were first measured on; on 2.10.0's 30/30 windows it had
  written 51.4 GB of lines, holding 166 GB of memory, when it was stopped), gives with it
  the 728 lines outside them in 3.3 s and 0.05 GB (792 in 3.4 s and 0.06 GB on 30/30), and
  the variant index's, 33,858,344 lines in 241 s and 26.6 GB, gives 2,110 in 16.7 s and
  0.17 GB. `--leaf-dump` writes, for every
  leaf a skip drops, the letters the index stores for its window on the + strand (`-` where it
  stores a placeholder: a bulge-pad base past `pad_valid`, a flank base off the contig), which
  a reader needs to rewrite a line of a skipped site from an index downloaded without its
  genomes: on IGL they equal the enriched and fake-indel FASTAs at every letter. Every search
  given a BED prints `SKIPPED LEAVES <n> KEPT LEAVES <n>` (and `LEAF DUMP <n>`), and
  `crispritz.py search` stops when its `searchTST` prints none (one without the hooks, which
  ignores the arguments); `crispritz.py version` lists
  `search hooks: skip-bed keep-bed leaf-dump`, `crispritz.py` states `SEARCH_HOOKS`, and
  `index-info` lists the index's `contigs`, the only ones a BED may name. A malformed BED, a
  contig the index has no tree for, a skip BED with a keep BED, a dump without a skip BED, or
  one of the flags on a brute-force search is refused before the search; `crispritz.py` passes
  the three variables to no search. With no BED the search is byte-identical to the one
  before: every `searchTST` call of the boundary, golden, four-allele, PAM-filter, flank-N,
  format-1, index-format and minimal tests, run beside `4579026`'s and compared byte for byte
  (373 calls and 21,739,834 target lines on amd64, 143 and 14,525,735 under ASan, the same on
  arm64), and every IGL and IGH search without a BED of docs/INDEX_FORMAT.md section 13.6;
  merged into 2.10.0, the same 373 and 143 calls beside c18c575's `searchTST`, identical, and
  every IGL search of section 13.6 on 30/30 windows (the INDELS 4+2+2 one aside, too large to
  keep).
  CRISPRme 2.10.1 uses the hooks (CRISPRme 2.10.0 does not); section 13.5 says what such a
  reader must hold to, among it the fake-indel window's width a scan of the dense fake-indel
  windows needs. With a search guide longer than the index's, `LEFT, RIGHT >= Gs + L - 1`,
  or `- 2` towards a dense window whose interval leaves out the separator `N` between the
  two (on a PAM-less index searched without `--prune-pam`, whose PAM can run across the
  separator, `max(Gs, G + mb) + L` in place of `Gs + L`): 2.10.0's 30/30 windows hold it up
  to 31 bases, for every site CRISPRme searches; 2.9.1's 25/27 fail from 27 bases, a 21-nt
  SaCas9 or a 23-nt Cas12a guide on a 20-nt index, where the dense windows' intervals hold
  the separator after them, as the haplotype-scan design writes them (from 28 where they
  leave the separators out).
- `test/run_region_hooks_test.sh` (`test/region_hooks/`): synthetic genomes whose BED puts a
  base on, and just beside, each planted site's core (its pad, its guide region's and PAM's
  edges, the distal flank of an 8-wide window) -- dedicated 3' and 5' indexes at 0 and 2+2
  bulges, PAM-less K=8 and 5' K=4 indexes under `--prune-pam`, guides longer (21-, 23-nt) and
  shorter (19-nt) than the index's, fake-indel genomes of 25/27 and 30/30 windows, their dense
  windows' intervals with and without the separator after them -- and the format-1 fixtures
  built by 2.9.1: 24 searches, 2,280,504 lines, the skip, the keep and the skip without pruning
  checked against the search without a BED line by line, by the core computed from each line
  (`test/region_hooks/check_hooks.py`), among them the lines on which the prototype's window
  rule and the core disagree (149,662 it would drop on the K=8 index, 4,267 it would keep on
  the 21-nt route); the leaf dump letter by letter against the FASTA, and complete; 4 threads
  as 1; the rule of docs/INDEX_FORMAT.md section 13.5 for a scan of the dense fake-indel
  windows alone (`test/region_hooks/hr27_bound.py`), on planted fake-indel genomes with the
  dense windows' intervals in four forms (the separators out, the one after each, the one
  before, both): site by site on the pad route (25/27 and 30/30 windows, 21- to 26-nt guides on
  a 20-nt index, 48 skips), and where a PAM runs across a separator (a PAM-less K=8 index
  searched without `--prune-pam`, windows of 28 to 30, 12 skips); 26 refusals; an empty BED,
  and one with comments, CRLF, extra columns and overlaps; the environment; a `searchTST`
  without the hooks (`crispritz.py` stops); the capability. With `SEARCHTST_REF` every search
  without a BED is compared byte for byte with another engine's; `QUICK=1` (CI's ASan job)
  searches one mismatch fewer.

### Changed
- **New indexes are format 2, which CRISPRitz 2.9.1 and earlier cannot search** (2.9.1
  refuses one: `not a valid CRISPRitz index (.bin) file (numLeaves=-28747197)`). Indexes
  built by earlier versions keep working (see Fixed, IDX-14) but lack, until rebuilt, the
  sites B2-B4 add and the sites behind the index-0 leaves IDX-14 drops.
- **A search guide of another length than the index was built for** (e.g. CRISPRme's 21-nt
  SaCas9 or 23-nt Cas12a guides off a 20-nt index) is still searched, now with a `NOTE:` on
  stdout naming both lengths and what that limits; its rows are exact for what the index
  holds, and on a format-2 index a longer guide is never paired with a base off the contig
  or on an `N`, but with bulges they can differ from the rows of an index built for that
  length, and a shorter guide misses targets at contig ends and `N` runs
  ([docs/INDEX_FORMAT.md](docs/INDEX_FORMAT.md) section 5.1: build the index with the guide
  length you search with). The guide-length offset is now kept per index file and search
  thread (it was one global that every thread wrote), and a format-2 header whose guide
  length differs from its body's is refused as malformed.
- **Indexes built with 2.9.1 or earlier from a normalised panel hold, at a split site, only
  its last SNP record's ALT; rebuild them**, with their dicts and Tier-0/Tier-1 stores, from
  one 2.10.0 enrichment (`docs/MULTIALLELIC_ENRICHMENT.md` §8 lists the published ones).
  Cost of the fix: one byte per base; on chr22, 986 s → 1,004 s and 6.72 → 6.86 GB (2^30
  bytes; measured with the four-base rule `Q` replaced, which kept ~17 bytes more per SNP
  ALT; 2.10.0 peaks at 6.83 GB).
- **The fake-indel window is 30 reference bases on each side of an indel (25 before and 27
  after up to 2.9.1).** `add-variants` writes each indel allele's contig as 30 reference
  bases, the record's REF (from POS on) with the ALT in its place, then 30 reference bases
  (fewer where the chromosome ends), in both enricher ports (`INDEL_LEFT`, `INDEL_RIGHT`;
  `crispritz.py` states them, for a reader that must know the window of an enrichment not
  made yet). A target of L bases (guide and PAM) with b DNA bulges holds a record's change on
  its window only where L + b <= LEFT + 2 and L + b <= RIGHT + 1 (LEFT + 1 for a record
  whose REF and ALT differ at their first base, an MNV or complex record): up to 2.9.1
  L + b <= 27 (26 at an MNV), which a 27-nt site (SaCas9's 21-nt guide and `NNGRRT`,
  Cas12a's 23-nt guide and `TTTV`) with a DNA bulge, and a 28-nt one (CEP290's 22-nt
  SaCas9 guide) at any budget, exceeded, whatever the index; now L + b <= 31: LEFT = RIGHT =
  Lmax - 1 + bmax, with Lmax = 29 (a 21-nt guide and a PAM of 8 nt, the longest the general
  K=8 index serves) and bmax = 2 (the most DNA bulges CRISPRme's web offers). The log's `CHR` field records each window's
  start and end as before, so a reader that takes the left flank as `POS - 1 - start`
  (CRISPRme 2.10.0) reads the windows of either release. Measured on hg38 chr22 with
  `add-variants` (compiled port; against the 881f761 build, whose 25/27 log and fake-indel
  FASTA are 2.9.1's byte for byte): on HGDP + 1000 Genomes the same 240,557 log rows in the
  same order with the same carriers, rsID and AF, each `CHR` start 5 less and end 3 more,
  each contig the 25/27 one with the 5 reference bases before it and the 3 after it (0
  exceptions, none cut by the chromosome's end), `refseq` the reference from the new start,
  one `N` between contigs; the fake-indel FASTA 14,174,243 -> 16,098,699 bytes (+13.6%, 8
  bases a window); the enriched genome, the dict and `change_version.txt` byte-identical;
  `enricher.py` gives the compiled port's five files and run log. On HPRC (195,813 rows) the
  same, +1.6% (its windows hold long ALTs). An indel at POS 30 or
  less is now skipped (up to POS 25 before), with the stdout `WARNING` naming 31 bases.
  Every `_INDELS` index is built from these windows, so it is rebuilt with the rest, and its
  log with it (`docs/MULTIALLELIC_ENRICHMENT.md` §2.4, §8). `verify_enricher.sh` checks every
  window of its fixture against the reference (its `CHR` field, contig, `FAKEPOS`, `refseq`
  and the skips at the chromosome's start), and `test/test_enricher_multiallelic.py` both
  ports' windows and `crispritz.py`'s `INDEL_LEFT`/`INDEL_RIGHT`; both fail on 2.9.1's.
- **One allele on two records (the same `REF>ALT` twice) is one dict entry with the
  carriers of both**, where the later record's entry replaced the earlier one: a sample
  listed by both, phased alike, carries the ALT on every slot where either record puts it
  (`0|1` and `1|0` give `1|1`), otherwise the genotype with more copies of it; the rsID and
  AF are the more frequent record's, and a merged entry's genotypes use `1` for its ALT (a
  record listing several ALTs is renumbered so). 362 such alleles genome-wide in the
  published panels; on HPRC chr22 the merge puts two ALTs on one haplotype slot for 37
  sample-positions, which CRISPRme must take as PUTATIVE (owner decision B4; no CRISPRme
  release does so yet, see below).
- **A reference `N` stays `N`**: the enriched genome holds an `N` exactly where the
  reference does. A record listing two or three SNP ALTs at a reference `N` (`N>A,C`) was
  coded over it (`M`), making a gap searchable, and merging split records (see Fixed) would
  have coded them over it too. A reference byte that is no base or IUPAC code (such as the
  `X` RepeatMasker `-x` masks with), which `index-genome` and `search` read as `N` (see
  Fixed), is no longer coded over either (a record listing two or more SNP ALTs was, `N` for
  all four). The dict still lists such records. (As before, the enricher reads the reference
  upper-cased, a soft-masked base being a base, and writes each chromosome it enriches in
  upper case; `add-variants` copies a chromosome without a VCF as it is.)
- `add-variants` refuses a reference FASTA holding `Q`, in either case (in a chromosome it
  enriches, in both enricher ports, and in one it copies), with
  `ERROR: <chr> holds 'Q' at <pos>: ...` and exit status 1: the search would read it as a
  four-allele position.
- **A record listing several ALTs writes each ALT's dict entry in split form** (`1` its own
  ALT, `0` the reference or another ALT; a missing allele and the separators as they were):
  a `1|2` carrier of `C>A,T` reads `1|0` under A and `0|1` under T, where its raw `1|2` stood
  under both. Every entry of a position now numbers its own ALT `1`, as a split record's and a
  merged one's already did, so a reader no longer guesses an entry's allele number from the
  position's other entries: beside a split record of the position (`C>G`, `1|0`) CRISPRme's
  readers tied A of `C>A,CT` to allele `2` and read a `1|2` carrier's A on the insertion's
  slot, losing its target of A (the CRISPRme review's world X, which CRISPRme 2.10.0 read
  without phase on the raw numbers), and a carrier of two ALTs of one record with no other
  carrier to tell their numbers is no longer read without phase. (This revisits owner decision B3,
  which left such records' genotypes raw for CRISPRme to read; approved 2026-11-12.) A
  position whose records each list one ALT is byte-identical: on hg38 chr22 built by
  CRISPRme 2.10.0's `build-index-only` (NRG, 1 + 1, `--samplesID`, 16 threads on both sides)
  with this CRISPRitz and with c18c575, HGDP + 1000 Genomes (24 files) and HPRC (29) are
  byte-identical, the enriched and fake-indel genomes, the dict, the indel log, CRISPRme's
  stores and every index `.bin` among them. A record listing several ALTs changes only where
  an entry lists a carrier holding an allele number above 1: on hg38 chr22 with the raw HPRC
  vcfwave VCF, which keeps a site's ALTs on one record (623,085 records), `add-variants`
  writes 2.9.1's dict (c18c575 writes it too) but at 7,088 of its 529,585 keys, every one
  such a key (of the 8,616 holding a record of several ALTs), and c18c575's enriched genome,
  indel log and fake-indel FASTA byte for byte. Through CRISPRme's real CLI on the review's
  synthetic worlds (per sample, 4 + 1 + 1), an independent haplotype oracle finds every
  carrier's target, phase and variants right where records list several ALTs (two such
  records at one position, worlds WMM and WMM2; beside a split record, worlds X and Y),
  where with c18c575's dicts 8 or 9 rows a world listed a carrier for a target its
  haplotypes do not hold, and WMM2 lost a CFD-1.0 target that 2.8.1 reported; the other 20
  world runs (18 without such a record, and two whose record's genotypes tell its numbers
  either way) give c18c575's rows byte for byte. Both ports; `test/test_enricher_multiallelic.py`
  and `verify_enricher.sh` (both fail on c18c575's ports); the unit tests also run through
  both ports the split form's examples the README and `docs/MULTIALLELIC_ENRICHMENT.md` §2.3
  give.
- **The dict now shows the phase of split records as it was phased, conflicts included.**
  Each entry keeps its own record's genotypes as they were phased (verbatim for a split
  record, in split form for a record listing several ALTs, and for a merged entry the union
  of its records', both above), and records phased one by one can put two ALTs of a
  position on one haplotype slot: 672 sample-positions at 584 positions on HGDP + 1000
  Genomes chr22 (e.g. 10,661,629 NA18869: `0|1` under G and `0|1` under T). On HPRC chr22
  the records' own phase puts none there, and the merged carriers of one allele's two
  records (above) put 37, all at 38,152,965. Up to 2.9.1 the dict held one record's ALTs per
  position, so no reader met one. CRISPRme's per-sample analysis builds a CONFIRMED
  haplotype holding both ALTs from it, measured through `complete-search --per-sample` of
  CRISPRme `a5379e0` on the multi-allelic branch's enrichment of chr22:1-10,700,000,
  indexed by CRISPRitz 2.9.1 in format 1, as no CRISPRme release reads a 2.10.0 index yet;
  2.10.0 enriches those inputs into the same dict and the same bases but one, a `Q` 13 kb
  from the measured site. CRISPRme must handle these before it uses a 2.10.0 index
  (`docs/MULTIALLELIC_ENRICHMENT.md` §6.3).

### Fixed
- **A site whose guide's far end (away from the PAM) lay within `max_bulges` nt of a contig
  end or an `N` run was not indexed, even a perfect match (IDX-13).** The whole
  pad + guide + PAM window had to fit with no `N` outside the PAM, though a 0-bulge
  alignment needs none of the pad, so an index with a larger `-bMax` lost more (a guide's
  own on-target could vanish from the report). A site is now kept when its guide and PAM
  fit, its guide holding no `N`; pad bases off the contig or on an `N` are recorded as
  unavailable, and an alignment that would need them is not reported.
- **`--prune-pam` off a PAM-less index whose flank `K` is wider than the PAM `L` lost PAMs
  within `K - L` nt of a contig end (IDX-4)**, because the distal flank bases the search never
  tests had to fit, and a window starting before base 0 could not be stored (the sign of a
  leaf index was its strand). A PAM-less window may now run off the contig on the flank side
  (those bases stored as `N`), and format 2 stores a signed window start with the strand in
  its own bit. With this, `--prune-pam` equals a dedicated index built for the same guide
  length and `-bMax` at contig ends too (not one built for another guide length: section
  5.1 of docs/INDEX_FORMAT.md).
- **A PAM on a contig's last bases was never found (IDX-15)**, in any index or by the
  brute-force scan: the scan stopped one short on a 3' PAM, and on a 5' PAM bounded on the
  target length, losing every reverse-strand 5' PAM within a guide length of the end. Fixed in
  `buildTST` and `searchBruteForce`; the brute-force scan gains the sites at once.
- **A reverse-strand PAM at a contig's first base was reported on the `+` strand (IDX-14)**:
  the builder wrote leaf index 0 as minus strand for a 3' PAM, the search read it as plus.
  The search now reads it as written, so existing indexes report it correctly with no
  rebuild; a leaf the 2.9.1 builder stored at index 0 from the wrong bases (a `+` window
  starting at 0), which could surface as a match without a PAM, is dropped, and the site it
  stood for (on the `+` strand, a 5' PAM at base 0 or a 3' site whose guide starts at base
  `max_bulges`) is found once the index is rebuilt. Every other row on a format-1 index is
  byte-identical to 2.9.1's.
- `searchBruteForce` segfaulted on a header-only FASTA (exit 139): its unsigned scan bound
  wrapped. It now returns at once on a contig shorter than one target.
- `searchBruteForce` wrote and read past its per-guide profile table for a PAM longer than
  3 nt (TTTV, NNGRRT); a `-t` run printed garbage in the profile's mismatch columns (43,739
  three-mismatch targets where the targets file has 48). Targets were not affected.
- **A multi-allelic site split into one record per ALT kept only its last ALT** (the 2.9.1
  known issue). `bcftools norm -m -any`, which builds every CRISPRme panel, writes a site
  with several ALTs as one record per ALT; each record rebuilt the position's IUPAC code
  from the reference base and overwrote its SNP dict entry, so only the last record in the
  file reached the enriched genome (hence the index) and the dict (hence CRISPRme's Tier-0
  and Tier-1 stores). Now each record adds its single-base ALTs (A, C, G or T) to what
  earlier records coded — the code is that of the site written as one record, in any record
  order — and the dict value lists every ALT of one character, one `$`-joined entry each
  with its own record's carriers, genotypes, rsID and AF (the form a record listing several
  ALTs always had; a repeated allele merges into its own entry, see Changed). Both ports,
  byte-identical. Measured with 2.10.0 on chr22 of HGDP + 1000 Genomes through
  `add-variants`: 43,389 positions carry two or three SNP records; 45,837 SNP records
  (2.6 %, 5,104 with AC ≥ 41) were not the last and reached neither output; after the fix
  exactly those 43,389 bases and 43,389 dict values change, each equal to an independent
  bcftools reconstruction, and the indel log and fake-indel FASTA are byte-identical;
  `enricher.py` gives the same outputs md5 for md5. The other published panels are affected
  too (checked the same way on chr22 with the multi-allelic branch, before `Q` and the
  merge): HPRC (normalised as CRISPRme builds it) 4,765 multi-record positions, mega 16,103.
  This fix changes nothing where no position holds more than one SNP record, but other
  changes of this release do (docs/MULTIALLELIC_ENRICHMENT.md section 2.6): a record listing
  three SNP ALTs over an A, C, G or T base, or four over a reference IUPAC code, is now `Q`,
  one listing two or three SNP ALTs at a reference `N` (or two or more at another byte that
  is no base or IUPAC code) leaves it, one listing `AF` at another place in its INFO than
  the file's first record gets its AF, one listing several ALTs writes its entries in split
  form (see Changed), a SNP outside the chromosome is an error, and a reference holding `Q`
  is refused; and the wider fake-indel window (see Changed) moves every indel's log row and
  contig. So an enrichment of a VCF that keeps a site's ALTs on one record can change too:
  on chr22 with the raw HPRC vcfwave VCF (no position with two SNP records), against 2.9.1,
  the genome differs at 61 bases, each `N` → `Q`, the dict at 7,088 of its 529,585 keys
  (the split form), and the indel log and fake-indel FASTA by the window alone (202,404
  rows, 0 exceptions); before the split form and the window, the dict, indel log and
  fake-indel FASTA were byte-identical.
- **A record listing three SNP ALTs over an A, C, G or T base, or four over a reference
  IUPAC code, wrote `N`**, which `index-genome` and `search` read as an assembly gap: every
  target through the position, the reference one included, left the index (shown through
  the CLI: 0 targets found, against 1 with the fix). Merging split records would have made
  it common (2,448 chr22 positions). A position whose ALTs, with its reference base when
  that is A, C, G or T, cover all four bases, at a reference base or IUPAC code, is now `Q`
  (see Added): no ALT is left out, the dict lists them all, and one stdout line counts such
  positions.
- **`add-variants` read each record's `INFO/AF` at the place it had in the file's first
  record**, so a record whose INFO held `AF` at another place (every record, if the first
  had no `AF`) got no AF, in the dict and the indel log: the retired 1000G2021_HGDP index (a
  merged VCF) lost the AF of 224,981 of 232,384 chr22 indels, and HPRC chrY 91,950 of its
  102,665 indel AFs. Each record's AF is now found by its key (exactly `AF`), in both ports.
- **`add-variants` read a record's REF and ALTs in the VCF's case**, where VCF 4.x reads their
  bases whatever their case. A SNP ALT written in lower case (`C>g`) was no A, C, G or T to the
  enricher, so the enriched genome kept the reference base there and no search read the
  allele, while the dictionary listed it and its carriers (CRISPRme's synthetic world WLOW:
  two carriers' perfect target never reported); a record that differed from another only by
  case (`C>G` and `c>g`) was another allele, two dictionary entries where one `REF>ALT` on two
  records is one (owner decision B2); and an indel or MNV in lower or mixed case got its
  fake-indel window with the ALT as written and its log row named it so (the `_INDELS` index
  was the same: `buildTST` reads a FASTA upper-cased). Both ports now upper-case REF and every
  ALT as they read a record (ASCII `a`-`z`; an ALT holding `<`, a symbolic allele, as written;
  `*` and `.` hold no letter): such a SNP ALT is coded, such records merge as two records
  written alike do, and the window, the log row and the dictionary entry are those of the
  record written in upper case. rsID, AF and genotypes are read as before, and a VCF without a
  lower-case allele gives the same outputs, byte for byte. Lower case comes in through
  normalisation (`bcftools norm -f` writes a soft-masked reference's case into the indels it
  left-aligns): HPRC normalised as CRISPRme builds it holds 108,178 such records over its 24
  chromosomes (1,546 on chr22), every one an indel or MNV, where HF's raw vcfwave VCFs hold
  none; the mega's merged VCFs hold 13, on chr3, chr4, chr5, chr10, chr12, chr18 and chr21
  (deletions of 238 to 317 bases whose REF is in mixed case and whose ALT is one lower-case
  base, none an SNV; docs/MULTIALLELIC_ENRICHMENT.md section 2.9 lists them); and the
  published HGDP + 1000 Genomes VCFs, normalised the same way, hold 6 of their 144,858,810
  records, on chr1, chr2, chr10 (two), chr17 and chr21 and none on chr22 (deletions of 228 to
  282 bases of FILTER `.`, their REF in mixed case and their ALT one lower-case base, each
  carried by 1 to 5 of the 4,091 samples; section 2.9 lists them too). Measured through
  `crisprme.py build-index-only` of CRISPRme `6bff65b` (16 threads) against `26ffe46`: HGDP +
  1000 Genomes chr22, HF's raw HPRC chr22 and the mega chr22 byte-identical, every file
  (enriched genome, dictionary, indel log, fake-indel genome, every `.bin`); HPRC chr22 as
  CRISPRme normalises it and the mega chr21 differ only where the VCF writes lower case: the
  log's `indel` field at the 1,546 (2) rows naming such a record and the fake-indel genome's
  1,818 (2) lower-case bases, upper-cased, every `.bin` byte-identical, as do the mega's
  whole recipe on each of its seven chromosomes, at its 13 rows and 13 bases, and HGDP + 1000
  Genomes on each of its five such chromosomes (`add-variants` alone; chr21 through
  `build-index-only` too), at its 6 rows and 6 bases (docs/MULTIALLELIC_ENRICHMENT.md section
  2.9). CRISPRme readers that read the VCF itself must read it in upper case too, while those
  that read the dictionary stay as they are (section 6.7): with this engine CRISPRme 6bff65b's
  mega recipe fails its check on each of chr3, chr4, chr5, chr10, chr12, chr18 and chr21
  (measured on all seven) until three of its VCF readers do.
  Tests: five unit tests, run against both ports, and `verify_enricher.sh` case 34; each
  fails on `26ffe46`.
- **`searchBruteForce` reported a window with an `N` in its guide region when the mismatch
  budget exceeded 100**: an `N` counted as 101 mismatches. Such a window is now dropped
  whatever the budget, as an index never holds one.
- **`searchTST` read a stray `.bin` as part of the index when the directory listed it after
  the index's own files**: a macOS `._` (AppleDouble) file copied along with an index was
  then read as a tree of 118,883,584 leaves, failing after 22.5 s (listed first, it was
  refused). A directory whose `.bin` files are named for more than one PAM is now refused
  before any file is read, whatever the order.
- **A `.bin` whose leaf count the file cannot hold** (damaged; a format-1 file has no end
  marker) made `searchTST` allocate that many leaves and fail in `std::bad_alloc`; it is now
  refused before the search, as by `index-info`.
- **`add-variants` took every file of the genome folder for a chromosome**: a README or
  `.DS_Store` beside the FASTA files was copied as one, failing with `cp: cannot stat` on
  stderr (which fails a CRISPRme search), a macOS `._<name>.fa` (the AppleDouble file macOS
  tar writes for a file's extended attributes) was copied into the enriched genome, and a
  README, `.DS_Store` or `.fai` listed first gave its extension to every chromosome without
  a VCF, none of which was then copied (exit 0). Only the `.fa` and `.fasta` files not
  hidden are chromosomes now, and a folder holding none is refused.
- **A genome byte that is no base, IUPAC code, `N` or `Q`**, in either case (such as the `X`
  RepeatMasker `-x` masks with), **broke the index around it.** The PAM scans of `buildTST`
  and `searchBruteForce` skipped it, so every later PAM position of the contig was read one
  base off and a site after it lost, and `buildTST` stored it in the guide region, pad or
  flank of every site whose window held it, with a message on stderr (which fails a CRISPRme
  stage), as the nibble `0`, which the search reads as an empty slot of the tree: the sites
  after it were lost, and rows the contig does not hold were reported. Such a byte is now
  read as `N`, a gap, by both: a guide region holding one is no target, a pad ends at it and
  a flank stores it as `0000`; `buildTST` says on stdout how many bytes of a chromosome it
  read so. On 20 random 400-bp contigs with 1 % of `X`, `-` and `.` and 80 planted guide
  sites, the 2.9.1 `buildTST` wrote 297 lines to stderr and its index gave 17 rows of a
  3-mismatch search without bulges; it now writes none and gives 58, the rows of
  `searchBruteForce`.
- **A SNP outside the chromosome** (POS 0 or past the end, e.g. a VCF of another assembly)
  now fails with an error naming the position, in both ports: `enricher.py` raised an
  `IndexError` past the end and coded the last base for POS 0, `enricher.cpp` read and wrote
  out of bounds.
- **`searchTST` called directly with more bulges than an index's `max_bulges` crashed**
  (`crispritz.py search` lowers `-bDNA` and `-bRNA` to the index directory's `max_bulges`
  and never makes the call; CRISPRme neither). It now refuses, exit 1 with one line on stderr
  naming the bulges asked for and the `max_bulges` they exceed, before it opens any output
  file: more DNA or RNA bulges than a format-2 header's `max_bulges` (`'<file>' was built
  with max_bulges 2 and holds no alignment with more bulges of either kind, but the search
  asks for 0 DNA and 3 RNA bulges ...`), or than its own `max_bulges` argument (the 10th),
  which sizes its per-bulge tables and is the one a format-1 index is read with. On the
  golden genome's NGG and TTTV indexes built with `-bMax` 0, 1 and 2, searched with every
  `max_bulges` argument 0-3 and every bDNA and bRNA 0-3 (384 calls), 174 of the 304 calls
  asking for more ended in a segfault once the targets file was open, every one with more
  RNA bulges than the argument among them, and the other 130 exited 0 without the
  alignments the index cannot hold (a `-bMax 0` index searched with 3 DNA bulges: 4 rows,
  where a `-bMax 3` index gives 82, 78 of them DNA-bulged); the 2.9.1-built format-1 fixture
  segfaulted in 4 of its 5 such calls. Now those 309 calls are refused and the 93 others
  give the exit status, targets and stdout they gave. Pre-existing (4579026, and 2.9.1 for
  the argument); `test/run_index_format_test.sh` makes four such calls, which fail with
  c18c575's `searchTST` (docs/INDEX_FORMAT.md section 5).

## [2.9.1] - 2026-10-30

### Added
- **`search --prune-pam`: one generous-flank PAM-less index serves ANY PAM up to
  the stored flank width (N-padding).** Build a single PAM-less index with a wide
  all-N flank of width `K` (the declared PAM length; recommended generous default
  `K = 8`), then serve *any* PAM of length `L ≤ K` on the fly — no dedicated
  per-PAM index. At TST leaf emit (`saveIndices` + `saveIndices2`), each
  candidate's stored concrete flank is bit-AND-compared against the IUPAC-expanded
  search PAM and skipped **before it is written or scored** if incompatible. A 3′
  search PAM shorter than `K` is **N-padded to width `K` at its distal end**, because
  a real PAM is always adjacent to the protospacer (proximal positions carry the real
  letters; the pad N's cover the outer stored positions the short PAM does not
  constrain). This yields the **same off-target loci** as a dedicated per-PAM index
  (proven in `test/run_pam_filter_test.sh` for NGG, NRG, NRR, NNNN and a 6-nt PAM off
  a `K = 8` index, and for an equal-width 5′ TTTV: an 8-flank NGG search == a
  dedicated 3-flank NGG index as a set). A 5′ PAM must equal `K`: a shorter one is
  refused, because the emit path trims the pad from the end of the reported sequence
  while a 5′ pad sits at its front. `L > K` is a **hard error** (the index never
  stored those bases). Enable via `crispritz.py search … --prune-pam`, the
  12th positional arg to `searchTST`, or the `CRISPRITZ_PAM_PRUNE` env var (any value but
  empty or `0`).
  **Default (flag absent) is a strict no-op — byte-identical to upstream.**
  Workflow + build recipe in README §3.3.3.
- **Index-size cost of a wider flank is linear, not multiplicative.** The flank is
  4-bit-packed (2 bases/byte) per leaf, so widening the flank adds a *constant*
  per-leaf cost. Measured on a 5 Mb hg38 chr22 reference slice (~10 M leaves,
  mismatch-only index): `K = 3` → 216.0 MB, `K = 8` → 236.0 MB, i.e. **+2 bytes/
  leaf, +9.3%**. Extrapolated to a full hg38 reference index (~6.2 × 10⁹ leaves):
  ≈134 GB (`K = 3`) → ≈146 GB (`K = 8`), **+≈12 GB (+9.3%)**; a variant/enriched
  index has more leaves but the same +2 bytes/leaf, so the percentage is unchanged.
  A generous `K = 8` flank is a ~10% premium that removes the need to store one
  index per PAM.
- `test/run_pam_filter_test.sh`: 31-assertion regression covering default-off
  byte-identity, the dedicated-index equivalence (3′, and 5′ at equal width), the
  generous-flank N-padding (NGG/NRG/NRR/NNNN(4)/6-nt PAM off a K = 8 index with exact
  PAM partitioning), the 8-flank==3-flank NGG equivalence, the `L > K` hard error, and
  the refusals of a shorter 5′ PAM off a wider flank and of a width mismatch with the
  filter off.
- `verify_enricher.sh` grows from 8 to 17 synthetic cases: the SNP+indel combinations
  (SNPs inside an indel's window on the same haplotype, the other haplotype and another
  sample; a deletion spanning SNPs; overlapping indel windows; multiallelic records
  mixing SNP and indel ALTs; unphased, missing and haploid genotypes; a REF that does
  not match the FASTA), the chromosome edges, and multiallelic SNPs whose ALTs include
  `*`. Besides the four outputs it now compares the set of files each port writes and
  their run logs.
- `test/run_add_variants_status_test.sh`: `add-variants` fails, naming the chromosome,
  when an enricher exits non-zero or is killed (see Fixed).

### Changed
- **The compiled `enricher` is now installed, so it is the default enrichment path.**
  `conda/build.sh` copies it to both places `crispritz.py`'s `_enricher_command()`
  checks — `$PREFIX/bin` (PATH) and beside the reference implementation — matching what
  CRISPRme's `install_from_source.sh` already did. Measured **1.82x** faster than
  `enricher.py` through the real `add-variants` CLI on the FULL chr22 of HGDP + 1000
  Genomes -- 1,998,332 records, 4,091 samples, 2.7 GB of output: 1,776 s vs 974 s -- on the
  stage that dominates a from-scratch genome-wide build. Byte-identity holds at that same
  scale: all four enricher outputs match `enricher.py` hash for hash (as does
  `INDELs_genome/change_version.txt`, fixed text that `crispritz.py` writes). It holds
  off Linux too: built natively on macOS arm64 (Apple clang, libc++ rather than
  libstdc++), the port reproduces those full-chr22 hashes exactly, and so does
  `enricher.py` under pandas 3
  (the recipe leaves pandas unpinned, and the indel log goes through a DataFrame). `zlib` was
  already in the recipe's build/host/run deps and `enricher.cpp` carries no OpenMP
  pragmas, so this ships on every platform the other three binaries do, through the same
  compiler invocation. Gated on the byte-identity fix above and on
  `verify_enricher.sh` covering the shape that broke it.

### Fixed
- **The C++ `enricher` dropped a SNP wherever an indel shared its position,** so it was
  not the byte-identical port of `enricher.py` that `add_variants.sh` claimed. Python's
  `iupac_code` dict holds only the 2-, 3- and 4-base combinations (12 + 24 + 24 entries,
  none of length 1), so when a record's allele set collapses to a single base the lookup
  finds no key and the enriched genome is left untouched. `iupacFromMask` instead returned
  that base, and the caller wrote it — clobbering the IUPAC code a SNP record at the same
  position had already placed. The set collapses on **every** indel record: `REF` is
  length 1 so `SNPsProcess` does not bail out, and the longer `ALT` is skipped, leaving
  only the reference base. A position carrying both a SNP and an indel as **separate
  biallelic records** — how every CRISPRme-supplied VCF is normalised — therefore lost its
  SNP from the enriched genome, which is precisely the co-occurrence the SNP+indel feature
  exists to find. Measured on a real 120,000-record chr22 slice with 4,091 samples: 176
  such positions were the **only** bytes differing from `enricher.py` in 50.8 MB
  (`R`→`G` ×39, `Y`→`T` ×38, `M`→`C` ×37, `W`→`T` ×30, `K`→`T` ×20, `S`→`G` ×12); after the
  fix the two ports are byte-identical across all four enricher outputs.
  **Exposure.** A *packaged* install never ran the compiled port: `conda/build.sh`
  compiled `enricher` and did not install it, and CRISPRme's Dockerfile built only
  `buildTST`/`searchTST`/`searchBruteForce`, so `_enricher_command()` fell back to
  `enricher.py`. Verified inside the shipped image: no `enricher` on `PATH` and none
  beside the reference. A **source** install is different — CRISPRme's
  `install_from_source.sh` has compiled and installed the binary to both places
  `_enricher_command()` checks, so such an install *did* use the C++ port, and an
  enriched genome it produced is missing the SNP at every position that also carries an
  indel. **No published index is affected:** the enriched genomes behind them were built
  through the image, confirmed by direct measurement on the artifacts — chr22:10684299
  carries a SNP (`C>A`) and an indel (`C>CTATA`) as separate records, and every enriched
  `chr22.enriched.fa` on the build host holds `M` there (the Python enricher's output),
  not `C`.
- **`verify_enricher.sh` could not have caught that** (and did not: it passed). Its
  SNP+indel case put both alleles in ONE multiallelic record, while production VCFs are
  split — 0 of 1,998,332 real chr22 records carry a comma in `ALT`. New case 8 plants a
  SNP and an indel at the same position as two separate biallelic records, which is the
  shape that fails; it reproduces the defect on the synthetic fixture
  (`GGCTAWCACATCG` vs `GGCTAACACATCG`) and passes after the fix.
- **A non-`.bin` file in the index directory crashed the search.** The directory scan
  skipped only `.` and `..`, so a sidecar such as `.pam_build` or `index_meta.json` was
  stat'd, pushed into the file list, and handed to `loadTST` as if it were a tree. The
  first four bytes became a leaf count, driving a huge allocation, and a later
  `numNodes = 0` wrote out of bounds in `deSerialize`. The scan now accepts only regular
  files whose name ends in `.bin`, and `loadTST` refuses a stream whose leaf count is not
  positive rather than trusting it. Covered by assertion 12 of
  `test/run_pam_filter_test.sh`.
- **A genomic `N` in the stored PAM/flank region dropped the target.** An `N` converts to
  a code that intersects nothing, so it failed a PAM match at every position even where
  the PAM stated no requirement — and a PAM-less index stores a wider flank than a
  dedicated index for a shorter PAM, so it was more exposed. Three gates acted on this,
  which is why fixing any one alone changed nothing: the build-time PAM scan tested every
  position, the build discarded any window containing an `N` over `pamlen + max_bulges`,
  and the prune test covered the N-pad that widens a search PAM to the stored width. The
  first is now relaxed only for a fully unconstrained PAM, the second is narrowed to the
  guide region (which cannot hold an `N`, the tree having no encoding for one), and the
  third tests only the positions the real PAM covers. A stored `N` is the `0x0` nibble,
  which the loader's `case 0x0`, unreachable in a well-formed stream until now, is free to
  carry. The writer's nibble shift was also rewritten to depend on the slot instead of the
  value; that changes no byte (shifting a zero nibble is a no-op: identical output for all
  6,725,600 PAM strings of length 1-8 over ACGTRYN), and fixed nothing. Impact is
  concentrated on the synthetic indel genome, which is `N`-separated roughly every 56
  bases: a width-31 window admitted 34,690,315 start positions on one contig against a
  width-26 window's 41,454,105, so 16.3% were invisible to a wider index. Real assembly
  `N` sits in megabase gaps where nothing matches, so reference and SNP searches were
  unaffected — a dedicated index built by the patched builder is bit-identical to one
  built before, and the patched engine returns identical sorted rows on pre-fix indexes.
  Covered by `test/run_flank_n_test.sh`, which plants targets on both strands.
- **A multiallelic site compounded its indel alleles into one contig.** `indel_to_fasta`
  assigned each allele's substitution back into the shared reference window, so allele 2
  was applied on top of allele 1: its contig carried both insertions — a haplotype no
  sample need carry — and the logged `refseq` and FAKEPOS drifted with its length. On a
  synthetic `G -> GA, GAA` the second contig inserted three bases for a two-base allele.
  Not reachable through split (biallelic) input, where the loop runs once per record, but
  `add-variants` also accepts VCFs merged into multiallelic sites. Each allele is now
  derived from the reference window; for a single-allele record the result is
  byte-identical, so no existing index moves. Also `start_position = POS - 26` had no
  floor: the window is anchored with the REF base at index 25, which a negative start
  cannot satisfy, and the two ports did not even agree — Python sliced from the far end of
  the chromosome and returned an empty window while C++ threw from `substr`. Both now skip
  the site with a warning on stdout. `verify_enricher.sh` gained a site with two carried
  indel alleles (its previous "multiallelic indel" case was one SNP allele plus one
  insertion, so the per-allele loop ran once and could not catch this) and now runs in CI.
- **A PAM/flank width mismatch mis-parsed the index, and a crashed engine read as an
  empty result.** `loadTST` consumes exactly `ceil(len(search PAM)/2)` stored PAM bytes
  per leaf, so the search PAM must have the width the index stored. Only the
  `--prune-pam` path reconciled the two; with the filter off, a width-8 PAM-less index
  searched with a 3-nt PAM desynchronised the leaf stream at the first leaf and ended in
  a SIGSEGV. A segfault writes nothing to stderr and leaves a 0-byte targets file, and
  nothing inspected the engine's exit status, so the caller saw a completed step and
  reported a clean, empty off-target result — indistinguishable from a real negative.
  The width check now runs whether or not the filter is on and refuses a mismatch with
  an error naming both widths; and `buildTST`, `searchTST` and `searchBruteForce` run
  through `_run_checked`, which exits non-zero on a bad status and names the signal.
  `--prune-pam` also refuses a shorter 5′ PAM off a wider flank: its loci were right,
  but every reported row carried the pad N's inside its PAM region, so a downstream
  PAM filter would drop them all.
- **`enricher.py` dropped a real SNP at a multiallelic site whose ALTs include `*`.** It
  built a SNP's IUPAC code from the set {FASTA base} + {1-character ALTs} and looked the
  joined string up among A/C/G/T spellings, so one member that is not A/C/G/T — the `*`
  of a spanning deletion, an `N`, a lowercase base — left the base unchanged. At a
  GATK-style record such as `C>A,*` the SNP `A` was therefore never coded, and an
  off-target needing it could not be found, although the SNP dict listed its carriers.
  The C++ port coded the A/C/G/T subset, so the two ports also disagreed on such input.
  `enricher.py` now follows the C++ rule; wherever it wrote a code before, it writes the
  same one. Split (biallelic) input, which every CRISPRme-supplied VCF is, cannot reach
  this. Pinned by `verify_enricher.sh` case 17.
- **`add-variants` exited 0 when an enricher failed.** It runs one enricher per
  chromosome in a multiprocessing pool and never read the results, so an enricher that
  crashed, was killed (by the OOM killer, say) or ran out of disk left `add-variants`
  successful, with that chromosome's enriched genome, SNP dict, indel log and contigs
  missing or truncated; CRISPRme's index build, which trusts that status, then indexed
  whatever was there. Each worker now returns its enricher's exit status, and
  `add-variants` exits 1 naming every failed chromosome with its exit code or signal.
  Covered by `test/run_add_variants_status_test.sh`, which also fails rather than
  hangs if a worker stops returning. That test then found two failures the old code
  had been hiding, both in the `enricher.py` path, and both are fixed: the script ran
  through its `#!/usr/bin/env python` line, so with no `python` on PATH (Ubuntu 24.04
  ships only `python3`) every chromosome failed with exit 127, and it now runs under
  `crispritz.py`'s own interpreter; and two workers finishing together could race to
  create the shared `SNPs_genome/<genome>_enriched/` directory, the loser raising
  `FileExistsError` before writing its chromosome (3 of 6 runs on two small
  chromosomes), which it now creates race-safely, as `enricher.cpp` already did.
- **`buildTST` ran away on a contig too short to hold one target, for an all-`N` PAM.**
  Its PAM scan's bounds are `length - k` on an unsigned length, so a sequence shorter than
  the target window wraps them to about 2^32. For a constrained PAM the scan then read
  past the sequence and segfaulted, as it always had; the genomic-`N` fix above makes an
  all-`N` PAM skip that read, so on such a contig the scan instead ran about 2^31
  iterations pushing indices, reaching 20-28 GB within 15 s until it was killed for
  memory, on a shared host possibly along with other users' jobs. A header-only FASTA is
  enough, and the enricher writes one for every chromosome with no carried indel.
  `searchPAMonGenome` now returns no target for a sequence shorter than
  `pam_len + max_bulges`, which holds none on either strand, so `buildTST` finishes at
  once, writes no `.bin` and exits 0 for any PAM, and `index-genome` skips such a contig
  instead of stopping at it now that `buildTST`'s exit status is checked. Found by the
  release review. Covered by `test/run_short_contig_test.sh` (header-only, 10-bp and
  one-base-short contigs; 3′ and 5′ all-`N` and NGG PAMs; under a timeout and a 4 GB
  memory cap), which the previous code fails at once with `std::bad_alloc`.

### Known issues
- **Under `--prune-pam` the profile outputs (`-p`/`-t`) differ from a dedicated index's,**
  though the targets list is identical as a set. The profiles count the engine's duplicate
  rows, which only the targets list removes, and how many there are depends on the index
  tree's shape; their guide label also carries the search PAM padded to the stored width
  (`…GAGNNNNNNNN` for NGG off a `K = 8` index). CRISPRme does not read the profiles.
- **At a position carrying more than one SNP record, only the last record's ALT is
  kept.** Normalising a VCF splits a multiallelic site into biallelic records, and each
  record rebuilds the position's IUPAC code from the reference base and overwrites its
  SNP dict entry, so in both ports an earlier record's ALT reaches neither the enriched
  genome nor the dict, and an off-target that needs it cannot be found. Pre-existing in
  every release. On HGDP + 1000 Genomes chr22, 43,389 positions carry more than one SNP
  record, and 45,837 of 1,757,775 SNP records (2.6%) lose their ALT this way. A fix
  changes the enriched genomes, and with them every index built from one, so it is
  being investigated separately.

## [2.9.0] - 2026-10-08

### Removed
- **CRISPRitz is now a pure enumeration engine.** All scoring, annotation and
  reporting have been removed; CRISPRitz keeps only `index-genome`, `search` and
  `add-variants` (plus `version`/help). Removed subcommands: **`scores`** (the
  vendored azimuth/Doench-2016 on-target model + CFD), **`annotate-results`**, and
  **`generate-report`**; and the **`-scores`** option of `search`. Deleted the
  vendored `Scores/` (incl. azimuth), `Annotator/`, `Plot/` and `ProcessData/`
  trees. Downstream scoring/annotation/reporting now live in **CRISPRme+**, which
  already implements them (CFD + CRISPR-Bulge scoring, interval-tree annotation,
  and its own report/plots) and never used CRISPRitz's versions.

### Changed
- **Dependencies unpinned — the numerical stack now floats free.** The exact pins
  `scikit-learn==1.1.3 / numpy==1.24.4 / pandas==2.0.3 / scipy==1.10.1` existed
  **solely** for the azimuth model. With azimuth gone, `scikit-learn`, `scipy`,
  `matplotlib`, `biopython`, `statsmodels` and `intervaltree` are dropped, and
  `numpy`/`pandas` are unpinned (the `add-variants` Python-fallback enricher uses
  pandas). This lets packages that require a modern numpy/pandas (e.g. **CRISPRme+**)
  depend on `crispritz` from Bioconda without an unsolvable conflict.
- **Recipe test** now exercises every engine subcommand (`index-genome`/`search`/
  `add-variants` help) instead of only `crispritz.py` — so a broken import or a bad
  dependency bump fails the build (the full index+search golden-output regression
  remains in the GitHub Actions CI).

### Migration
- Users who relied on `crispritz scores` / `annotate-results` / `generate-report`
  (or `search -scores`) should use **CRISPRme+**, which provides scoring (CFD +
  CRISPR-Bulge), annotation and reporting. The search/index/enrich engine is
  unchanged and byte-identical.

## [2.8.4] - 2026-10-08

### Fixed
- **`add-variants` enricher assigned samples to the wrong ALT allele on
  multiallelic records.** The enricher matched a sample to an ALT allele by
  searching the raw genotype text for the allele number and stopping at the
  first hit (`gt.find(str(allele))` + `break`), which failed two ways: (1) the
  substring test let a record with ten or more ALTs bind a sample to an allele
  it does not carry (GT `0|12` substring-matches allele `1`), and (2) the
  `break` dropped the second allele of a sample carrying two different ALTs
  (`2|1`), leaving that allele with no carrier. The genotype is now parsed into
  integer alleles (splitting on `|` and `/`, so phased, unphased, haploid and
  partially-missing calls are all handled, skipping `.` and the reference),
  membership is tested on those integers, and the sample is recorded under every
  ALT it carries. Applied to the multiallelic SNP path and the indel path in
  both the C++ enricher and the Python fallback. The biallelic path is unchanged
  (its `gt.find('1')` is correct because a biallelic GT holds only `0`, `1` or
  `.`) and the enriched FASTA is byte-identical in every test, so the search
  space does not move — only per-sample carrier assignment on multiallelic sites
  is corrected. (#44)

## [2.8.3] - 2026-09-24

### Fixed
- **`add-variants` no longer feeds tabix/bcftools index sidecars to the enricher.**
  The VCF-directory listing filtered out `.tbi` files with
  `for file in listChrs: if file.endswith(".tbi"): listChrs.remove(file)`, which
  mutates the list while iterating over it and therefore **skips elements** — so some
  `.tbi` files survived the filter and were passed to the enricher, which `gzip.open()`s
  each input and reads it as text. Decompressing a binary `.tbi` index then raised
  `UnicodeDecodeError: 'utf-8' codec can't decode byte 0x9d`. The listing is now built
  with a single comprehension that excludes both `.tbi` and `.csi` index sidecars, so
  only real VCFs reach the enricher. (Surfaced building a genotyped 1000G index whose
  VCF folder was tabix-indexed; the stray index workers were harmless — they died before
  writing — but produced confusing tracebacks and were a latent collision risk.)

### Changed
- **License changed from AGPL-3.0 (dual) to the MGB Open Access License 1.0 (MGBOA 1.0)**,
  matching CRISPRme+ and CRISPRme. Non-commercial, non-revenue-generating academic use is
  permitted; commercial use requires a commercial license (contact lpinello@mgh.harvard.edu /
  rosalba.giugno@univr.it). See `LICENSE`.

## [2.8.2] - 2026-08-14

### Fixed
- `index` (`-index`) search now ignores hidden metadata sidecar files (dotfiles) in the
  index directory (#42). Released to Bioconda as 2.8.2; the in-repo `VERSION`/recipe pins
  were not bumped at the time (version drift) and are reconciled in 2.8.3.

## [2.8.1] - 2026-08-08

### Added
- `--max-edits` option for `searchTST`: cap the TOTAL number of edits
  (mismatches + DNA bulges + RNA bulges) considered during the ternary-search-tree
  search. The cap prunes the search tree in C++ as it descends, rather than
  filtering results after the fact, so bounded searches complete dramatically
  faster (e.g. an otherwise-intractable pamless whole-chromosome search finishes
  in seconds). `-1` — the default, and the value used when the argument is
  omitted — disables the cap, preserving prior behaviour. This backs CRISPRme's
  total-edits control (CRISPRme issue #107).

### Fixed
- Brute-force (unindexed) `search` used the ENTIRE FASTA header line — including
  spaces and description — as the chromosome name, so non-human assemblies (e.g.
  `>NW_012020308.1 Macaca nemestrina …`) put a spaced description in the
  Chromosome column and broke `bedtools getfasta` in the `-scores` step (empty
  sequence -> azimuth assertion). Now the sequence name is the first
  whitespace-delimited token, matching the samtools/bedtools/UCSC convention.
  No-op for single-token headers like `>chr1` (CRISPRitz issue #13).
- Variant-enrichment (`enricher`) allele-frequency / FILTER robustness: read
  `INFO/AF` by exact key (so a pooled `AF` is selected over source-specific
  `AF_*` / `*_AF` fields that may appear earlier in the record), guard the
  allele-frequency count against the multiallelic records that `bcftools merge`
  can emit, and accept both `PASS` and `.` in the FILTER column. This makes
  enrichment of merged multi-source population panels correct (e.g. CRISPRme's
  combined 1000G + HGDP panel).

- `generate-report`: `radar_chart.py` compared a literal with `is`
  (`if '-pop' is sys.argv[:]`), which raised a `SyntaxWarning` and — always being
  False — left the `-pop` population-plot option dead. Use membership (`in`),
  matching the adjacent `-ws` check (CRISPRitz issue #14).

### Changed
- `add-variants` enrichment progress is written to stdout instead of stderr, so
  it no longer trips CRISPRme's `[ -s logerror ]` stage-failure gate.

## [2.8.0] - 2026-08-04

### Added
- Python 3.11 support; stale committed bytecode removed and ignored going
  forward (#17, #18).
- On-target (azimuth) scoring restored on Python 3.11 via CRISPR-HAWK.
- Golden-output search regression tests for Cas9 (NGG) and Cas12a (TTTV).
- Continuous integration: build + minimal search test + AddressSanitizer
  (guarding the #105 heap-overflow regression) + a native arm64 / Apple-Silicon
  build.
- Multi-architecture (amd64 + arm64) Docker image published to Docker Hub on
  release (#32).

### Changed
- `enricher`: byte-identical C++ port of the variant-enrichment step, replacing
  the previous Python implementation.
- `add-variants`: memory-bounded cross-chromosome parallelism built on the
  compiled enricher.
- `azimuth`: matplotlib imports made lazy/optional (scoring no longer requires
  matplotlib).
- Conda recipe: added the `zlib` dependency (the enricher links `-lz` and
  includes `zlib.h`).
- `conda/meta.yaml`: `source` and `home` now point at pinellolab (were stale
  InfOmics URLs) (#33).
- README: documented platform support (Linux native; Docker for macOS + Windows).

### Fixed
- `add-variants`: normalize multi-rsID (`;`-separated) variant IDs — the root
  cause of CRISPRme issue #143 (#35).
- `add-variants`: worker diagnostic messages written to stdout, not stderr (#34).

## [2.7.1] - 2026-08-01

### Fixed
- `searchTST` heap-buffer-overflow in the guide buffer when a degenerate PAM
  motif is combined with bulges, which crashed the search (CRISPRme issue #105).

## [2.7.0] - 2025-09-17

### Added
- Comprehensive, dual-direction bulge reporting: bulges are now reported in both
  the guide RNA (gRNA) and the target sequence, giving complete bulge
  characterization.
- Support for identifying and reporting target sites that begin with a bulge,
  including edge-case coverage at sequence boundaries.

### Fixed
- Resolved scanning issues in the non-bulge search path, improving off-target
  detection completeness.

### Changed
- Broader target coverage and richer bulge information across off-target analysis
  while maintaining API/command-line compatibility.

## [2.6.6] - 2022-09-14

### Changed
- Release cut to connect the project to Zenodo for archival/DOI.

## [2.6.5] - 2022-03-23

### Changed
- Removed the constraint over the PAM-length sequence.

## [2.6.4] - 2022-03-21

### Fixed
- Fixed a bug related to PAM length.

## [2.6.3] - 2021-09-10

### Fixed
- Fixed handling of `N` bases in targets when no bulges are present in the input.

## [2.6.2] - 2021-09-06

### Fixed
- Fixed a bug in PAM reading and searching.

## [2.6.1] - 2021-09-06

### Changed
- Improved brute-force search and added support for longer PAMs in brute force.

## [2.6.0] - 2021-08-24

### Added
- Support for longer PAMs and mismatches within PAMs (beta).

[Unreleased]: https://github.com/pinellolab/CRISPRitz/compare/v2.9.1...HEAD
[2.10.0]: https://github.com/pinellolab/CRISPRitz/compare/v2.9.1...HEAD
[2.9.1]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.9.1
[2.9.0]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.9.0
[2.8.4]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.8.4
[2.8.3]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.8.3
[2.8.2]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.8.2
[2.8.1]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.8.1
[2.8.0]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.8.0
[2.7.1]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.7.1
[2.7.0]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.7.0
[2.6.6]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.6
[2.6.5]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.5
[2.6.4]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.4
[2.6.3]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.3
[2.6.2]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.2
[2.6.1]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.1
[2.6.0]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.0
