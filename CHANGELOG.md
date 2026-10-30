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
  12th positional arg to `searchTST`, or the `CRISPRITZ_PAM_PRUNE` env var.
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
  third tests only the positions the real PAM covers. Carrying `N` through also needed two
  encoding defects fixed: the writer shifted its nibble conditionally on the *value*, so a
  nucleotide encoding to `0x0` lost its slot and displaced every following nibble, and the
  loader's `case 0x0` was an unreachable no-op now free to carry `N`. Impact is
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
  hangs if a worker stops returning.

### Known issues
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
