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

[Unreleased]: https://github.com/pinellolab/CRISPRitz/compare/v2.8.1...HEAD
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
