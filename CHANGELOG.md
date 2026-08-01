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

### Fixed
- `searchTST` heap buffer overflow when a degenerate PAM motif is combined with
  bulges, which crashed the search (issue #105).

### Added
- Continuous-integration workflow that builds CRISPRitz and runs the smoke/e2e
  tests on Linux, catching build and search regressions before release.

### Changed
- Removed stale committed `.pyc` bytecode from the source tree and ignore it
  going forward, so the packaged sources no longer ship out-of-date compiled
  artifacts.

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

[Unreleased]: https://github.com/pinellolab/CRISPRitz/compare/v2.7.0...HEAD
[2.7.0]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.7.0
[2.6.6]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.6
[2.6.5]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.5
[2.6.4]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.4
[2.6.3]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.3
[2.6.2]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.2
[2.6.1]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.1
[2.6.0]: https://github.com/pinellolab/CRISPRitz/releases/tag/v2.6.0
