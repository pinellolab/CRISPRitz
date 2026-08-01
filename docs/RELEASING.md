# Releasing CRISPRitz

This guide explains how a maintainer cuts a new CRISPRitz release, publishes it to
**Bioconda**, and then repins the downstream **CRISPRme** package to it. For an
automated, step-by-step assistant version, see the `release-crispritz` Claude Code
skill (`.claude/skills/release-crispritz/SKILL.md`); this document is the
human-readable counterpart.

CRISPRitz is a **compiled** (C++/OpenMP) package, which is what makes its release
process different from a pure-Python one: the Bioconda recipe builds it from
source, and that build must succeed on **both** `linux-64` and `linux-aarch64`
before the package can ship.

## Using the release skill (Claude Code)

Day to day, the easiest way to cut a release is the `release-crispritz` skill from
a Claude Code session opened at the repo root:

1. **Prerequisites:** `master` is green and up to date, and `gh` is authenticated.
2. **Invoke it:** type `/release-crispritz`, or just ask in plain language, e.g.
   *"release CRISPRitz 2.7.1"*. The skill then walks the whole flow:
   pre-flight version-consistency check -> changelog -> in-repo version bump ->
   GitHub release/tag -> Bioconda recipe update (compiled: linux-64 + aarch64
   build CI must pass) -> install verification -> CRISPRme repin + release.
3. It **pauses before anything irreversible** (creating the public tag/release, or
   the cross-repo CRISPRme repin) so you can confirm first.

The rest of this document is the manual version of exactly what the skill
automates.

## The release chain

There is one source of truth for a version — the number `X.Y.Z` — and it flows
outward through GitHub, Bioconda, and finally CRISPRme. The order of operations
matters because each stage consumes an artifact produced by the previous one (the
Bioconda recipe needs the sha256 of the GitHub release tarball; CRISPRme needs the
published Bioconda package).

```
                    ┌──────────────────────────────────────────────┐
                    │             version  X.Y.Z                    │
                    └──────────────────────────────────────────────┘
                                        │
         ┌──────────────────────────────┼──────────────────────────────┐
         ▼                              ▼                              ▼
   crispritz.py                   conda/meta.yaml                 git tag vX.Y.Z
   VERSION = "X.Y.Z"              {% set version = "X.Y.Z" %}           │
                                  (in-repo stub; advisory)              ▼
                                                          GitHub release tarball
                                                    .../archive/vX.Y.Z.tar.gz
                                                                        │
                                                            sha256 of that tarball
                                                                        │
                                                                        ▼
                            Bioconda  recipes/crispritz/meta.yaml   (COMPILED)
                            {% set version = "X.Y.Z" %}
                            source.sha256: <hash>
                            build.number: 0
                            ── BUILD CI must PASS on linux-64 AND linux-aarch64 ──
                            ── a Bioconda maintainer merges ──
                                                                        │
                                          crispritz X.Y.Z live on Bioconda
                                                                        │
                    ┌───────────────────────────────────────────────────┘
                    ▼   (CROSS-REPO repin, in CRISPRme)
   CRISPRme Dockerfile:            ARG crispritz_version=X.Y.Z
   CRISPRme Bioconda recipe:       - crispritz X.Y.Z   (exact pin)
                    │
                    ▼
              cut a CRISPRme release (its own release tooling)
```

### The compiled-package caveat (arm64 must build)

Because the Bioconda recipe compiles CRISPRitz from source:

- The recipe declares `build.skip: True  # [osx]` and
  `extra.additional-platforms: [linux-aarch64]` — so CI builds **linux-64 and
  linux-aarch64**, not macOS.
- The autobump bot can open the version-bump PR, but **it will not merge with a
  red build**. The most common blocker is the `linux-aarch64` job (arch-specific
  compiler/OpenMP issues). A Bioconda maintainer merges once both builds are green.
- This is why a fresh version is not "released" the moment the tag exists — it is
  released when the compiled package for both arches is published on the channel.
  CRISPRme must not be repinned until then.

### The in-repo `conda/meta.yaml` is NOT the published recipe

The authoritative Bioconda recipe lives in
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes) at
`recipes/crispritz/meta.yaml`. The `conda/meta.yaml` in *this* repo is a
historical stub and has drifted from what Bioconda actually ships (different
version, different `source.url` host). Treat the in-repo file as advisory: the
helper script bumps it for tidiness, but the real recipe edit happens in
bioconda-recipes. Likewise, historically the `crispritz.py` `VERSION` string has
lagged the released tag — reconcile all of these in the pre-flight step.

## The version-sync points

| Location | What to change |
|---|---|
| `crispritz.py` | `VERSION = "X.Y.Z"` (surfaced by `crispritz.py version` / `--version`) |
| `conda/meta.yaml` | `{% set version = "X.Y.Z" %}` (in-repo stub; advisory) |
| Git tag / GitHub release | tag `vX.Y.Z`, notes from the changelog |
| `recipes/crispritz/meta.yaml` (in `bioconda/bioconda-recipes`) | `{% set version = "X.Y.Z" %}`, `source.sha256`, `build.number: 0` |
| **CRISPRme** `Dockerfile` | `ARG crispritz_version=X.Y.Z` |
| **CRISPRme** `recipes/crisprme/meta.yaml` (Bioconda) | `- crispritz X.Y.Z` (exact pin) |

The Bioconda `source.url` is templated
(`.../archive/v{{ version }}.tar.gz` — note `v{{ version }}`, no `refs/tags/`) and
must not be edited by hand — bumping `version` repoints it automatically.

## Step-by-step

### 1. Pre-flight consistency check

Confirm the repo, the latest tag, and the live recipe agree on the *current*
(previous) version. If a file lags (a known historical drift), fix it first so you
start from a consistent state.

```bash
grep -E '^VERSION *= *' crispritz.py
grep -E 'set version' conda/meta.yaml
git tag --sort=-v:refname | grep -E '^v?[0-9]' | head -1
curl -sL https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/crispritz/meta.yaml \
  | grep -E 'set version|sha256|number:|additional-platforms|linux-aarch64'
```

### 2. Bump the in-repo version pins

The helper edits only `crispritz.py` and `conda/meta.yaml` (and the `Dockerfile`
if it ever grows a crispritz pin — currently it has none, so it is skipped), is
idempotent, and prints the Bioconda-recipe, changelog, and cross-repo scaffolds:

```bash
python scripts/prepare_release.py X.Y.Z --no-download
git diff crispritz.py conda/meta.yaml     # sanity-check the change
```

(`--no-download` skips the sha256 fetch, which cannot succeed until the tag
exists. Re-run without it in step 5.)

### 3. Update the changelog

CRISPRitz follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Move
entries from `[Unreleased]` into a new `## [X.Y.Z] - YYYY-MM-DD` section, grouped
under Added / Changed / Deprecated / Removed / Fixed / Security. Leave a fresh
empty `[Unreleased]` at the top and update the link-reference footer.

### 4. Create the GitHub release and tag

```bash
git add crispritz.py conda/meta.yaml CHANGELOG.md
git commit -m "Release vX.Y.Z"
git push origin HEAD

gh release create vX.Y.Z \
  --repo pinellolab/CRISPRitz \
  --title "CRISPRitz vX.Y.Z" \
  --notes-file <changelog-section-for-X.Y.Z> \
  --target master
```

Confirm the tarball is live (HTTP 200) before touching Bioconda:

```bash
curl -sIL https://github.com/pinellolab/CRISPRitz/archive/vX.Y.Z.tar.gz | head -1
```

### 5. Update the Bioconda crispritz recipe (compiled)

Recompute the sha256 from the *published* tarball:

```bash
python scripts/prepare_release.py X.Y.Z        # prints sha256 + the exact meta.yaml diff
# equivalently (mirror the recipe URL form):
curl -sL https://github.com/pinellolab/CRISPRitz/archive/vX.Y.Z.tar.gz | shasum -a 256
```

Get it into `bioconda/bioconda-recipes` (bump `version`, set `sha256`, reset
`build.number` to `0`, keep `additional-platforms: [linux-aarch64]` and
`skip: True # [osx]`):

- **Autobump (preferred).** The bot watches releases and opens
  `Update crispritz to X.Y.Z` PRs, often within a day. Verify it carries the right
  version/sha256/`build.number: 0`, then **watch the build CI on both arches** —
  arm64 is the usual blocker for a compiled package. Once lint + `linux-64` +
  `linux-aarch64` are green and approved, a maintainer comments
  `@BiocondaBot please merge`. Nudge a scan with `@BiocondaBot autobump crispritz`.
- **Manual PR.** Fork, branch from an up-to-date `master`, edit
  `recipes/crispritz/meta.yaml`, and open `Update crispritz to X.Y.Z`. CI builds
  and tests on both arches. `@BiocondaBot please fetch artifacts` confirms the
  arm64 build produced a package.

> Build-number rule: reset `build.number` to `0` for a **new version**. Only
> **increment** it (keeping the version) when you rebuild an existing version
> because its dependencies changed.

### 6. Verify the new version installs on linux-64 + linux-aarch64

After the Bioconda PR merges and the package publishes (can take an hour+):

```bash
conda search -c bioconda 'crispritz=X.Y.Z' --info | grep -E 'version|subdir|platform'
# expect linux-64 AND linux-aarch64 rows

conda create -n crispritz_xyz -c conda-forge -c bioconda crispritz=X.Y.Z -y \
  && conda run -n crispritz_xyz crispritz.py version   # expect CRISPRitz v X.Y.Z
```

### 7. Cross-repo: repin CRISPRme and release it

Only after crispritz `X.Y.Z` is live on Bioconda for the arches CRISPRme builds
on, bump BOTH crispritz pins in CRISPRme and cut a CRISPRme release with its own
tooling:

- CRISPRme `Dockerfile`: `ARG crispritz_version=X.Y.Z`
- CRISPRme Bioconda recipe `recipes/crisprme/meta.yaml`: `- crispritz X.Y.Z`

See CRISPRme's own `release-crisprme` skill / `docs/RELEASING.md` for cutting that
release. Do not repin early: CRISPRme's Docker build and Bioconda solve both fail
on `crispritz=X.Y.Z` until the package is actually published.

## Rollback / troubleshooting

- **Bad in-repo bump:** `git checkout crispritz.py conda/meta.yaml CHANGELOG.md`.
- **Premature/incorrect GitHub release:** `gh release delete vX.Y.Z --repo pinellolab/CRISPRitz --cleanup-tag`, then redo.
- **sha256 mismatch in the Bioconda PR:** GitHub occasionally regenerates archive tarballs; recompute the hash from the live tarball and update `meta.yaml`.
- **arm64 build fails on Bioconda:** the defining compiled-package blocker. Read the `linux-aarch64` log (arch-specific compiler flags / OpenMP linkage are common), fix in the recipe build deps or upstream source, push, re-run CI. Never merge with a red arch build.
- **CRISPRme build fails on `crispritz=X.Y.Z`:** the crispritz Bioconda package for that version/arch has not published yet — wait for step 6, then repin/rebuild.
- **Version lag between files:** always run the step-1 pre-flight; historically `crispritz.py` `VERSION` and the in-repo `conda/meta.yaml` have trailed the released tag.
