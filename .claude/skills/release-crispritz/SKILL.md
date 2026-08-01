---
name: release-crispritz
description: Cut a new CRISPRitz release, publish it to Bioconda, and repin the downstream CRISPRme package to it. Use when the user wants to release, publish, tag, or bump the version of CRISPRitz; update the Bioconda crispritz recipe; update the CRISPRitz changelog; or propagate a new CRISPRitz version into CRISPRme. Covers the version-consistency pre-flight, CHANGELOG.md, the GitHub release/tag, the COMPILED Bioconda recipe update (linux-64 + linux-aarch64 build CI), and the cross-repo CRISPRme repin + release.
---

# Release CRISPRitz

CRISPRitz is a **compiled** (C++/OpenMP) Bioconda package. CRISPRme pins it to an
EXACT version. So a CRISPRitz release drives a chain:

```
git tag vX.Y.Z ─▶ GitHub release tarball ─▶ Bioconda recipes/crispritz/meta.yaml
      │                                        (sha256 of tarball; COMPILED: build
      │                                         CI must pass on linux-64 AND
      │                                         linux-aarch64; maintainer merges)
      │                                                     │
      └── version pinned in crispritz.py + conda/meta.yaml  │  once crispritz X.Y.Z
                                                            ▼  is LIVE on Bioconda:
                                          CRISPRme repin (both pins) ─▶ CRISPRme release
                                            Dockerfile: ARG crispritz_version=X.Y.Z
                                            recipe:     - crispritz X.Y.Z
```

Version-sync points inside CRISPRitz (all must equal the same `X.Y.Z`):
- `crispritz.py` -> `VERSION = "X.Y.Z"` (surfaced by `crispritz.py version` / `--version`)
- `conda/meta.yaml` -> `{% set version = "X.Y.Z" %}` (in-repo stub; the AUTHORITATIVE recipe is on Bioconda)
- git tag / GitHub release -> `vX.Y.Z`
- Bioconda `recipes/crispritz/meta.yaml` -> `{% set version = "X.Y.Z" %}` + fresh `sha256` + `build.number: 0`

The Bioconda crispritz `source.url` is
`https://github.com/pinellolab/CRISPRitz/archive/v{{ version }}.tar.gz`
(note: `v{{ version }}`, NOT `refs/tags/`), so the recipe pulls the **GitHub
release tarball** — its sha256 is the hash of that exact tarball.

> Heads-up: historically CRISPRitz version pins have drifted. On `master` the
> `crispritz.py` `VERSION` and the in-repo `conda/meta.yaml` version can lag the
> actual released tag / Bioconda version. The in-repo `conda/meta.yaml` is NOT
> what Bioconda ships — the real recipe lives in `bioconda/bioconda-recipes`.
> Reconcile these in Step 0 before releasing.

Run all commands from the CRISPRitz repo root. Replace `X.Y.Z` with the real
version (no leading `v` except in tags/URLs).

---

## Step 0 — Pre-flight: version-consistency check

Confirm where the repo currently stands and that nothing is half-bumped.

```bash
echo "crispritz.py: $(grep -E '^VERSION *= *' crispritz.py)"
echo "conda/meta.yaml: $(grep -E 'set version' conda/meta.yaml)"
echo "latest tag:  $(git tag --sort=-v:refname | grep -E '^v?[0-9]' | head -1)"
echo "latest release: $(gh release list --repo pinellolab/CRISPRitz -L1)"
echo "--- live Bioconda crispritz recipe ---"
curl -sL https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/crispritz/meta.yaml \
  | grep -E 'set version|sha256|number:|additional-platforms|linux-aarch64'
```

Expected before a release: `crispritz.py`, `conda/meta.yaml`, the latest tag, and
the Bioconda recipe should all agree on the **previous** version. If they
disagree (e.g. `crispritz.py` `VERSION` lags the released tag — a known drift),
fix the lag as part of this release so you start from a consistent state.
**Do not proceed with inconsistent state you don't understand.**

## Step 1 — Bump in-repo version pins (helper script)

The helper edits ONLY `crispritz.py` and `conda/meta.yaml` (and the `Dockerfile`
if it carries a crispritz pin — the current one does not, so it is skipped). It is
idempotent and prints the Bioconda-recipe + changelog + cross-repo scaffolds. The
sha256 fetch only works AFTER the tag exists, so run it now with `--no-download`,
then again after Step 3 to get the hash.

```bash
python scripts/prepare_release.py X.Y.Z --no-download
git --no-pager diff --stat        # expect ONLY crispritz.py (+ conda/meta.yaml)
git --no-pager diff crispritz.py conda/meta.yaml
```

Rollback if wrong: `git checkout crispritz.py conda/meta.yaml`.

## Step 2 — Update CHANGELOG.md

Move items from `## [Unreleased]` into a new dated section, using the script's
scaffold as a template. Keep-a-Changelog categories: Added / Changed / Deprecated
/ Removed / Fixed / Security.

```bash
python scripts/prepare_release.py X.Y.Z --no-download | sed -n '/CHANGELOG.md scaffold/,/releases\/tag/p'
```

Edit `CHANGELOG.md`:
- Add `## [X.Y.Z] - YYYY-MM-DD` with the curated entries.
- Leave a fresh empty `## [Unreleased]` at the top.
- Update the link-reference footer:
  `[X.Y.Z]: https://github.com/pinellolab/CRISPRitz/releases/tag/vX.Y.Z`
  and repoint `[Unreleased]` to `.../compare/vX.Y.Z...HEAD`.

## Step 3 — Commit, then create the GitHub release + tag

Commit the version bump + changelog, push, then create the release.
`gh release create` creates the annotated tag `vX.Y.Z` for you.

```bash
git add crispritz.py conda/meta.yaml CHANGELOG.md
git commit -m "Release vX.Y.Z"
git push origin HEAD

# Notes come straight from the changelog section you just wrote:
awk '/^## \[X.Y.Z\]/{f=1;next} /^## \[/{f=0} f' CHANGELOG.md > /tmp/crispritz_notes.md

gh release create vX.Y.Z \
  --repo pinellolab/CRISPRitz \
  --title "CRISPRitz vX.Y.Z" \
  --notes-file /tmp/crispritz_notes.md \
  --target master
```

Verify the tarball is now published (this is exactly the tarball Bioconda hashes):

```bash
curl -sIL https://github.com/pinellolab/CRISPRitz/archive/vX.Y.Z.tar.gz | head -1  # expect 200
```

Rollback: `gh release delete vX.Y.Z --repo pinellolab/CRISPRitz --cleanup-tag`
(only if the release must be withdrawn before Bioconda picks it up).

## Step 4 — Update the Bioconda crispritz recipe (COMPILED)

First compute the sha256 of the now-published tarball (this is the recipe's
`source.sha256`):

```bash
python scripts/prepare_release.py X.Y.Z    # no --no-download; prints sha256 + meta.yaml diff
# or directly (mirror the recipe's URL form, no refs/tags/):
curl -sL https://github.com/pinellolab/CRISPRitz/archive/vX.Y.Z.tar.gz | shasum -a 256
```

> **CRISPRitz is compiled.** Unlike a pure-Python package, the Bioconda recipe
> builds CRISPRitz from source. The build CI must PASS on **both** `linux-64` and
> `linux-aarch64` (the recipe declares `extra.additional-platforms: linux-aarch64`
> and `skip: True # [osx]`). The autobump bot opens the PR, but a **Bioconda
> maintainer** merges it after that build CI is green.

You have two routes. **Prefer autobump; fall back to a manual PR.**

### Route A — Autobump bot (preferred)

Bioconda runs an autobump bot that periodically scans GitHub releases and opens
`Update crispritz to X.Y.Z` PRs on its own. Wait ~a day, then check:

```bash
gh pr list --repo bioconda/bioconda-recipes --search "crispritz in:title" --state all -L 5
```

If a PR exists, review that it has the right version + sha256 + `build.number: 0`,
then WATCH THE BUILD CI on both arches. Because the package is compiled, arm64
build failures are the most likely blocker — check the `linux-aarch64` job
specifically. Once lint + both builds are green and it is approved, a Bioconda
member (or you, if a recipe maintainer) comments on the PR:

```
@BiocondaBot please merge
```

To nudge the bot to scan immediately instead of waiting, comment on any recipe
PR/issue:

```
@BiocondaBot autobump crispritz
```

### Route B — Manual PR (when you don't want to wait, or arm64 needs a fix)

```bash
# one-time: fork bioconda/bioconda-recipes to your account
gh repo fork bioconda/bioconda-recipes --clone --remote
cd bioconda-recipes
git checkout master && git pull upstream master
git checkout -b crispritz-X.Y.Z
```

Edit `recipes/crispritz/meta.yaml`:
- `{% set version = "X.Y.Z" %}`  <- bump
- `source.sha256:` <- the hash from above
- `build.number: 0`  <- RESET to 0 (only INCREMENT, keeping version, for a same-version rebuild)
- Do NOT edit `source.url` — it uses `v{{ version }}` and resolves automatically.
- Keep `extra.additional-platforms: [linux-aarch64]` and `build.skip: True # [osx]`.

```bash
git commit -am "Update crispritz to X.Y.Z"
git push -u origin crispritz-X.Y.Z
gh pr create --repo bioconda/bioconda-recipes \
  --title "Update crispritz to X.Y.Z" \
  --body "Bump crispritz to vX.Y.Z. sha256 recomputed from the GitHub release tarball. build.number reset to 0. Compiled package: needs green build on linux-64 + linux-aarch64."
```

Useful @BiocondaBot PR comments:
- `@BiocondaBot please update` — sync your PR branch with upstream master.
- `@BiocondaBot please add label` — request the `please review & merge` label.
- `@BiocondaBot please fetch artifacts` — links to the built test packages/containers (use to confirm the arm64 build produced a package).
- `@BiocondaBot please merge` — after lint + BOTH builds are green and the PR is approved.

Wait for lint + `linux-64` + `linux-aarch64` builds to pass. Fix any build/lint
findings before asking for a merge.

## Step 5 — Verify the new version installs on linux-64 + linux-aarch64

Do this AFTER the Bioconda PR merges and the package is on the channel (can take
an hour+). Confirm the compiled package resolves and installs on both arches.

```bash
# The version is now discoverable and carries both arch builds:
conda search -c bioconda 'crispritz=X.Y.Z' --info | grep -E 'version|subdir|platform'
# expect linux-64 AND linux-aarch64 rows

# Actually install + smoke-test (native or emulated):
#   linux-64:
conda create -n crispritz_xyz_x64  -c conda-forge -c bioconda crispritz=X.Y.Z -y \
  && conda run -n crispritz_xyz_x64 crispritz.py version
#   linux-aarch64 (native arm64 host, or under emulation):
conda create -n crispritz_xyz_arm  -c conda-forge -c bioconda --platform linux-aarch64 crispritz=X.Y.Z -y \
  && conda run -n crispritz_xyz_arm crispritz.py version
```

Both should report `CRISPRitz v X.Y.Z`.

## Step 6 — Cross-repo: repin CRISPRme and cut a CRISPRme release

Once crispritz `X.Y.Z` is LIVE on Bioconda (Step 5 passed), CRISPRme must be
repinned. CRISPRme pins crispritz in **two** places — bump BOTH:

1. **CRISPRme `Dockerfile`** — `ARG crispritz_version=X.Y.Z`
   (installed via `micromamba install ... crispritz=$crispritz_version ...`).
2. **CRISPRme Bioconda recipe** — `recipes/crisprme/meta.yaml` `run:` requirement
   `- crispritz X.Y.Z` (an exact pin).

```bash
# In the CRISPRme repo:
grep -nE 'ARG crispritz_version=' Dockerfile          # bump to X.Y.Z
# In bioconda/bioconda-recipes:
grep -nE '^\s*- crispritz ' recipes/crisprme/meta.yaml # bump to X.Y.Z, reset crisprme build.number: 0
```

Then cut a CRISPRme release using CRISPRme's own tooling (its `release-crisprme`
skill / `docs/RELEASING.md`). The CRISPRme release order is the usual
GitHub tag -> Bioconda recipe (sha256) -> Docker; the crispritz repin rides along
in that same CRISPRme release.

> Order matters: do NOT repin CRISPRme to crispritz `X.Y.Z` until crispritz
> `X.Y.Z` is actually published on Bioconda for the arches CRISPRme builds on —
> otherwise CRISPRme's Docker build and Bioconda solve will fail on
> `crispritz=X.Y.Z` (unresolved package).

---

## Failure / rollback quick reference
- Wrong in-repo bump: `git checkout crispritz.py conda/meta.yaml CHANGELOG.md`.
- Wrong/premature GitHub release: `gh release delete vX.Y.Z --repo pinellolab/CRISPRitz --cleanup-tag` and redo.
- sha256 mismatch on the Bioconda PR: GitHub can regenerate archive tarballs; recompute the sha256 from the live tarball and update `meta.yaml`, then push.
- **arm64 build fails on Bioconda:** this is the compiled-package blocker. Inspect the `linux-aarch64` build log; common causes are arch-specific compiler flags / OpenMP linkage. Fix in the recipe (build deps) or upstream (source), push, re-run CI. Do NOT merge with a red arch build.
- Deps changed but version unchanged: keep `version`, INCREMENT `build.number` (do not reset to 0).
- CRISPRme build fails on `crispritz=X.Y.Z`: the crispritz Bioconda package for that version/arch isn't published yet — wait for Step 5 to pass, then repin/rebuild.
