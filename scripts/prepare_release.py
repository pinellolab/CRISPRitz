#!/usr/bin/env python3
"""prepare_release.py - Prepare a CRISPRitz release across all version-pinned files.

Given a target version (e.g. 2.7.1), this script:

  1. Bumps ``VERSION = "X.Y.Z"`` in ``crispritz.py``.
  2. Bumps ``{% set version = "X.Y.Z" %}`` in ``conda/meta.yaml`` (the in-repo
     recipe stub -- see the note below).
  3. Bumps a version pin in the ``Dockerfile`` if one is present
     (``ARG crispritz_version=X.Y.Z`` or ``crispritz=X.Y.Z``); if the Dockerfile
     carries no explicit crispritz version, it is left untouched.
  4. Downloads the GitHub *release* tarball for the tag ``vX.Y.Z`` and computes
     its sha256 (this is what the **Bioconda** crispritz recipe pins as
     ``source.sha256``).
  5. Prints the exact ``meta.yaml`` changes needed for the Bioconda recipe.
  6. Prints a ``CHANGELOG.md`` scaffold entry for the new version.

CRISPRitz is a COMPILED (C++/OpenMP) Bioconda package. The authoritative recipe
lives in ``bioconda/bioconda-recipes`` at ``recipes/crispritz/meta.yaml`` -- NOT
in this repo. The in-repo ``conda/meta.yaml`` is a historical stub and may be out
of sync with what Bioconda actually ships; this script bumps it for tidiness but
the real recipe must be edited in the bioconda-recipes repo (or via the autobump
bot). See ``docs/RELEASING.md`` and the ``release-crispritz`` skill.

The script is IDEMPOTENT: re-running it for the same version makes no further
edits (it reports "already at X.Y.Z"). It only ever touches the in-repo version
files; it never edits the Bioconda recipe or the changelog for you -- those are
printed for you to apply/verify manually (the recipe lives in a different repo,
and changelog prose needs a human).

It does NOT create git tags, GitHub releases, or Bioconda PRs.

Usage:
    python scripts/prepare_release.py 2.7.1
    python scripts/prepare_release.py 2.7.1 --no-download   # skip sha256 fetch
    python scripts/prepare_release.py 2.7.1 --repo-root /path/to/CRISPRitz
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
import urllib.request

REPO = "pinellolab/CRISPRitz"
# NOTE: the Bioconda crispritz recipe uses the .../archive/vX.Y.Z.tar.gz form
# (no "refs/tags/"). We match it exactly so the sha256 we compute is the sha256
# Bioconda will see. Both URL forms resolve to the same bytes, but we mirror the
# recipe to be unambiguous.
TARBALL_URL = "https://github.com/{repo}/archive/v{version}.tar.gz"

VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")


def die(msg: str) -> "None":
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def bump_crispritz_py(path: str, version: str) -> bool:
    """Set ``VERSION = "X.Y.Z"`` in crispritz.py. Returns True if changed."""
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    pat = re.compile(r'^(VERSION\s*=\s*)"(\d+\.\d+\.\d+)"', re.MULTILINE)
    m = pat.search(text)
    if not m:
        die(f'could not find `VERSION = "..."` in {path}')
    current = m.group(2)
    if current == version:
        print(f"  crispritz.py       already at {version} (no change)")
        return False
    new_text = pat.sub(lambda _m: f'{_m.group(1)}"{version}"', text, count=1)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new_text)
    print(f"  crispritz.py       {current} -> {version}")
    return True


def bump_meta_yaml(path: str, version: str) -> bool:
    """Set ``{% set version = "X.Y.Z" %}`` in conda/meta.yaml. Returns True if changed."""
    if not os.path.isfile(path):
        print(f"  conda/meta.yaml    not found -- skipped")
        return False
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    pat = re.compile(r'(\{%\s*set\s+version\s*=\s*")(\d+\.\d+\.\d+)("\s*%\})')
    m = pat.search(text)
    if not m:
        print(f"  conda/meta.yaml    no `set version` line -- skipped")
        return False
    current = m.group(2)
    if current == version:
        print(f"  conda/meta.yaml    already at {version} (no change)")
        return False
    new_text = pat.sub(lambda _m: f"{_m.group(1)}{version}{_m.group(3)}", text, count=1)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new_text)
    print(f"  conda/meta.yaml    {current} -> {version}  (in-repo stub; see Bioconda note)")
    return True


def bump_dockerfile(path: str, version: str) -> bool:
    """Set a crispritz version pin in the Dockerfile if one exists.

    Matches either ``ARG crispritz_version=X.Y.Z`` or an inline ``crispritz=X.Y.Z``.
    If the Dockerfile carries no explicit crispritz version, it is left untouched
    (the CRISPRitz Dockerfile historically installs the latest package without a
    pin). Returns True if changed.
    """
    if not os.path.isfile(path):
        print(f"  Dockerfile         not found -- skipped")
        return False
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    arg_pat = re.compile(r"^(ARG\s+crispritz_version=)(\d+\.\d+\.\d+)", re.MULTILINE)
    pin_pat = re.compile(r"(crispritz=)(\d+\.\d+\.\d+)")
    for pat, label in ((arg_pat, "ARG crispritz_version"), (pin_pat, "crispritz=")):
        m = pat.search(text)
        if m:
            current = m.group(2)
            if current == version:
                print(f"  Dockerfile         already at {version} (no change)")
                return False
            new_text = pat.sub(lambda _m: f"{_m.group(1)}{version}", text, count=1)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(new_text)
            print(f"  Dockerfile         {current} -> {version}  (via {label})")
            return True
    print(f"  Dockerfile         no crispritz version pin found -- skipped")
    return False


def compute_sha256(version: str) -> str:
    url = TARBALL_URL.format(repo=REPO, version=version)
    print(f"\nFetching release tarball to compute sha256:\n  {url}")
    try:
        with urllib.request.urlopen(url) as resp:  # noqa: S310 (trusted host)
            data = resp.read()
    except Exception as exc:  # noqa: BLE001
        die(
            f"could not download tarball for v{version}: {exc}\n"
            "       (Has the GitHub release/tag been created yet? "
            "The sha256 can only be computed AFTER the tag exists.)"
        )
    digest = hashlib.sha256(data).hexdigest()
    print(f"  bytes: {len(data)}")
    print(f"  sha256: {digest}")
    return digest


def print_meta_yaml_changes(version: str, sha256: str | None) -> None:
    sha_line = sha256 if sha256 else "<run without --no-download to compute>"
    print(
        "\n"
        "==============================================================\n"
        "  Bioconda recipe changes  (recipes/crispritz/meta.yaml)\n"
        "  Repo: https://github.com/bioconda/bioconda-recipes\n"
        "  NOTE: CRISPRitz is COMPILED -- the recipe builds from source on\n"
        "        linux-64 AND linux-aarch64. The autobump bot opens the PR,\n"
        "        but the Bioconda build CI must pass on BOTH arches and a\n"
        "        Bioconda maintainer merges it.\n"
        "==============================================================\n"
        f'  {{% set version = "{version}" %}}      # bump version\n'
        f"  source:\n"
        f"    sha256: {sha_line}\n"
        f"  build:\n"
        f"    number: 0                           # RESET to 0 for a new version\n"
        "\n"
        "  NOTE: source.url uses v{{ version }}, so it auto-resolves to\n"
        f"        .../archive/v{version}.tar.gz -- do not edit it.\n"
        "  NOTE: if this is a REBUILD of the same version (deps changed only),\n"
        "        keep version the same and INCREMENT build.number instead.\n"
    )


def print_changelog_scaffold(version: str) -> None:
    import datetime

    today = datetime.date.today().isoformat()
    print(
        "\n"
        "==============================================================\n"
        "  CHANGELOG.md scaffold  (move items out of [Unreleased])\n"
        "==============================================================\n"
        f"## [{version}] - {today}\n"
        "### Added\n-\n"
        "### Changed\n-\n"
        "### Fixed\n-\n"
        "\n"
        "  Also update the link-reference footer, e.g.:\n"
        f"  [{version}]: https://github.com/{REPO}/releases/tag/v{version}\n"
    )


def print_crossrepo_reminder(version: str) -> None:
    print(
        "\n"
        "==============================================================\n"
        "  Cross-repo repin  (AFTER crispritz {v} is live on Bioconda)\n"
        "==============================================================\n"
        "  CRISPRme pins crispritz EXACTLY in two places -- bump both, then\n"
        "  cut a CRISPRme release:\n"
        "    - CRISPRme Dockerfile:            ARG crispritz_version={v}\n"
        "    - CRISPRme Bioconda recipe:       - crispritz {v}\n"
        "  (recipes/crisprme/meta.yaml in bioconda-recipes)\n".format(v=version)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("version", help="target version, e.g. 2.7.1 (no leading v)")
    parser.add_argument(
        "--repo-root",
        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        help="path to the CRISPRitz repo root (default: parent of scripts/)",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="skip downloading the tarball / computing sha256",
    )
    args = parser.parse_args()

    version = args.version.lstrip("v")
    if not VERSION_RE.match(version):
        die(f"version must look like X.Y.Z (got {args.version!r})")

    crispritz_py = os.path.join(args.repo_root, "crispritz.py")
    meta_yaml = os.path.join(args.repo_root, "conda", "meta.yaml")
    dockerfile = os.path.join(args.repo_root, "Dockerfile")
    if not os.path.isfile(crispritz_py):
        die(f"expected file not found: {crispritz_py}")

    print(f"Preparing CRISPRitz release v{version}\n")
    print("Bumping in-repo version pins:")
    changed = False
    changed |= bump_crispritz_py(crispritz_py, version)
    changed |= bump_meta_yaml(meta_yaml, version)
    changed |= bump_dockerfile(dockerfile, version)
    if not changed:
        print("  (all in-repo files already at target version)")

    sha256 = None if args.no_download else compute_sha256(version)

    print_meta_yaml_changes(version, sha256)
    print_changelog_scaffold(version)
    print_crossrepo_reminder(version)

    print(
        "Next steps (see docs/RELEASING.md or the release-crispritz skill):\n"
        "  1. Review the git diff of the in-repo version files.\n"
        "  2. Move [Unreleased] changelog items into the new version section.\n"
        "  3. Commit, then create the GitHub release/tag vX.Y.Z.\n"
        "  4. Update the Bioconda crispritz recipe (autobump PR or manual) with the\n"
        "     sha256 above; wait for the linux-64 + linux-aarch64 build CI to pass.\n"
        "  5. Once crispritz vX.Y.Z is live on Bioconda, repin it in CRISPRme\n"
        "     (Dockerfile ARG + Bioconda recipe) and cut a CRISPRme release.\n"
    )


if __name__ == "__main__":
    main()
