# Azimuth (Doench 2016) on-target scoring — provenance & attribution

This `azimuth/` package provides CRISPRitz's on-target (Rule Set 2 / Doench 2016)
efficiency scoring.

## Why this was re-vendored

The azimuth package previously shipped with CRISPRitz depended on scikit-learn
~0.21. Its trained model pickle (`V3_model_nopos.pickle`) referenced sklearn
module paths that were removed in scikit-learn 0.22 — e.g.
`sklearn.ensemble.gradient_boosting` and `sklearn.tree.tree`. As a result the
model failed to unpickle on modern Python/sklearn with a `ModuleNotFoundError`,
and on-target scoring was broken on Python 3.11 (issues #17/#18). PR #21
temporarily guarded the failure so that CFD (off-target) scoring kept working
while Doench scores degraded to 0.

## What was adopted

This package (code + trained models) is vendored from **CRISPR-HAWK**
(https://github.com/pinellolab/CRISPR-HAWK), the same lab's tool, which
maintains a modern, working azimuth port. Its `V3_model_*.pickle` models were
re-serialized against a modern scikit-learn and unpickle cleanly on
`scikit-learn==1.1.3` (trained with 1.1.1) / `numpy==1.24.4`, referencing the
modern module paths `sklearn.ensemble._gb`, `sklearn.tree._classes`, etc.

Reused with permission (same lab). The trained model artifacts
(`saved_models/*.pickle`) are distributed by CRISPR-HAWK via Zenodo
(record 19680449, `azimuth.zip`).

## Upstream license

The azimuth algorithm and original source are Copyright (c) 2015, Microsoft
Research — see `LICENSE.txt` in this directory. CRISPR-HAWK preserves that
license; CRISPRitz preserves it here unchanged.

## Runtime pins

The trained models load correctly with:

- `scikit-learn == 1.1.3`
- `numpy == 1.24.4`
- `scipy == 1.10.1`
- `pandas == 2.0.3`

These are the versions verified to unpickle the model and reproduce
CRISPR-HAWK's predictions on Python 3.11. Note: scikit-learn 1.2 removed
`sklearn.ensemble._gb_losses`, which the model references, so do not bump
scikit-learn past 1.1.x without re-serializing the model.
