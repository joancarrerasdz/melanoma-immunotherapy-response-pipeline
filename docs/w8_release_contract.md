# W8 final reproducibility certification contract

## Purpose

Week 8 certifies the final reproducibility, integrity and release readiness
of the P1 melanoma immunotherapy response pipeline.

W8 is a certification and release phase. It does not introduce a new model,
new feature selection, threshold retuning, new cohort or new biological
interpretation.

## Canonical baseline

Baseline branch: main

Baseline commit:

2652c75a3f53a6f94a743245e420a3975e037bdd

This commit represents the canonical post-W7 analytical state.

## Canonical artifact manifest

Number of baseline files protected by the manifest:

222

SHA-256 fingerprint of the canonical manifest:

79de17fe373e0664f1e06165dfd073d0960585bb3b3b9b2091bdb1a023644841

The canonical file-level hashes are stored in:

results/w8_canonical_artifact_manifest.sha256

A machine-readable representation is stored in:

results/w8_canonical_artifact_manifest.csv

## Protected baseline

The manifest covers the canonical pre-W8 files in:

- data_processed
- results
- figures
- docs
- scripts
- renv.lock
- DESCRIPTION
- .Rprofile
- .github/workflows

Pre-W8 baseline files are treated as immutable during certification unless
a reproducibility defect requires an explicitly documented remediation.

New W8 certification artifacts may be added without altering the frozen
W1-W7 analytical conclusions.

## W8 certification gates

The final release requires:

1. Canonical SHA-256 manifest verification.
2. Fresh-clone and clean-room environment restoration.
3. End-to-end pipeline regeneration.
4. Reproducibility of protected analytical artifacts.
5. Repository and documentation consistency.
6. Repository QA success.
7. Clean Git worktree.
8. Local and remote main synchronization.
9. Final tagged portfolio release.

Any failed gate requires review before release certification.

## W8 reproducibility remediation

During the W8 clean-room end-to-end regeneration audit, the protected
Week 2 exploratory PCA numerical outputs reproduced exactly, but the
historical PNG raster did not reproduce byte-for-byte.

The historical `figures/week2_pca_exploratory.png` predates the formal
locked `renv` reproducibility environment. Investigation showed that the
difference was restricted to graphics rendering; PCA scores, explained
variance and other protected analytical outputs remained unchanged.

A controlled graphics-backend audit demonstrated deterministic repeated
rendering with the Cairo PNG backend at the existing 2400 x 1800 pixel,
300 dpi RGB output contract.

Therefore W8 applies one explicit reproducibility remediation:

- `scripts/04_exploratory_analysis.R` explicitly fixes the PCA PNG device
  to `grDevices::png` with `type = "cairo"`.
- `figures/week2_pca_exploratory.png` is regenerated under that locked
  graphics contract.
- No model, cohort, preprocessing logic, PCA values, predictive result,
  threshold, feature selection or biological conclusion is changed.

This remediation is a reproducibility fix rather than an analytical
revision and is permitted by the protected-baseline exception defined
in this contract.
