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
