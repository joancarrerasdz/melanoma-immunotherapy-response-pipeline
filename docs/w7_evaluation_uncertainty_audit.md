# W7 — Evaluation and uncertainty audit

## Status

W7 is an additive post-hoc audit of frozen Week 3 and Week 4 predictions and metrics.

## Objective

Quantify predictive-performance variability and uncertainty without retraining models, reselection of genes, threshold optimization, cohort modification or retrospective alteration of frozen analytical outputs.

## Frozen inputs

- `results/week3_nested_cv_summary_strict.csv`
- `results/week3_outer_fold_metrics_strict.csv`
- `results/week3_outer_predictions_strict.csv`
- `results/week4_external_validation_metrics_strict.csv`
- `results/week4_external_validation_predictions_strict.csv`

Their SHA256 digests are recorded in `results/w7_frozen_input_manifest.csv`.

## Non-negotiable constraints

- No model retraining.
- No feature or gene reselection.
- No modification of the frozen Week 3 nested-CV design.
- No modification of the locked Week 4 deployment model.
- No threshold retuning.
- No addition, removal or relabelling of samples.
- No retrospective modification of W1-W4 analytical outputs.

## W7.3 outer-fold variability

The five strict outer folds are summarized using mean, sample standard deviation, median, minimum, maximum and range.

These quantities describe empirical between-fold variability only. They are not confidence intervals and are not treated as five independent estimates of a population parameter.

The pooled out-of-fold metric is reported separately because a metric computed after pooling all OOF predictions is not generally identical to the unweighted mean of fold-specific metrics.

## Planned W7 analyses

1. Frozen-input provenance and checksum manifest.
2. Outer-fold performance variability.
3. No-skill and prevalence-aware reference audit.
4. Prediction-level metric uncertainty.
5. Internal-versus-external transportability audit.
6. Final W7 QA and clean-room reproduction.

## Interpretation boundary

W7 characterizes uncertainty in already-frozen results. It does not attempt to improve predictive performance.


## W7.4 reference baselines

Frozen predictive performance is compared with simple non-informative or prevalence-aware references.

- Accuracy: majority-class classifier.
- Balanced accuracy: 0.50.
- ROC AUC: 0.50.
- Brier score: constant predictor equal to observed responder prevalence.

For Brier score, lower values are better; for the remaining metrics, higher values are better.

Results are recorded in `results/w7_reference_baseline_audit.csv`.

## W7.4 prediction-level uncertainty

Non-parametric bootstrap confidence intervals are calculated from already-frozen prediction records.

Week 3 bootstrap samples are stratified by outcome. Week 4 cohort-specific samples are stratified by outcome, and the pooled Week 4 analysis is stratified jointly by dataset and outcome.

The bootstrap uses 5,000 deterministic replicates with seed 20261016.

These intervals characterize evaluation uncertainty conditional on the frozen predictions. They do not include model-training, feature-selection or hyperparameter-selection uncertainty.

Results are recorded in `results/w7_metric_uncertainty.csv`.

