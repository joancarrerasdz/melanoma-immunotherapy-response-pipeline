# W7 — Final evaluation and uncertainty synthesis

## Purpose

Week 7 is a post-hoc evaluation audit of the frozen Week 3 strict internal validation and Week 4 locked external validation.

No model was retrained, no feature set was reselected, no threshold was retuned and no W1-W4 analytical output was modified.

## Internal performance

- W3 OOF ROC AUC: 0.890 (95% bootstrap CI 0.763–0.977).
- W3 OOF balanced accuracy: 0.812 (95% bootstrap CI 0.669–0.929).

These results support strong internal predictive discrimination.

## External transportability

- Pooled external AUC degradation: 0.318 (95% bootstrap CI 0.145–0.478).
- Pooled external balanced-accuracy degradation: 0.264 (95% bootstrap CI 0.098–0.418).
- Pooled external Brier degradation: 0.237 (95% bootstrap CI 0.163–0.310).

The pooled comparison therefore supports degradation of discrimination, balanced classification and probability quality when moving from internal to external evaluation.

## Cohort-specific interpretation

GSE91061 shows supported degradation for balanced accuracy, ROC AUC and Brier score.

GSE78220 shows a poorer external point estimate for all three metrics, but the balanced-accuracy and ROC-AUC gap intervals include zero; Brier degradation remains supported.

## Interpretation boundary

These results demonstrate limited transportability of the frozen predictive performance across the external datasets evaluated here.

They do not establish that dataset shift alone caused the performance differences, and they should not be interpreted as a causal comparison between cohorts.

Differences in cohort composition, responder prevalence, assay context and study-specific characteristics remain plausible contributors to the observed external degradation.

## Portfolio conclusion

The final model shows strong internal discrimination but substantially weaker and more uncertain external performance.

The principal methodological conclusion is therefore not simply that the model performs well internally, but that external transportability is the dominant limitation revealed by the reproducible evaluation.

## Generated artifacts

- `results/w7_final_synthesis.csv`
- `results/w7_key_conclusions.csv`
- `figures/w7_transportability_uncertainty.png`
- `docs/w7_final_synthesis.md`

