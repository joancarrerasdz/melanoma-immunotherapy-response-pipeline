# U07 / P1-10 — Leakage, confounding and batch assessment

## Purpose

This document formalizes the assessment of leakage, clinical confounding,
cohort effects and technical batch information for the strict Week 3 and
Week 4 melanoma anti-PD-1 workflow.

It addresses:

- U07 — leakage/confounding risk evaluated;
- P1-10 — batch effects/confounders considered.

This QA remediation does not modify the frozen Week 3 or Week 4 analytical
results.

## Internal cohort

The internal modelling cohort contains 36 pre-treatment anti-PD-1 samples:

- 22 Responders;
- 14 NonResponders.

The prediction target is treatment response.

Relevant clinical variables available in the authoritative clinical table
include:

- Gender;
- Age at treatment;
- BRAF mutation;
- LDH high;
- Biopsy site;
- Biopsy type;
- prior MAPKi exposure;
- prior ipilimumab exposure;
- treatment-history variables.

These variables are treated as possible sources of clinical heterogeneity
or confounding. They are not silently assumed to be balanced.

## Existing descriptive assessment

Clinical distributions were assessed before predictive modelling and are
recorded in:

- `results/day4_clinical_by_response.csv`
- `results/day4_age_by_response.csv`
- `results/day4_clinical_missingness.csv`
- `docs/day4_qc_anomalies.md`

### Gender

Gender is exactly balanced within both response groups:

- NonResponders: 7 female / 7 male;
- Responders: 11 female / 11 male.

No obvious response-class imbalance is therefore observed for Gender.

### Age

Age is similar descriptively between response groups:

- NonResponders: mean 59.64 years, median 60;
- Responders: mean 61.41 years, median 67.

The cohort is small and these summaries are not interpreted as evidence
that age cannot confound expression or response.

### BRAF mutation

BRAF mutation is present in:

- 5/14 NonResponders;
- 6/22 Responders.

Missing BRAF information is also present in both response groups.

The available descriptive evidence does not indicate a major separation by
BRAF status, but residual confounding cannot be excluded.

### LDH

High LDH is present in:

- 6/14 NonResponders;
- 4/22 Responders.

This represents a clinically relevant imbalance and is retained as a
potential source of confounding.

LDH is not retrospectively introduced into the frozen Week 3 model because
doing so after model development and external evaluation would redefine the
analytical procedure.

### Biopsy site and biopsy type

Biopsy site and biopsy type are heterogeneous.

Several biopsy-type categories contain very small numbers of samples and
raw category naming also contains capitalization or terminology variation.

With only 36 samples, stable multivariable adjustment across these sparse
categories would be poorly supported and could introduce substantial model
instability.

These variables are therefore documented as sources of residual clinical
heterogeneity rather than retrospectively incorporated into the frozen
model.

## Technical batch information

No explicit technical batch variable is available in the processed
clinical or GEO metadata tables used by the strict workflow.

In particular, the available metadata do not provide a validated modelling
variable for sequencing batch, lane, plate or processing centre.

Therefore technical batch effects cannot be directly estimated or adjusted
from the available project metadata.

Absence of a recorded batch variable must not be interpreted as evidence
that technical batch effects are absent.

This is retained as a methodological limitation.

## Leakage controls

The strict workflow controls the main identified sources of information
leakage.

Week 3 ensures that:

- outer-test samples do not drive supervised preprocessing;
- feature selection is performed inside resampling;
- hyperparameter tuning is restricted to training data;
- outer-test predictions are generated only after model development inside
  the corresponding training partition.

Week 4 ensures that:

- the 12-gene candidate set is frozen before external evaluation;
- the cross-dataset representation is selected without external outcomes;
- the deployment model is frozen before external response labels are opened;
- all 76 external predictions are frozen before unblinding;
- no post-unblinding gene reselection is allowed;
- no post-unblinding representation change is allowed;
- no post-unblinding model refitting or hyperparameter tuning is allowed;
- no external-label threshold optimization is allowed.

## External-cohort effects

The two external cohorts are treated as separate primary validation cohorts:

- GSE91061;
- GSE78220.

Dataset-specific performance is reported as the primary external-validation
result.

The pooled 76-sample result is retained only as a secondary descriptive
summary because the external cohorts differ in response prevalence and
cohort composition.

The Week 4 `samplewise_rank` representation reduces dependence on absolute
cross-dataset expression scale and is selected without outcome labels.

It must not be interpreted as proving removal of all platform, cohort or
technical batch effects.

## Methodological boundary

This assessment is retrospective QA documentation of the already completed
strict workflow.

It does not:

- change the Week 3 preprocessing procedure;
- alter nested cross-validation;
- reselect genes;
- change the frozen 12-gene candidate set;
- refit the Week 4 deployment model;
- change `mtry`;
- change the classification threshold;
- recompute external predictions using outcome information.

Potential residual confounding is reported as a limitation rather than
removed by post-hoc outcome-driven model modification.

## QA conclusion

The principal leakage risks are explicitly controlled by the strict Week 3
and locked Week 4 designs.

Relevant available clinical variables have been inspected descriptively.

Technical batch cannot be directly assessed because a validated batch
variable is not available in the project metadata.

LDH, biopsy characteristics, treatment history and other unmeasured or
sparsely observed variables remain possible sources of residual confounding.

The external-cohort effect is addressed by reporting cohort-specific
validation as primary and pooled validation only as secondary.

Therefore U07 and P1-10 are considered satisfied as an explicit
risk-assessment and documentation requirement.

This conclusion does not claim that confounding or batch effects have been
eliminated.
