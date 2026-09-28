# Melanoma immunotherapy response — reproducible bioinformatics pipeline

A reproducible transcriptomic machine-learning workflow for evaluating
anti-PD-1 treatment response in melanoma.

This repository reconstructs an original MSc thesis analysis under a stricter
methodological framework focused on:

- traceable clinical endpoint reconstruction;
- leakage-aware preprocessing;
- nested cross-validation;
- feature-selection stability;
- frozen model specification;
- blind external prediction;
- locked external validation;
- reproducible computing environments;
- explicit methodological QA gates.

The original academic repository state is preserved at tag `tfm-original`.

The strict Week 1–4 methodological rebuild and the subsequent Week 3–4 QA
remediation are complete and merged into `main`.

---

## Project overview

### Development cohort

**GSE160638**

- 36 anti-PD-1 melanoma samples;
- 22 Responders;
- 14 NonResponders;
- 17,002 validated expression features.

The clinical treatment-response endpoint is reconstructed from the official
supplementary clinical Table S2A.

Response mapping:

- `CR` / `PR` → `Responder`
- `SD` / `PD` → `NonResponder`

Full data provenance is documented in [`DATA.md`](DATA.md).

### External cohorts

The frozen Week 4 procedure is evaluated independently on:

- **GSE91061:** 49 pre-treatment samples;
- **GSE78220:** 27 pre-treatment samples.

External response labels are not used to redefine the gene set, expression
representation, model, hyperparameters or classification threshold.

---

## Methodological workflow

~~~text
Validated raw counts + clinical endpoint
                |
                v
Week 2 exploratory preprocessing
(descriptive / exploratory only)
                |
                v
Week 3 strict nested cross-validation
                |
                +--> predictive preprocessing from training data
                |
                +--> feature selection inside resampling
                |
                +--> hyperparameter tuning inside training data
                |
                v
Outer out-of-fold validation
                |
                v
Post-validation gene-stability analysis
                |
                v
Frozen 12-gene consensus candidate set
                |
                v
Week 4 label-free representation audit
                |
                v
Frozen samplewise_rank representation
                |
                v
Frozen model + blind external predictions
                |
                v
External outcome unblinding
                |
                v
Locked external evaluation
~~~

The formal preprocessing boundary is documented in
[`docs/preprocessing_contract.md`](docs/preprocessing_contract.md).

---

## Week 3 — strict internal validation

Week 3 uses nested cross-validation to separate model development from
outer-fold performance estimation.

The strict workflow uses:

- 5 outer folds;
- 3 inner folds;
- feature selection restricted to training partitions;
- hyperparameter tuning restricted to training partitions;
- one outer out-of-fold prediction per sample.

### Internal nested-CV performance

| Metric | Value |
|---|---:|
| Accuracy | 0.8333 |
| Sensitivity | 0.9091 |
| Specificity | 0.7143 |
| Balanced Accuracy | 0.8117 |
| AUC | 0.8896 |

These values are **internal nested-cross-validation estimates**.

They are not external-validation performance estimates.

### Feature-selection stability

| Selection frequency | Genes |
|---|---:|
| ≥1/5 outer folds | 312 |
| ≥2/5 outer folds | 100 |
| ≥3/5 outer folds | 51 |
| ≥4/5 outer folds | 25 |
| 5/5 outer folds | 12 |

The 12 genes selected in all five outer folds are treated as
**stability-derived consensus candidates**.

They are not described as a clinically validated molecular signature.

See
[`docs/week3_historical_vs_strict_stability_audit.md`](docs/week3_historical_vs_strict_stability_audit.md).

---

## Week 4 — locked external validation

Week 4 evaluates a frozen deployment procedure on two independent cohorts.

The locked specification uses:

- the 12 Week 3 consensus candidates;
- `samplewise_rank` representation;
- Random Forest;
- `mtry = 1`;
- 500 trees;
- classification threshold `0.48`.

The deployment model and all 76 blind external predictions were frozen before
external response labels were opened.

Blind-prediction freeze commit:

`3e427be1b885b9514cb8acc857e3eb3ca90e4ceb`

### Primary cohort-specific results

| Dataset | n | AUC | Accuracy | Sensitivity | Specificity | Balanced Accuracy |
|---|---:|---:|---:|---:|---:|---:|
| GSE91061 | 49 | 0.6641 | 0.3265 | 0.9000 | 0.1795 | 0.5397 |
| GSE78220 | 27 | 0.6444 | 0.6667 | 0.8000 | 0.5000 | 0.6500 |

Dataset-specific results are the primary external-validation assessment.

### Secondary pooled summary

| Metric | Value |
|---|---:|
| AUC | 0.5718 |
| Accuracy | 0.4474 |
| Sensitivity | 0.8400 |
| Specificity | 0.2549 |
| Balanced Accuracy | 0.5475 |

The pooled result is retained only as a **secondary descriptive summary**
because the two external cohorts differ in prevalence and composition.

Detailed methodology is documented in
[`docs/week4_locked_external_validation.md`](docs/week4_locked_external_validation.md).

---

## Interpretation

The reconstructed analysis shows an important distinction between
**internal predictive signal** and **external transportability**.

The internal nested-CV AUC is:

`0.8896`

The two independent external-cohort AUCs are:

- GSE91061: `0.6641`
- GSE78220: `0.6444`

The difference is preserved rather than optimized away after external
unblinding.

Accordingly:

- Week 3 performance is reported as internal nested-CV performance;
- external validation is reported primarily by cohort;
- the pooled external estimate remains secondary;
- the 12-gene set remains a stability-derived deployment candidate;
- no claim of clinical biomarker validation is made.

---

## Leakage safeguards

Predictive preprocessing, feature selection and tuning are restricted to
training data during model development.

After external unblinding, the workflow explicitly prohibits:

- gene reselection;
- expression-representation changes;
- model refitting;
- hyperparameter retuning;
- classification-threshold changes;
- threshold optimization against external outcomes.

---

## Confounding and cohort heterogeneity

Clinical and methodological heterogeneity were formally assessed during QA
remediation.

The assessment covers:

- response-class composition;
- age;
- gender;
- BRAF mutation;
- LDH;
- biopsy heterogeneity;
- technical batch-metadata availability;
- residual confounding;
- external cohort effects.

No validated technical batch variable was available in the project metadata.

This limitation is documented explicitly and is not interpreted as evidence
that batch effects are absent.

See
[`docs/confounding_assessment.md`](docs/confounding_assessment.md).

---

## Reproducibility

The strict Week 3–4 workflow has a project-level reproducible R environment.

Key files:

- `renv.lock`
- `DESCRIPTION`
- `.Rprofile`
- `renv/`
- [`docs/reproducible_environment.md`](docs/reproducible_environment.md)

The historical frozen Week 4 environment records:

| Component | Version |
|---|---:|
| R | 4.5.1 |
| edgeR | 4.6.3 |
| caret | 7.0.1 |
| randomForest | 4.7.1.2 |
| pROC | 1.19.0.1 |

### Restore

From the repository root, using the R version recorded in `renv.lock`:

~~~bash
Rscript -e 'renv::restore(prompt = FALSE)'
~~~

### Execute the strict workflows

~~~bash
Rscript scripts/run_week3.R
Rscript scripts/run_week4.R
~~~

A successful Week 4 reproduction terminates with:

~~~text
WEEK 4 REPRODUCTION: PASS
~~~

Week 1 and Week 2 historical foundation workflows remain available as:

~~~bash
Rscript --vanilla scripts/run_week1.R
Rscript --vanilla scripts/run_week2.R
~~~

---

## Quality assurance

Three formal remediation gates were introduced after the portfolio QA review
without changing the previously frozen Week 3 or Week 4 analytical results.

| QA gate | Result |
|---|---:|
| U08 / P1-04 preprocessing contract | 35 PASS / 0 FAIL |
| U11 reproducible environment | 24 PASS / 0 FAIL |
| U07 / P1-10 confounding assessment | 15 PASS / 0 FAIL |
| **Consolidated** | **74 PASS / 0 FAIL** |

The remediation did not:

- reselect genes;
- change preprocessing retrospectively;
- alter nested cross-validation;
- change the frozen 12-gene set;
- modify `samplewise_rank`;
- refit the Week 4 model;
- retune `mtry`;
- change the classification threshold;
- recompute blind predictions using response labels.

The completed remediation is preserved at tag:

`qa-remediation-w3-w4-complete`

---

## Repository structure

~~~text
.
├── data_raw/          # Public source data
├── data_processed/    # Validated and derived analytical objects
├── scripts/           # Analysis, reproduction and QA scripts
├── results/           # Analytical outputs and QA tables
├── figures/           # Generated figures
├── docs/              # Methodological and audit documentation
├── DATA.md            # Data provenance and analytical contract
├── DESCRIPTION        # R dependency contract
├── renv.lock          # Reproducible environment snapshot
└── README.md
~~~

---

## Key documentation

- [`DATA.md`](DATA.md) — data provenance and endpoint reconstruction
- [`docs/preprocessing_contract.md`](docs/preprocessing_contract.md) — preprocessing boundaries
- [`docs/reproducible_environment.md`](docs/reproducible_environment.md) — reproducible environment
- [`docs/confounding_assessment.md`](docs/confounding_assessment.md) — confounding and cohort heterogeneity
- [`docs/week3_historical_vs_strict_stability_audit.md`](docs/week3_historical_vs_strict_stability_audit.md) — historical versus strict stability
- [`docs/week4_locked_external_validation.md`](docs/week4_locked_external_validation.md) — locked external validation
- [`docs/historical_thesis_results.md`](docs/historical_thesis_results.md) — historical thesis outputs retained for provenance

---

## Historical thesis results

The repository originally accompanied the MSc thesis:

**Validació i optimització d’una signatura molecular predictiva de resposta
a immunoteràpia en melanoma mitjançant un pipeline bioinformàtic reproduïble**

Historical results remain available for provenance but are separated from the
current strict methodological conclusions.

See
[`docs/historical_thesis_results.md`](docs/historical_thesis_results.md).

Original academic repository state:

`tfm-original`

---

## Data sources

Public transcriptomic datasets:

- GSE160638
- GSE91061
- GSE78220

The authoritative GSE160638 treatment-response endpoint used in the rebuild is
derived from supplementary clinical Table S2A associated with:

*MYC Induces Immunotherapy and IFNγ Resistance Through Downregulation of JAK2*

DOI: `10.1158/2326-6066.CIR-22-0184`

Full provenance is available in [`DATA.md`](DATA.md).

---

## Scope and limitations

This repository is a methodological bioinformatics and machine-learning
portfolio project.

It demonstrates:

- reproducible cohort construction;
- transcriptomic preprocessing;
- nested model validation;
- feature-stability analysis;
- frozen external validation;
- reproducible computational environments;
- explicit methodological QA.

It does **not** establish a clinically validated biomarker or molecular
signature.

External validation shows limited and cohort-dependent transportability of the
current frozen deployment procedure.

---

## Author

**Joan Carreras Díaz**

MSc Bioinformatics & Biostatistics

BSc Biomedical Engineering
