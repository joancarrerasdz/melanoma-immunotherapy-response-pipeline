# Data dictionary

## Purpose

This document is the canonical data dictionary for the reproducible melanoma
anti-PD-1 response pipeline.

It documents the identifiers, clinical variables, response codings, units,
expression representations and harmonized metadata used by the analytical
workflow.

The dictionary complements `DATA.md`, which describes data provenance,
cohort construction and the analytical contract.

The primary modelling cohort is GSE160638. GSE91061 and GSE78220 are
independent external-validation cohorts.

---

## Identifier hierarchy

| Field | Meaning | Format / coding | Unit |
|---|---|---|---|
| `dataset` | GEO dataset/cohort identifier | `GSE160638`, `GSE91061`, `GSE78220` | N/A |
| `sample_id` | GEO sample accession | `GSM<number>` | N/A |
| `sample_name` | Source or study-specific sample label | e.g. `PD1_1`, `Pt1_Pre_AD101148-6`, `Pt1` | N/A |
| `patient_id_raw` | Patient identifier as represented in the source clinical table | Source string | N/A |
| `patient_id` | Normalized patient identifier used for GSE160638 alignment | `PD1_<number>` | N/A |

For the primary GSE160638 cohort, `patient_id` is the canonical identifier
used to align the clinical table, GEO metadata and expression-count columns.

`sample_id` remains the GEO accession and is not interchangeable with
`patient_id`.

---

## Primary GSE160638 analytical metadata

Canonical file:

`data_processed/metadata_GSE160638.csv`

| Field | Meaning | Coding / values | Unit |
|---|---|---|---|
| `sample_id` | GEO accession for the RNA-seq sample | `GSM<number>` | N/A |
| `sample_name` | Study sample identifier | `PD1_<number>` | N/A |
| `patient_id` | Normalized patient identifier | `PD1_<number>` | N/A |
| `title` | GEO sample title | Free text | N/A |
| `source_name` | Biological source reported by GEO | `melanoma` in the primary cohort | N/A |
| `characteristics` | GEO free-text disease, tissue and treatment metadata | Free text; not used to infer the response endpoint | N/A |
| `response_recist` | Source treatment-response category retained for the analytical cohort | `CR`, `PR`, `SD`, `PD` | N/A |
| `response` | Binary modelling endpoint derived from `response_recist` | `Responder`, `NonResponder` | N/A |

### Primary response coding

The authoritative response endpoint is reconstructed from the official
supplementary Table S2A.

Binary coding is:

- `CR` -> `Responder`
- `PR` -> `Responder`
- `SD` -> `NonResponder`
- `PD` -> `NonResponder`

GEO free-text metadata are not used to infer the response endpoint.

---

## GSE160638 source clinical variables

Canonical source-derived file:

`data_processed/GSE160638_TableS2A_PD1_clinical.csv`

| Field | Meaning | Coding / scale | Unit |
|---|---|---|---|
| `patient_id_raw` | Original patient identifier read from Table S2A | Source identifier | N/A |
| `patient_id` | Normalized patient identifier | `PD1_<number>` | N/A |
| `TIL ID` | Source TIL identifier when present | Source identifier or missing | N/A |
| `Cohort` | Treatment cohort in the supplementary table | `PD1` for the primary cohort | N/A |
| `Gender` | Patient sex/gender field as provided by source | `M`, `F` | N/A |
| `Age at treatment` | Age when treatment was administered | Numeric | years |
| `Biopsy site` | Anatomical biopsy location | Source categorical/free text | N/A |
| `Biopsy type` | Source biopsy/tissue category | Source categorical/free text | N/A |
| `BRAF mutation` | Source BRAF mutation indicator | `1` present, `0` absent, missing when unavailable | N/A |
| `LDH high` | Source indicator of elevated LDH | `1` high, `0` not high, missing when unavailable | N/A |
| `Response` | Original clinical response category | `CR`, `PR`, `SD`, `PD` | N/A |
| `OR` | Objective-response indicator consistent with CR/PR response | `1` response, `0` no objective response | N/A |
| `PFS T` | Source progression-free-survival time field | Numeric source value; preserved without analytical reinterpretation | Source time scale |
| `PFS E` | Source progression-free-survival event field | Binary source event indicator | N/A |
| `OS T` | Source overall-survival time field | Numeric source value; preserved without analytical reinterpretation | Source time scale |
| `OS E` | Source overall-survival event field | Binary source event indicator | N/A |
| `Prior MAPKi` | Prior MAPK-inhibitor treatment indicator | `1` yes, `0` no | N/A |
| `Prior Ipi` | Prior ipilimumab treatment indicator | `1` yes, `0` no | N/A |
| `Prior treatments` | Other treatments administered before the indexed therapy | Source free text or missing | N/A |
| `Post treatments` | Treatments administered after the indexed therapy | Source free text or missing | N/A |
| `Immune cluster` | Source immune-state grouping | `Low`, `Med`, `High` | N/A |
| `ImmuneScore` | Continuous immune score supplied by the source table | Numeric score | source-defined score |
| `response_recist` | Preserved response category used by the rebuild | `CR`, `PR`, `SD`, `PD` | N/A |
| `response` | Derived binary treatment-response endpoint | `Responder`, `NonResponder` | N/A |

`PFS T`, `PFS E`, `OS T`, `OS E`, `Immune cluster` and `ImmuneScore` are
preserved source variables. They are not used to redefine the primary
binary classification endpoint.

Where a measurement unit is not encoded explicitly in the repository-level
source field, this dictionary records that limitation rather than assigning
an undocumented unit.

---

## Missing clinical values

Missing values in source clinical metadata represent unavailable
source information unless explicitly documented otherwise.

The Week 1 rebuild does not impute clinical missingness.

Important examples include:

- `BRAF mutation`
- `LDH high`
- `Prior treatments`
- `Post treatments`
- `TIL ID`

Missing source values must not be silently converted to biological
absence or clinical negativity.

---

## GSE160638 clinical-to-expression join audit

Canonical file:

`data_processed/GSE160638_clinical_join_audit.csv`

| Field | Meaning | Coding / scale | Unit |
|---|---|---|---|
| `patient_id` | Canonical PD1 patient/sample identifier | `PD1_<number>` | N/A |
| `in_counts` | Whether the identifier exists in the raw-count matrix | `TRUE`, `FALSE` | N/A |
| `in_geo_metadata` | Whether the identifier exists in reconstructed GEO metadata | `TRUE`, `FALSE` | N/A |
| `in_table_s2a` | Whether the identifier exists in clinical Table S2A | `TRUE`, `FALSE` | N/A |
| `count_order` | Position of the sample in the aligned count matrix | Positive integer | index |
| `metadata_order` | Position of the sample in aligned metadata | Positive integer | index |
| `clinical_order` | Position of the sample in the clinical table after alignment | Positive integer | index |
| `status` | Result of the cross-source identifier audit | `matched` when all required sources agree | N/A |

The validated primary cohort contains 36 matched PD1 samples.

---

## External-validation harmonized metadata

Canonical files include:

- `data_processed/metadata_GSE91061_harmonized_kept.csv`
- `data_processed/metadata_GSE78220_harmonized_kept.csv`

Both use the following harmonized schema.

| Field | Meaning | Coding / scale | Unit |
|---|---|---|---|
| `dataset` | External GEO cohort | `GSE91061` or `GSE78220` | N/A |
| `sample_id` | GEO sample accession | `GSM<number>` | N/A |
| `sample_name` | Study-specific patient/sample label | Source string | N/A |
| `response_raw` | Response label before project harmonization | Cohort-specific source category | N/A |
| `response` | Harmonized binary external-validation endpoint | `Responder`, `NonResponder` | N/A |
| `treatment_time` | Timing of biopsy relative to treatment | `PreTreatment` for the locked evaluation samples | N/A |
| `keep` | Inclusion flag for the locked external evaluation | `TRUE`, `FALSE` | N/A |
| `reason_excluded` | Explicit reason for exclusion when `keep = FALSE` | Text or missing | N/A |

### GSE91061 response coding

Observed locked pre-treatment metadata use:

- `PRCR` -> `Responder`
- `SD` -> `NonResponder`
- `PD` -> `NonResponder`

### GSE78220 response coding

Observed locked pre-treatment metadata use:

- `Complete Response` -> `Responder`
- `Partial Response` -> `Responder`
- `Progressive Disease` -> `NonResponder`

External response labels are used only during the locked Week 4 evaluation
after the model specification and blind predictions have been frozen.

---

## Expression representations and units

### Raw expression counts

Canonical primary object:

`data_processed/GSE160638_raw_counts_PD1_aligned.rds`

Representation:

- rows: genes/features;
- columns: 36 aligned PD1 samples;
- values: non-negative integer-like RNA-seq counts;
- scale: raw count scale;
- missing values: none in the validated object.

Raw counts are not interpreted as normalized expression values.

### Exploratory logCPM

Canonical exploratory object:

`data_processed/GSE160638_logCPM_exploratory.rds`

Representation:

- filtered expression features;
- TMM-normalized exploratory expression;
- logCPM scale;
- continuous transformed values.

This object is for exploratory analyses and is not a globally preprocessed
input to predictive modelling.

### Strict predictive preprocessing

Predictive preprocessing is estimated within training partitions during
strict model validation.

No supervised feature selection is fitted globally before cross-validation.

### Locked external representation

The Week 4 deployment representation is `samplewise_rank`.

For each external sample, the frozen deployment genes are ranked within
that sample.

The resulting values are:

- ordinal;
- sample-relative;
- unitless.

This representation is fixed before external response labels are opened.

---

## Biological interpretation boundary

Clinical and expression fields describe observational patient/sample data.

The project does not interpret:

- a response-associated gene as causally determining treatment response;
- a selected gene as a clinically validated biomarker;
- missing clinical metadata as evidence of biological absence;
- cohort-specific expression distributions as directly exchangeable across
  studies without the locked representation and validation procedure.

The 12-gene strict consensus is a candidate research signature and not a
clinically validated biomarker panel.

---

## Canonical metadata sources

Primary cohort:

- `data_processed/metadata_GSE160638.csv`
- `data_processed/GSE160638_TableS2A_PD1_clinical.csv`
- `data_processed/GSE160638_clinical_join_audit.csv`
- `data_processed/GSE160638_raw_counts_PD1_aligned.rds`

External evaluation:

- `data_processed/metadata_GSE91061_harmonized_kept.csv`
- `data_processed/metadata_GSE78220_harmonized_kept.csv`

Project-level provenance and analytical context:

- `DATA.md`
- `data_raw/clinical/README.md`
- `docs/preprocessing_contract.md`
- `docs/week4_locked_external_validation.md`

---

## QA interpretation

This dictionary documents:

- variable meaning;
- measurement scale or unit where applicable;
- explicit acknowledgement when a source unit is not encoded;
- identifier roles;
- clinical endpoint codings;
- external-response harmonization;
- missing-value semantics;
- expression-data scales;
- biological interpretation boundaries.

It is documentation only.

Creating this dictionary does not alter raw data, processed data, cohort
membership, response labels, preprocessing, model fitting, predictions,
metrics or the released `v1.0.0` analytical baseline.
