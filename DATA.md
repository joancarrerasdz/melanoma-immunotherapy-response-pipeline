# Data documentation

This document describes the data foundation used for the Week 1 rebuild of the GSE160638 melanoma anti–PD-1 cohort.

The purpose of Week 1 is to establish a reproducible, auditable clinical and expression-data foundation before any predictive modelling is performed.

## Primary dataset

- GEO accession: `GSE160638`
- Disease context: melanoma
- Treatment context: anti–PD-1 immunotherapy
- Primary analytical cohort: 36 PD1 samples
- Expression input: raw RNA-seq count matrix
- Clinical endpoint source: official supplementary Table S2, worksheet `S2A`

The GEO record provides the transcriptomic data and sample metadata. The authoritative patient-level treatment-response endpoint is obtained from the supplementary clinical table rather than inferred from GEO free-text metadata.

## Clinical source provenance

The official clinical source is:

- Publication: *MYC Induces Immunotherapy and IFNγ Resistance Through Downregulation of JAK2*
- Journal: *Cancer Immunology Research*
- Year: 2023
- DOI: `10.1158/2326-6066.CIR-22-0184`
- Supplementary source: AACR Figshare
- Local file: `data_raw/clinical/cir-22-0184_table_s2_suppst2.xlsx`
- Worksheet used: `S2A`

Additional source and reuse information is recorded in:

`data_raw/clinical/README.md`

## Population

Week 1 restricts the analysis to the primary anti–PD-1 cohort.

After clinical reconstruction and identifier validation:

- 36 unique PD1 patients/samples are retained.
- 22 are classified as Responders.
- 14 are classified as NonResponders.
- No `TIL_*` samples are included in the primary cohort.

TIL-ACT/TIL samples present in the original expression dataset are therefore excluded from the primary anti–PD-1 analytical cohort.

No sample is removed on the basis of the initial descriptive QC performed during Week 1.

## Response endpoint

The original RECIST-like response categories in Table S2A are:

- `CR`
- `PR`
- `SD`
- `PD`

The binary endpoint used in the rebuild is derived explicitly as:

- `CR` or `PR` → `Responder`
- `SD` or `PD` → `NonResponder`

The resulting distribution is:

- Responders: 22
- NonResponders: 14

The original response is retained as `response_recist`, while the binary endpoint is stored as `response`.

No response label is inferred from uncontrolled GEO free-text fields.

## Identifier processing and alignment

Patient identifiers are normalized to the form:

`PD1_<number>`

The Week 1 pipeline verifies that:

- exactly 36 unique clinical identifiers are present;
- no duplicated clinical identifiers are present;
- no duplicated count-matrix identifiers are present;
- all 36 PD1 clinical samples are present in the count matrix;
- clinical, GEO metadata and expression-count identifiers correspond completely;
- expression columns and clinical rows are stored in identical sample order.

A dedicated join-audit table records the correspondence between the three sources.

## Expression data

The aligned expression object contains:

- 17,002 features
- 36 PD1 samples

The Week 1 QC verifies that the raw count matrix contains:

- no missing values;
- no negative values;
- no non-integer-like values.

No normalization, gene filtering or modelling transformation is applied during Week 1. The aligned object remains a raw-count matrix for downstream preprocessing.

## Clinical transformations

Week 1 performs only the transformations required to construct and audit the analytical cohort:

1. Read the official Table S2A clinical source.
2. Normalize PD1 patient identifiers.
3. Restrict the clinical data to the primary anti–PD-1 cohort.
4. Preserve the official `CR/PR/SD/PD` endpoint.
5. Derive the binary `Responder/NonResponder` endpoint.
6. Match the clinical cohort to GEO metadata.
7. Restrict the expression matrix to the 36 PD1 samples.
8. Align clinical rows and count-matrix columns.
9. Generate an explicit join audit.
10. Perform descriptive QC.

Week 1 does not:

- impute clinical missing values;
- silently harmonize inconsistent clinical categories;
- remove QC-flagged samples;
- perform statistical association testing;
- perform feature selection;
- train predictive models.

## Initial QC observations

The Week 1 descriptive QC identified one sample for monitoring:

`PD1_9`

It was flagged by the 1.5 × IQR descriptive rule because it has:

- 16,069 detected genes;
- 5.488% zero counts;
- no library-size flag.

This is a descriptive QC flag only. `PD1_9` remains in the analytical cohort.

Clinical missingness observed during Week 1 includes:

- Age at treatment: 0/36
- Gender: 0/36
- BRAF mutation: 6/36
- LDH high: 5/36
- Biopsy site: 0/36
- Biopsy type: 0/36

BRAF and LDH missingness is retained as source-data missingness and is not treated as a pipeline error.

Differences in capitalization, punctuation and abbreviation in biopsy categories are documented rather than silently harmonized.

## Week 1 outputs

### Validated analytical data

- `data_processed/metadata_GSE160638.csv`
- `data_processed/GSE160638_TableS2A_PD1_clinical.csv`
- `data_processed/GSE160638_raw_counts_PD1_aligned.rds`
- `data_processed/GSE160638_clinical_join_audit.csv`

### QC tables

- `results/day4_sample_qc.csv`
- `results/day4_clinical_missingness.csv`
- `results/day4_clinical_by_response.csv`
- `results/day4_age_by_response.csv`

### QC figures

- `figures/day4_response_distribution.png`
- `figures/day4_library_size.png`
- `figures/day4_genes_detected.png`
- `figures/day4_percent_zeros.png`

### QC documentation

- `docs/day4_qc_anomalies.md`

## Reproduction

The complete Week 1 rebuild is executed from the repository root with:

```bash
Rscript --vanilla scripts/run_week1.R
```

Required Week 1 dependencies can be installed with:

```bash
Rscript --vanilla scripts/install_week1_dependencies.R
```

## Environment and dependencies

The environment used to validate the Week 1 rebuild was:

- R 4.5.1
- GEOquery 2.76.0
- dplyr 1.2.0
- readxl 1.4.5

The installation script installs compatible available versions. The versions above document the environment in which the Week 1 gate was validated; they are not an exact dependency lock file.