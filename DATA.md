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

## Week 2 exploratory preprocessing and analysis

Week 2 builds directly on the validated Week 1 analytical foundation.

The Week 2 exploratory analyses use:

- 36 validated PD1 samples;
- 17,002 raw expression features;
- 22 Responders;
- 14 NonResponders;
- the sample ordering established and validated during Week 1.

No samples are added, removed or reordered during Week 2.

### Exploratory gene filtering

Expression filtering is performed with `edgeR::filterByExpr`.

Starting from 17,002 genes:

- 15,097 genes are retained;
- 1,905 genes are removed.

The filtering performed here is intended for exploratory analysis only.

It must not be interpreted as a globally fitted feature-selection step for predictive modelling.

### TMM normalization and exploratory logCPM

The filtered exploratory expression data are normalized using edgeR TMM normalization.

The resulting exploratory expression matrix contains:

- 15,097 genes;
- 36 samples;
- no missing values;
- only finite expression values.

The exploratory logCPM matrix is stored as:

`data_processed/GSE160638_logCPM_exploratory.rds`

This matrix is intended for descriptive and exploratory analyses such as PCA and sample-correlation assessment.

It is not a valid globally preprocessed input for predictive modelling.

Predictive preprocessing must be fitted independently within the training partitions of the modelling workflow.

### Exploratory PCA

Principal component analysis is performed on the exploratory filtered, TMM-normalized logCPM matrix.

The PCA includes all 36 validated PD1 samples.

The first two principal components explain:

- PC1: 17.24% of total variance;
- PC2: 10.14% of total variance.

The corresponding exploratory outputs are:

- `results/week2_pca_scores_exploratory.csv`
- `results/week2_pca_variance_exploratory.csv`
- `figures/week2_pca_exploratory.png`

### Sample-correlation analysis

Pairwise Pearson correlations are calculated using the exploratory logCPM matrix.

The resulting sample-correlation matrix has dimensions:

`36 x 36`

The corresponding outputs are:

- `results/week2_sample_correlation_exploratory.csv`
- `figures/week2_sample_correlation_exploratory.png`

This analysis is descriptive and is used to inspect global sample-expression similarity.

### Exploratory differential expression

Differential-expression analysis is performed with edgeR using the validated Week 1 response endpoint.

The model contrast is explicitly defined as:

`Responder - NonResponder`

Therefore:

- positive logFC indicates higher expression in Responders;
- negative logFC indicates higher expression in NonResponders.

The analysis uses 15,097 genes after `filterByExpr`.

The exploratory differential-expression results contain:

- 2,044 genes with FDR < 0.05;
- 3,134 genes with FDR < 0.10;
- 6,603 genes with positive logFC;
- 8,494 genes with negative logFC.

The observed logFC range is approximately:

`-7.290606 to 5.186878`

The corresponding outputs are:

- `results/week2_differential_expression_exploratory.csv`
- `results/week2_differential_expression_summary.csv`

These differential-expression results are exploratory only.

They must not be used to globally preselect genes before predictive model validation.

Any predictive feature-selection procedure must be fitted independently within training folds.

### Week 2 audit outputs

#### Preprocessing

- `results/week2_gene_filtering_exploratory.csv`
- `results/week2_tmm_normalization_exploratory.csv`
- `results/week2_preprocessing_summary.csv`
- `data_processed/GSE160638_logCPM_exploratory.rds`

#### Exploratory analysis

- `results/week2_pca_scores_exploratory.csv`
- `results/week2_pca_variance_exploratory.csv`
- `results/week2_sample_correlation_exploratory.csv`
- `results/week2_exploratory_analysis_summary.csv`
- `figures/week2_pca_exploratory.png`
- `figures/week2_sample_correlation_exploratory.png`

#### Exploratory differential expression

- `results/week2_differential_expression_exploratory.csv`
- `results/week2_differential_expression_summary.csv`

### Week 2 methodological boundary

Week 2 establishes an exploratory preprocessing and expression-analysis layer.

The filtered/TMM/logCPM matrix, PCA results, sample correlations and differential-expression results are descriptive or exploratory outputs.

They are not globally fitted preprocessing or feature-selection inputs for predictive modelling.

To prevent information leakage, predictive preprocessing and feature selection must be estimated independently inside the training partitions of subsequent model-validation procedures.


## Week 3 strict nested-CV outputs

Week 3 uses the validated Week 1 clinical and raw-count foundation directly.

### Validated modelling inputs

- `data_processed/GSE160638_raw_counts_PD1_aligned.rds`
- `data_processed/GSE160638_TableS2A_PD1_clinical.csv`

Validated dimensions:

- 17,002 genes;
- 36 samples;
- 22 Responders;
- 14 NonResponders.

The Week 2 exploratory filtered/TMM/logCPM matrix is **not** used as a
predictive modelling input.

Predictive preprocessing and feature selection are estimated independently
inside training partitions.

### Strict nested cross-validation outputs

- `results/week3_inner_cv_metrics_strict.csv`
- `results/week3_inner_tuning_summary_strict.csv`
- `results/week3_outer_predictions_strict.csv`
- `results/week3_outer_fold_metrics_strict.csv`
- `results/week3_outer_selected_genes_strict.csv`
- `results/week3_gene_selection_frequency_strict.csv`
- `results/week3_nested_cv_summary_strict.csv`
- `figures/week3_nested_cv_roc_strict.png`

The nested-CV design uses 5 outer folds and 3 inner folds.

Each sample receives exactly one outer out-of-fold prediction.

Validation performance is calculated exclusively from the pooled outer
out-of-fold predictions.

Validated performance:

- Accuracy: 0.8333;
- Sensitivity: 0.9091;
- Specificity: 0.7143;
- Balanced Accuracy: 0.8117;
- AUC: 0.8896.

### Strict gene-stability outputs

- `results/week3_gene_stability_tiers_strict.csv`
- `results/week3_consensus_repeated_2of5_strict.csv`
- `results/week3_consensus_majority_3of5_strict.csv`
- `results/week3_consensus_high_4of5_strict.csv`
- `results/week3_consensus_core_5of5_strict.csv`
- `figures/week3_gene_selection_stability_strict.png`

Observed stability:

- at least 1/5 outer folds: 312 genes;
- at least 2/5 outer folds: 100 genes;
- at least 3/5 outer folds: 51 genes;
- at least 4/5 outer folds: 25 genes;
- 5/5 outer folds: 12 genes.

The 12 genes selected in all five outer folds are core consensus candidates,
not an independently validated final signature.

### Historical-versus-strict audit outputs

- `results/week3_historical_vs_strict_5of5.csv`
- `results/week3_historical_vs_strict_5of5_summary.csv`
- `results/week3_historical_13_genes_in_strict_pipeline.csv`
- `results/week3_historical_13_genes_in_strict_pipeline_summary.csv`
- `docs/week3_historical_vs_strict_stability_audit.md`

Historical audit:

- historical 5/5 genes: 13;
- strict 5/5 genes: 12;
- shared 5/5 genes: 0;
- historical genes present in validated input: 13/13;
- historical genes selected in at least one strict outer fold: 0/13.

The zero overlap documents failure to reproduce the historical stability set
under the strict nested-CV procedure. It does not demonstrate biological
irrelevance of the historical genes.

### Week 3 methodological boundary

The strict nested-CV metrics are validation estimates.

The consensus-gene outputs are post-validation stability analyses.

Any model subsequently fitted to a consensus signature using the complete
36-sample cohort is exploratory and must not be reported as independent
validation on those same samples.


## Reproduction

The reproducible rebuild can be executed by week from the repository root.

Week 1:

```bash
Rscript --vanilla scripts/run_week1.R
```

Week 2:

```bash
Rscript --vanilla scripts/run_week2.R
```

Week 3:

```bash
Rscript --vanilla scripts/run_week3.R
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