# TFM immunotherapy melanoma pipeline

> [!WARNING]
> **Results under methodological revalidation**
>
> The metrics, gene signatures, and biological interpretations currently included in this repository are historical outputs from the original master's thesis and are being methodologically revalidated. A clinical-endpoint labelling issue was identified in GSE160638; these outputs must not be interpreted as clinically validated evidence. The original academic state is preserved at tag `tfm-original`, and the reconstruction is being developed outside `main`.

## Week 1 reproducible rebuild

The current methodological rebuild is being developed outside `main`, starting from the frozen academic baseline preserved at tag `tfm-original`.

Week 1 establishes the validated clinical and raw-expression foundation for the primary GSE160638 anti–PD-1 cohort before any predictive modelling is performed.

The validated cohort contains:

- 36 unique PD1 samples;
- 22 Responders;
- 14 NonResponders;
- no TIL-ACT/TIL samples in the primary analytical cohort.

The official treatment-response endpoint is reconstructed from supplementary Table S2A and explicitly mapped from `CR/PR/SD/PD` to the binary `Responder/NonResponder` outcome.

### Reproduce Week 1

From the repository root:

```bash
Rscript --vanilla scripts/run_week1.R
```

If the required Week 1 R packages are not installed:

```bash
Rscript --vanilla scripts/install_week1_dependencies.R
```

The Week 1 runner performs:

1. clinical and expression-data reconstruction and alignment;
2. automated integrity tests;
3. initial descriptive quality control;
4. validation of the expected Week 1 outputs.

Detailed data provenance, cohort definitions, transformations and outputs are documented in [`DATA.md`](DATA.md).

The authoritative clinical-source record is documented in [`data_raw/clinical/README.md`](data_raw/clinical/README.md).

QC findings and explicitly documented anomalies are recorded in [`docs/day4_qc_anomalies.md`](docs/day4_qc_anomalies.md).

## Week 2 reproducible preprocessing and exploratory analysis

Week 2 builds on the validated Week 1 clinical and raw-expression foundation for the primary GSE160638 anti–PD-1 cohort.

The validated input remains:

- 36 PD1 samples;
- 22 Responders;
- 14 NonResponders;
- 17,002 raw expression features;
- identical clinical and expression sample order.

Week 2 adds a reproducible exploratory preprocessing and analysis layer.

The exploratory preprocessing performs:

1. expression filtering with `edgeR::filterByExpr`;
2. TMM library-size normalization;
3. logCPM transformation for exploratory analysis.

The exploratory filtering retains 15,097 genes and removes 1,905 genes.

Exploratory analyses include:

- PCA of the filtered/TMM-normalized logCPM matrix;
- sample-to-sample expression correlation;
- exploratory differential-expression analysis between Responders and NonResponders.

For the exploratory PCA:

- PC1 explains 17.24% of the variance;
- PC2 explains 10.14% of the variance.

The exploratory differential-expression contrast is:

`Responder - NonResponder`

and identifies:

- 2,044 genes with FDR < 0.05;
- 3,134 genes with FDR < 0.10.

### Reproduce Week 2

From the repository root:

```bash
Rscript --vanilla scripts/run_week2.R
```

The Week 2 runner executes:

1. exploratory preprocessing;
2. exploratory PCA and sample-correlation analysis;
3. exploratory differential-expression analysis;
4. validation of the expected Week 2 outputs.

### Methodological boundary

The Week 2 filtered/TMM/logCPM matrix and differential-expression results are **exploratory only**.

They must not be used as globally preprocessed or globally preselected inputs for predictive modelling.

Predictive preprocessing and feature selection must be fitted independently inside the training partitions of the modelling workflow.



## Week 3 — Strict nested cross-validation and gene-stability rebuild

Week 3 reconstructs the predictive modelling layer using a strict nested
cross-validation architecture designed to prevent information leakage.

The validated Week 1 cohort is preserved:

- 36 anti-PD-1 samples;
- 17,002 input genes;
- 22 Responders;
- 14 NonResponders.

### Strict nested-CV design

The rebuilt validation procedure uses:

- 5 outer folds for unbiased out-of-fold evaluation;
- 3 inner folds for hyperparameter tuning;
- 100 genes selected independently inside each training partition;
- training-only expression filtering and feature selection;
- training-only hyperparameter tuning;
- outer-test samples excluded from preprocessing, feature selection,
  tuning and model fitting.

Each sample receives exactly one outer out-of-fold prediction.

Global validation performance is calculated exclusively from the pooled
outer out-of-fold predictions.

### Strict nested-CV performance

The validated pooled outer-fold performance is:

- Accuracy: 0.8333;
- Sensitivity: 0.9091;
- Specificity: 0.7143;
- Balanced Accuracy: 0.8117;
- AUC: 0.8896.

These values represent the validation performance of the strict nested-CV
procedure on the 36-sample GSE160638 cohort.

### Gene-selection stability

Feature-selection stability across the five outer training partitions is:

- selected in at least 1/5 folds: 312 genes;
- selected in at least 2/5 folds: 100 genes;
- selected in at least 3/5 folds: 51 genes;
- selected in at least 4/5 folds: 25 genes;
- selected in 5/5 folds: 12 genes.

The 12 genes selected in all five outer folds are treated as
**core consensus candidates**.

They do not constitute an independently validated final molecular signature.

### Historical-versus-strict audit

The historical workflow contained 13 genes reported as stable in 5/5 folds.

The reproducibility audit shows:

- historical 5/5 genes: 13;
- strict 5/5 genes: 12;
- shared 5/5 genes: 0;
- historical genes present in the validated 17,002-gene input: 13/13;
- historical genes selected in at least one strict outer fold: 0/13.

Therefore, the historical stable set is not reproduced by the strict
nested-CV selection procedure.

This result does not demonstrate biological irrelevance of the historical
genes. It documents a lack of reproducibility of their previously reported
selection stability under the stricter modelling architecture.

### Reproduce Week 3

From the repository root:

```bash
Rscript --vanilla scripts/run_week3.R
```

The Week 3 runner executes:

1. strict nested cross-validation;
2. strict gene-selection stability analysis;
3. historical-versus-strict stability audit;
4. validation of the expected Week 3 outputs and metrics.

### Week 3 methodological boundary

Validation performance and post-validation consensus-gene analyses are kept
strictly separate.

Only pooled outer out-of-fold predictions are used to calculate validation
performance.

Consensus sets derived after nested cross-validation are stability-derived
candidate signatures and must not be reported as independently validated
final signatures.

Historical full-dataset Random Forest importance analyses are exploratory and
are not used to define the strict consensus candidates.


This repository contains the R code, input public datasets, figures and result tables for the master's thesis project:

**Validació i optimització d’una signatura molecular predictiva de resposta a immunoteràpia en melanoma mitjançant un pipeline bioinformàtic reproduïble**

Author: Joan Carreras Díaz  
Programme: Màster Universitari en Bioinformàtica i Bioestadística  
Area: Desenvolupament de Programari i Aplicacions

## Project summary

The objective of this project is to develop and validate a reproducible bioinformatics pipeline for predicting response to anti-PD-1 immunotherapy in melanoma using transcriptomic data.

The final approach uses:

- Training cohort: **GSE160638**
- External validation cohorts: **GSE91061** and **GSE78220**
- Clinical harmonisation into `Responder` and `NonResponder`
- Pre-treatment sample restriction when applicable
- Gene filtering and TMM/logCPM transformation with `edgeR`
- Differential-expression-based feature selection inside cross-validation folds
- Random Forest modelling
- Internal validation by cross-validation without data leakage
- External validation on independent cohorts
- Gene stability and biological interpretation analyses

## Repository structure

```text
.
├── scripts/          # Ordered R scripts for the full pipeline
├── data_raw/         # Public input files used by the scripts
├── data_processed/   # Small metadata tables; large RDS caches are regenerated by scripts
├── results/          # Final CSV result tables
├── figures/          # Final figures generated by the pipeline
├── docs/             # Uploaded thesis PDF copy / delivery documents
├── script_manifest.csv
└── README.md
```

## Script execution order

Run the scripts from the project root or from inside the `scripts/` folder.

```r
source("scripts/01_load_data.R")
source("scripts/02_clinical_harmonization.R")
source("scripts/03_preprocessing_train.R")
source("scripts/04_exploratory_analysis.R")
source("scripts/05_differential_expression_exploratory.R")
source("scripts/06_nested_cv_final_100genes.R")
source("scripts/07_signature_size_comparison.R")
source("scripts/08_gene_stability_analysis.R")
source("scripts/09_prepare_external_validation.R")
source("scripts/10_define_common_genes_across_datasets.R")
source("scripts/11_biological_interpretation_final_100.R")
source("scripts/12_external_validation_final_100.R")
source("scripts/13_gene_stability_across_datasets.R")
```

## Main final results

### Internal validation, final 100-gene signature

From `results/rf_cv_metrics_summary_final_100_no_leakage.csv`:

| Metric | Value |
|---|---:|
| Accuracy | 0.753425 |
| Kappa | 0.506006 |
| Sensitivity | 0.810811 |
| Specificity | 0.694444 |
| Balanced Accuracy | 0.752628 |
| AUC | 0.827703 |
| Mean best mtry | 5.200000 |
| Unique genes selected | 292 |
| Genes selected in ≥4 folds | 29 |
| Genes selected in 5/5 folds | 13 |

### Signature-size comparison

From `results_reference/signature_size_comparison_metrics.csv` / final comparison output:

| Number of genes | AUC | Accuracy | Balanced Accuracy |
|---:|---:|---:|---:|
| 50 | 0.800300 | 0.712329 | 0.712087 |
| 100 | 0.832958 | 0.767123 | 0.766517 |
| 250 | 0.831456 | 0.780822 | 0.780405 |
| 500 | 0.827703 | 0.780822 | 0.780030 |

### External validation, final 100-gene model

From `results/external_validation_metrics_by_dataset_final_100.csv`:

| Dataset | AUC | Accuracy | Sensitivity | Specificity | Balanced Accuracy |
|---|---:|---:|---:|---:|---:|
| GSE91061 | 0.670513 | 0.612245 | 0.600000 | 0.615385 | 0.607692 |
| GSE78220 | 0.611111 | 0.444444 | 0.000000 | 1.000000 | 0.500000 |

### Gene-stability distribution, final 100-gene signature

From `results/gene_stability_summary_final_100.csv`:

| Folds selected | Number of genes |
|---:|---:|
| 1 | 180 |
| 2 | 58 |
| 3 | 25 |
| 4 | 16 |
| 5 | 13 |

## Notes on reproducibility

The pipeline avoids data leakage by performing feature selection and model training inside the training partitions. External datasets are used only for final model evaluation.

Large intermediate RDS files are not required in the repository because they can be regenerated by the ordered scripts. Final result tables and figures are included for traceability.

## Data sources

All datasets used in this project are publicly available from Gene Expression Omnibus (GEO):

- GSE160638
- GSE91061
- GSE78220

## Software requirements

The analysis was developed in R. Main packages include:

- `GEOquery`
- `dplyr`
- `readr`
- `readxl`
- `stringr`
- `edgeR`
- `randomForest`
- `caret`
- `pROC`
- `ggplot2`
- `pheatmap`
- `clusterProfiler`
- `org.Hs.eg.db`
- `AnnotationDbi`
- `ReactomePA` (optional)

## Important delivery note

The thesis PDF included in `docs/` is the uploaded copy. If the final written report is updated, make sure that its numerical results match the values listed above and the CSV files in `results/`.
