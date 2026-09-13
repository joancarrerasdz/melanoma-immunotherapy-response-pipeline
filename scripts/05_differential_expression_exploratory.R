#!/usr/bin/env Rscript

# ============================================================================
# Week 2 — Exploratory differential-expression analysis
# GSE160638 primary anti-PD-1 cohort
#
# IMPORTANT:
# This analysis is descriptive / biological exploration only.
# Its results MUST NOT be used as globally preselected predictive features.
# Predictive preprocessing and feature selection must be fitted independently
# inside the training folds.
# ============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
})

# ----------------------------------------------------------------------------
# 1. Paths
# ----------------------------------------------------------------------------

base_dir <- getwd()

data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

counts_file <- file.path(
  data_processed_dir,
  "GSE160638_raw_counts_PD1_aligned.rds"
)

clinical_file <- file.path(
  data_processed_dir,
  "GSE160638_TableS2A_PD1_clinical.csv"
)

required_files <- c(counts_file, clinical_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    paste0(
      "Missing validated Week 1 input files:\n  ",
      paste(missing_files, collapse = "\n  ")
    ),
    call. = FALSE
  )
}

# ----------------------------------------------------------------------------
# 2. Load validated Week 1 foundation
# ----------------------------------------------------------------------------

counts <- readRDS(counts_file)
counts <- as.matrix(counts)

clinical <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_clinical_columns <- c("patient_id", "response")

missing_columns <- setdiff(
  required_clinical_columns,
  colnames(clinical)
)

if (length(missing_columns) > 0L) {
  stop(
    paste(
      "Missing required clinical columns:",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

# ----------------------------------------------------------------------------
# 3. Input integrity checks
# ----------------------------------------------------------------------------

if (!is.numeric(counts)) {
  stop("Count matrix is not numeric.", call. = FALSE)
}

if (nrow(counts) != 17002L) {
  stop(
    paste0(
      "Unexpected raw feature count: ",
      nrow(counts),
      ". Expected 17002."
    ),
    call. = FALSE
  )
}

if (ncol(counts) != 36L) {
  stop(
    paste0(
      "Unexpected sample count: ",
      ncol(counts),
      ". Expected 36."
    ),
    call. = FALSE
  )
}

if (nrow(clinical) != 36L) {
  stop(
    paste0(
      "Unexpected clinical row count: ",
      nrow(clinical),
      ". Expected 36."
    ),
    call. = FALSE
  )
}

if (anyDuplicated(clinical$patient_id)) {
  stop("Duplicated clinical patient identifiers detected.", call. = FALSE)
}

if (anyDuplicated(colnames(counts))) {
  stop("Duplicated count-matrix identifiers detected.", call. = FALSE)
}

if (!identical(colnames(counts), clinical$patient_id)) {
  stop(
    "Count columns and validated clinical patient IDs do not have identical order.",
    call. = FALSE
  )
}

if (anyNA(counts)) {
  stop("Missing values detected in raw counts.", call. = FALSE)
}

if (any(counts < 0)) {
  stop("Negative values detected in raw counts.", call. = FALSE)
}

non_integer_like <- sum(
  abs(counts - round(counts)) > 1e-08,
  na.rm = TRUE
)

if (non_integer_like != 0L) {
  stop(
    paste(
      "Non-integer-like values detected in raw counts:",
      non_integer_like
    ),
    call. = FALSE
  )
}

if (anyNA(clinical$response) || any(clinical$response == "")) {
  stop("Missing response labels detected.", call. = FALSE)
}

if (!all(clinical$response %in% c("NonResponder", "Responder"))) {
  stop("Unexpected response labels detected.", call. = FALSE)
}

# ----------------------------------------------------------------------------
# 4. Response definition
# ----------------------------------------------------------------------------

group <- factor(
  clinical$response,
  levels = c("NonResponder", "Responder")
)

expected_response_counts <- c(
  NonResponder = 14L,
  Responder = 22L
)

observed_response_counts <- table(group)

if (!all(
  observed_response_counts[names(expected_response_counts)] ==
    expected_response_counts
)) {
  stop(
    "Response distribution does not match the validated Week 1 cohort.",
    call. = FALSE
  )
}

# ----------------------------------------------------------------------------
# 5. edgeR gene filtering
# ----------------------------------------------------------------------------

dge <- DGEList(
  counts = counts,
  group = group
)

keep_genes <- filterByExpr(
  dge,
  group = group
)

dge <- dge[
  keep_genes,
  ,
  keep.lib.sizes = FALSE
]

# ----------------------------------------------------------------------------
# 6. TMM normalization
# ----------------------------------------------------------------------------

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

# ----------------------------------------------------------------------------
# 7. Design matrix
#
# group levels:
#   reference = NonResponder
#   coefficient 2 = Responder - NonResponder
# ----------------------------------------------------------------------------

design <- model.matrix(~ group)

if (
  ncol(design) != 2L ||
  colnames(design)[2] != "groupResponder"
) {
  stop(
    "Unexpected design matrix: coefficient 2 is not groupResponder.",
    call. = FALSE
  )
}

# ----------------------------------------------------------------------------
# 8. Quasi-likelihood differential-expression model
# ----------------------------------------------------------------------------

dge <- estimateDisp(
  dge,
  design
)

fit <- glmQLFit(
  dge,
  design
)

qlf <- glmQLFTest(
  fit,
  coef = 2
)

de <- topTags(
  qlf,
  n = Inf,
  sort.by = "PValue"
)$table

de_output <- data.frame(
  gene_id = rownames(de),
  de,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# 9. Summary statistics
# ----------------------------------------------------------------------------

n_fdr_005 <- sum(de_output$FDR < 0.05, na.rm = TRUE)
n_fdr_010 <- sum(de_output$FDR < 0.10, na.rm = TRUE)

n_positive_logfc <- sum(de_output$logFC > 0, na.rm = TRUE)
n_negative_logfc <- sum(de_output$logFC < 0, na.rm = TRUE)

de_summary <- data.frame(
  metric = c(
    "input_genes",
    "input_samples",
    "genes_after_filterByExpr",
    "genes_removed_by_filterByExpr",
    "responders",
    "nonresponders",
    "fdr_lt_0.05",
    "fdr_lt_0.10",
    "positive_logFC",
    "negative_logFC",
    "minimum_logFC",
    "maximum_logFC"
  ),
  value = c(
    nrow(counts),
    ncol(counts),
    nrow(dge),
    nrow(counts) - nrow(dge),
    sum(group == "Responder"),
    sum(group == "NonResponder"),
    n_fdr_005,
    n_fdr_010,
    n_positive_logfc,
    n_negative_logfc,
    min(de_output$logFC, na.rm = TRUE),
    max(de_output$logFC, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# 10. Write exploratory outputs
# ----------------------------------------------------------------------------

write.csv(
  de_output,
  file.path(
    results_dir,
    "week2_differential_expression_exploratory.csv"
  ),
  row.names = FALSE
)

write.csv(
  de_summary,
  file.path(
    results_dir,
    "week2_differential_expression_summary.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 11. Final validation
# ----------------------------------------------------------------------------

if (nrow(de_output) != nrow(dge)) {
  stop(
    "Differential-expression output has unexpected dimensions.",
    call. = FALSE
  )
}

if (anyNA(de_output$PValue) || anyNA(de_output$FDR)) {
  stop(
    "Missing PValue/FDR values detected in DE output.",
    call. = FALSE
  )
}

if (
  any(de_output$FDR < 0) ||
  any(de_output$FDR > 1)
) {
  stop(
    "Invalid FDR values detected.",
    call. = FALSE
  )
}

# ----------------------------------------------------------------------------
# 12. Report
# ----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("WEEK 2 EXPLORATORY DIFFERENTIAL EXPRESSION: PASS\n")
cat("============================================================\n")

cat("Validated input samples:", ncol(counts), "\n")
cat("Validated input genes:", nrow(counts), "\n")
cat("Genes after filterByExpr:", nrow(dge), "\n")
cat("Genes removed:", nrow(counts) - nrow(dge), "\n")

cat("Responders:", sum(group == "Responder"), "\n")
cat("NonResponders:", sum(group == "NonResponder"), "\n")

cat("Contrast: Responder - NonResponder\n")

cat("FDR < 0.05:", n_fdr_005, "\n")
cat("FDR < 0.10:", n_fdr_010, "\n")

cat("Positive logFC:", n_positive_logfc, "\n")
cat("Negative logFC:", n_negative_logfc, "\n")

cat(
  "logFC range:",
  sprintf("%.6f", min(de_output$logFC, na.rm = TRUE)),
  "to",
  sprintf("%.6f", max(de_output$logFC, na.rm = TRUE)),
  "\n"
)

cat("\nIMPORTANT:\n")
cat("This differential-expression analysis is exploratory only.\n")
cat("Positive logFC = higher expression in Responders.\n")
cat("Negative logFC = higher expression in NonResponders.\n")
cat("These genes are not globally preselected for predictive modelling.\n")
cat("Predictive feature selection must occur independently inside training folds.\n")

cat("\nOutputs:\n")
cat("- results/week2_differential_expression_exploratory.csv\n")
cat("- results/week2_differential_expression_summary.csv\n")
