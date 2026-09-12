# ============================================================================
# Week 2 - Preprocessing foundation
# GSE160638 primary anti-PD-1 cohort
#
# Purpose:
#   - consume the validated Week 1 analytical cohort;
#   - perform cohort-level exploratory preprocessing;
#   - generate auditable preprocessing outputs;
#   - preserve raw aligned counts as the modelling source of truth.
#
# IMPORTANT:
#   The exploratory filtering/TMM/logCPM objects produced here MUST NOT be used
#   directly as globally preprocessed input for predictive cross-validation.
#   Model-development preprocessing must be fitted inside training folds.
# ============================================================================


# ----------------------------------------------------------------------------
# 0. Dependencies
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
})


# ----------------------------------------------------------------------------
# 1. Paths
# ----------------------------------------------------------------------------

data_processed_dir <- "data_processed"
results_dir <- "results"

counts_file <- file.path(
  data_processed_dir,
  "GSE160638_raw_counts_PD1_aligned.rds"
)

clinical_file <- file.path(
  data_processed_dir,
  "GSE160638_TableS2A_PD1_clinical.csv"
)

required_files <- c(
  counts_file,
  clinical_file
)

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

clinical <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(counts)


# ----------------------------------------------------------------------------
# 3. Input integrity checks
# ----------------------------------------------------------------------------

if (!is.numeric(counts)) {
  stop("Count matrix is not numeric.", call. = FALSE)
}

if (nrow(counts) != 17002L) {
  stop(
    paste0(
      "Unexpected feature count: ",
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
      "Unexpected clinical cohort size: ",
      nrow(clinical),
      ". Expected 36."
    ),
    call. = FALSE
  )
}

if (!"patient_id" %in% names(clinical)) {
  stop("Clinical table does not contain patient_id.", call. = FALSE)
}

if (!"response" %in% names(clinical)) {
  stop("Clinical table does not contain response.", call. = FALSE)
}

if (!identical(colnames(counts), clinical$patient_id)) {
  stop(
    "Count-matrix columns and clinical patient identifiers are not in identical order.",
    call. = FALSE
  )
}

if (anyDuplicated(colnames(counts)) > 0L) {
  stop("Duplicated count-matrix sample identifiers detected.", call. = FALSE)
}

if (anyDuplicated(clinical$patient_id) > 0L) {
  stop("Duplicated clinical patient identifiers detected.", call. = FALSE)
}

if (anyNA(counts)) {
  stop("Missing values detected in the raw count matrix.", call. = FALSE)
}

if (any(counts < 0)) {
  stop("Negative values detected in the raw count matrix.", call. = FALSE)
}

non_integer_like <- sum(
  abs(counts - round(counts)) > 1e-8,
  na.rm = TRUE
)

if (non_integer_like != 0L) {
  stop(
    paste0(
      "Non-integer-like values detected in raw counts: ",
      non_integer_like
    ),
    call. = FALSE
  )
}

expected_response_counts <- c(
  NonResponder = 14L,
  Responder = 22L
)

observed_response_counts <- table(
  factor(
    clinical$response,
    levels = names(expected_response_counts)
  )
)

if (!identical(
  as.integer(observed_response_counts),
  as.integer(expected_response_counts)
)) {
  stop(
    paste0(
      "Response distribution does not match the validated Week 1 cohort. ",
      "Observed: ",
      paste(
        names(expected_response_counts),
        as.integer(observed_response_counts),
        sep = "=",
        collapse = ", "
      ),
      ". Expected: ",
      paste(
        names(expected_response_counts),
        as.integer(expected_response_counts),
        sep = "=",
        collapse = ", "
      ),
      "."
    ),
    call. = FALSE
  )
}


# ----------------------------------------------------------------------------
# 4. Create exploratory DGE object
# ----------------------------------------------------------------------------
#
# No response/group variable is supplied to filterByExpr().
#
# This avoids using the treatment-response endpoint to define the global
# exploratory gene set.
#
# Even so, this cohort-level filtering is for exploratory/descriptive analysis
# only. Predictive modelling must repeat preprocessing within training folds.
# ----------------------------------------------------------------------------

group <- factor(
  clinical$response,
  levels = c("NonResponder", "Responder")
)

dge <- DGEList(
  counts = counts,
  group = group
)

keep_exploratory <- filterByExpr(
  dge,
  group = group
)

dge_exploratory <- dge[
  keep_exploratory,
  ,
  keep.lib.sizes = FALSE
]


# ----------------------------------------------------------------------------
# 5. TMM normalization - exploratory only
# ----------------------------------------------------------------------------

dge_exploratory <- calcNormFactors(
  dge_exploratory,
  method = "TMM"
)


# ----------------------------------------------------------------------------
# 6. logCPM transformation - exploratory only
# ----------------------------------------------------------------------------

logcpm_exploratory <- cpm(
  dge_exploratory,
  log = TRUE,
  prior.count = 1
)


# ----------------------------------------------------------------------------
# 7. Auditable preprocessing outputs
# ----------------------------------------------------------------------------

gene_filtering_audit <- data.frame(
  gene_id = rownames(counts),
  keep_exploratory = keep_exploratory,
  stringsAsFactors = FALSE
)

normalization_audit <- data.frame(
  sample_name = colnames(dge_exploratory$counts),
  library_size = dge_exploratory$samples$lib.size,
  norm_factor = dge_exploratory$samples$norm.factors,
  effective_library_size =
    dge_exploratory$samples$lib.size *
    dge_exploratory$samples$norm.factors,
  stringsAsFactors = FALSE
)

preprocessing_summary <- data.frame(
  metric = c(
    "input_genes",
    "input_samples",
    "exploratory_genes_retained",
    "exploratory_genes_removed",
    "responders",
    "nonresponders"
  ),
  value = c(
    nrow(counts),
    ncol(counts),
    sum(keep_exploratory),
    sum(!keep_exploratory),
    sum(clinical$response == "Responder"),
    sum(clinical$response == "NonResponder")
  ),
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------------
# 8. Save outputs
# ----------------------------------------------------------------------------

saveRDS(
  logcpm_exploratory,
  file.path(
    data_processed_dir,
    "GSE160638_logCPM_exploratory.rds"
  )
)

write.csv(
  gene_filtering_audit,
  file.path(
    results_dir,
    "week2_gene_filtering_exploratory.csv"
  ),
  row.names = FALSE
)

write.csv(
  normalization_audit,
  file.path(
    results_dir,
    "week2_tmm_normalization_exploratory.csv"
  ),
  row.names = FALSE
)

write.csv(
  preprocessing_summary,
  file.path(
    results_dir,
    "week2_preprocessing_summary.csv"
  ),
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 9. Final validation
# ----------------------------------------------------------------------------

if (ncol(logcpm_exploratory) != 36L) {
  stop(
    "Exploratory expression matrix does not contain 36 samples.",
    call. = FALSE
  )
}

if (!identical(colnames(logcpm_exploratory), clinical$patient_id)) {
  stop(
    "Exploratory expression matrix lost validated Week 1 sample order.",
    call. = FALSE
  )
}


# ----------------------------------------------------------------------------
# 10. Report
# ----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("WEEK 2 PREPROCESSING FOUNDATION: PASS\n")
cat("============================================================\n")

cat("Validated input samples:", ncol(counts), "\n")
cat("Validated input genes:", nrow(counts), "\n")

cat(
  "Exploratory genes retained:",
  sum(keep_exploratory),
  "\n"
)

cat(
  "Exploratory genes removed:",
  sum(!keep_exploratory),
  "\n"
)

cat(
  "Responders:",
  sum(clinical$response == "Responder"),
  "\n"
)

cat(
  "NonResponders:",
  sum(clinical$response == "NonResponder"),
  "\n"
)

cat("Sample order preserved: TRUE\n")
cat("TMM normalization generated: TRUE\n")
cat("Exploratory logCPM generated: TRUE\n")

cat("\n")
cat("IMPORTANT:\n")
cat(
  "The exploratory filtered/TMM/logCPM matrix is not a modelling input.\n"
)
cat(
  "Predictive preprocessing must be fitted independently inside training folds.\n"
)

cat("\nOutputs:\n")
cat(
  "- data_processed/GSE160638_logCPM_exploratory.rds\n"
)
cat(
  "- results/week2_gene_filtering_exploratory.csv\n"
)
cat(
  "- results/week2_tmm_normalization_exploratory.csv\n"
)
cat(
  "- results/week2_preprocessing_summary.csv\n"
)