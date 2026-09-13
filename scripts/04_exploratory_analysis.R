# =============================================================================
# Week 2 — Exploratory analysis
# =============================================================================
#
# Purpose
# -------
# Perform descriptive exploratory analysis of the validated GSE160638 PD1 cohort
# using the Week 2 exploratory TMM/logCPM expression matrix.
#
# IMPORTANT
# ---------
# This script is descriptive only.
#
# The filtered/TMM/logCPM matrix used here MUST NOT be used directly as the
# predictive-modelling input. Predictive preprocessing must be fitted
# independently inside the training folds.
#
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

data_processed_dir <- "data_processed"
results_dir <- "results"
figures_dir <- "figures"

expression_file <- file.path(
  data_processed_dir,
  "GSE160638_logCPM_exploratory.rds"
)

clinical_file <- file.path(
  data_processed_dir,
  "GSE160638_TableS2A_PD1_clinical.csv"
)

required_files <- c(
  expression_file,
  clinical_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    paste0(
      "Missing Week 2 exploratory input files:\n  ",
      paste(missing_files, collapse = "\n  ")
    ),
    call. = FALSE
  )
}

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. Load validated exploratory inputs
# -----------------------------------------------------------------------------

expr <- readRDS(expression_file)

clinical <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr <- as.matrix(expr)

# -----------------------------------------------------------------------------
# 3. Input integrity checks
# -----------------------------------------------------------------------------

if (!is.numeric(expr)) {
  stop("Exploratory expression matrix is not numeric.", call. = FALSE)
}

if (nrow(expr) != 15097L) {
  stop(
    paste0(
      "Unexpected exploratory feature count: ",
      nrow(expr),
      ". Expected 15097."
    ),
    call. = FALSE
  )
}

if (ncol(expr) != 36L) {
  stop(
    paste0(
      "Unexpected exploratory sample count: ",
      ncol(expr),
      ". Expected 36."
    ),
    call. = FALSE
  )
}

if (nrow(clinical) != 36L) {
  stop(
    paste0(
      "Unexpected clinical sample count: ",
      nrow(clinical),
      ". Expected 36."
    ),
    call. = FALSE
  )
}

if (!all(c("patient_id", "response") %in% colnames(clinical))) {
  stop(
    "Clinical file must contain patient_id and response columns.",
    call. = FALSE
  )
}

if (anyNA(expr)) {
  stop("Missing values detected in exploratory expression matrix.", call. = FALSE)
}

if (any(!is.finite(expr))) {
  stop(
    "Non-finite values detected in exploratory expression matrix.",
    call. = FALSE
  )
}

if (!identical(colnames(expr), clinical$patient_id)) {
  stop(
    "Exploratory expression matrix and clinical data have different sample order.",
    call. = FALSE
  )
}

expected_response_counts <- c(
  NonResponder = 14L,
  Responder = 22L
)

if (!all(clinical$response %in% names(expected_response_counts))) {
  stop(
    paste0(
      "Unexpected response labels detected: ",
      paste(
        setdiff(
          unique(clinical$response),
          names(expected_response_counts)
        ),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

observed_response_counts <- table(
  factor(
    clinical$response,
    levels = names(expected_response_counts)
  )
)

if (
  !identical(
    as.integer(observed_response_counts),
    unname(expected_response_counts)
  )
) {
  stop(
    paste0(
      "Response distribution does not match the validated Week 1 cohort. ",
      "Observed: NonResponder=",
      observed_response_counts["NonResponder"],
      ", Responder=",
      observed_response_counts["Responder"],
      "."
    ),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 4. Principal component analysis
# -----------------------------------------------------------------------------

# Samples are observations and genes are variables.
#
# The logCPM matrix is centred gene-wise. No additional unit-variance scaling is
# applied because the purpose is descriptive transcriptomic PCA and scaling
# low-variance genes to unit variance can amplify noise.

pca <- prcomp(
  t(expr),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_scores <- data.frame(
  patient_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  response = clinical$response,
  stringsAsFactors = FALSE
)

if (!identical(pca_scores$patient_id, clinical$patient_id)) {
  stop(
    "PCA scores lost the validated sample order.",
    call. = FALSE
  )
}

pca_variance <- data.frame(
  component = paste0("PC", seq_along(variance_explained)),
  variance_percent = variance_explained,
  cumulative_variance_percent = cumsum(variance_explained),
  stringsAsFactors = FALSE
)

write.csv(
  pca_scores,
  file.path(results_dir, "week2_pca_scores_exploratory.csv"),
  row.names = FALSE
)

write.csv(
  pca_variance,
  file.path(results_dir, "week2_pca_variance_exploratory.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. PCA figure
# -----------------------------------------------------------------------------

pca_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = response
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(
    title = "GSE160638 — exploratory PCA",
    subtitle = "Week 2 filtered and TMM-normalized logCPM matrix",
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2]),
    color = "Response"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "week2_pca_exploratory.png"
  ),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 6. Sample correlation
# -----------------------------------------------------------------------------

sample_correlation <- cor(
  expr,
  method = "pearson",
  use = "pairwise.complete.obs"
)

if (
  nrow(sample_correlation) != 36L ||
  ncol(sample_correlation) != 36L
) {
  stop(
    "Unexpected sample-correlation matrix dimensions.",
    call. = FALSE
  )
}

write.csv(
  sample_correlation,
  file.path(
    results_dir,
    "week2_sample_correlation_exploratory.csv"
  ),
  row.names = TRUE
)

# Use base R here so that Week 2 does not require an additional visualization
# dependency solely for this diagnostic figure.

png(
  filename = file.path(
    figures_dir,
    "week2_sample_correlation_exploratory.png"
  ),
  width = 2200,
  height = 2000,
  res = 250
)

heatmap(
  sample_correlation,
  symm = TRUE,
  margins = c(7, 7),
  main = "GSE160638 — exploratory sample correlation"
)

dev.off()

# -----------------------------------------------------------------------------
# 7. Exploratory-analysis summary
# -----------------------------------------------------------------------------

exploratory_summary <- data.frame(
  metric = c(
    "samples",
    "genes",
    "responders",
    "nonresponders",
    "missing_expression_values",
    "finite_expression_values",
    "pc1_variance_percent",
    "pc2_variance_percent"
  ),
  value = c(
    ncol(expr),
    nrow(expr),
    sum(clinical$response == "Responder"),
    sum(clinical$response == "NonResponder"),
    sum(is.na(expr)),
    all(is.finite(expr)),
    variance_explained[1],
    variance_explained[2]
  ),
  stringsAsFactors = FALSE
)

write.csv(
  exploratory_summary,
  file.path(
    results_dir,
    "week2_exploratory_analysis_summary.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8. Final report
# -----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("WEEK 2 EXPLORATORY ANALYSIS: PASS\n")
cat("============================================================\n")

cat("Samples:", ncol(expr), "\n")
cat("Genes:", nrow(expr), "\n")
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
cat("Missing expression values:", sum(is.na(expr)), "\n")
cat("Finite expression values:", all(is.finite(expr)), "\n")
cat(
  "PC1 variance:",
  sprintf("%.2f%%", variance_explained[1]),
  "\n"
)
cat(
  "PC2 variance:",
  sprintf("%.2f%%", variance_explained[2]),
  "\n"
)

cat("\nIMPORTANT:\n")
cat(
  "These PCA and correlation analyses are descriptive exploratory analyses only.\n"
)
cat(
  "Predictive preprocessing must be fitted independently inside training folds.\n"
)

cat("\nOutputs:\n")
cat("- results/week2_pca_scores_exploratory.csv\n")
cat("- results/week2_pca_variance_exploratory.csv\n")
cat("- results/week2_sample_correlation_exploratory.csv\n")
cat("- results/week2_exploratory_analysis_summary.csv\n")
cat("- figures/week2_pca_exploratory.png\n")
cat("- figures/week2_sample_correlation_exploratory.png\n")