# ==============================================================================
# Week 2 reproduction runner
# ==============================================================================
#
# Reproduces the Week 2 exploratory preprocessing and analysis workflow from
# the validated Week 1 data foundation.
#
# IMPORTANT:
# - Exploratory filtering, TMM normalization, logCPM transformation, PCA,
#   sample correlation and differential-expression analysis are descriptive.
# - The exploratory expression matrix and differential-expression results are
#   NOT predictive modelling inputs.
# - Predictive preprocessing and feature selection must be fitted independently
#   inside training folds.
#
# Run from the repository root with:
#
#   Rscript --vanilla scripts/run_week2.R
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Validate execution context
# ------------------------------------------------------------------------------

required_scripts <- c(
  "scripts/03_preprocessing_train.R",
  "scripts/04_exploratory_analysis.R",
  "scripts/05_differential_expression_exploratory.R"
)

missing_scripts <- required_scripts[!file.exists(required_scripts)]

if (length(missing_scripts) > 0L) {
  stop(
    paste0(
      "Week 2 runner must be executed from the repository root.\n",
      "Missing scripts:\n  ",
      paste(missing_scripts, collapse = "\n  ")
    ),
    call. = FALSE
  )
}


required_inputs <- c(
  "data_processed/GSE160638_raw_counts_PD1_aligned.rds",
  "data_processed/GSE160638_TableS2A_PD1_clinical.csv"
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0L) {
  stop(
    paste0(
      "Validated Week 1 inputs are missing:\n  ",
      paste(missing_inputs, collapse = "\n  ")
    ),
    call. = FALSE
  )
}


required_packages <- c(
  "edgeR",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing Week 2 R dependencies: ",
      paste(missing_packages, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# 2. Helper for isolated script execution
# ------------------------------------------------------------------------------

run_step <- function(script, label) {

  cat("\n")
  cat("============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")

  rscript <- Sys.which("Rscript")

  if (!nzchar(rscript)) {
    stop(
      "Rscript executable could not be found.",
      call. = FALSE
    )
  }

  status <- system2(
    command = rscript,
    args = c("--vanilla", script)
  )

  if (!identical(status, 0L)) {
    stop(
      paste0(
        "Week 2 reproduction failed during: ",
        label,
        "\nScript: ",
        script,
        "\nExit status: ",
        status
      ),
      call. = FALSE
    )
  }

  cat("\n[PASS] ", label, "\n", sep = "")
}


# ------------------------------------------------------------------------------
# 3. Execute Week 2
# ------------------------------------------------------------------------------

run_step(
  "scripts/03_preprocessing_train.R",
  "STEP 1/3 - Exploratory preprocessing foundation"
)

run_step(
  "scripts/04_exploratory_analysis.R",
  "STEP 2/3 - Exploratory PCA and sample correlation"
)

run_step(
  "scripts/05_differential_expression_exploratory.R",
  "STEP 3/3 - Exploratory differential expression"
)


# ------------------------------------------------------------------------------
# 4. Validate expected outputs
# ------------------------------------------------------------------------------

expected_outputs <- c(
  "data_processed/GSE160638_logCPM_exploratory.rds",

  "results/week2_gene_filtering_exploratory.csv",
  "results/week2_tmm_normalization_exploratory.csv",
  "results/week2_preprocessing_summary.csv",

  "results/week2_pca_scores_exploratory.csv",
  "results/week2_pca_variance_exploratory.csv",
  "results/week2_sample_correlation_exploratory.csv",
  "results/week2_exploratory_analysis_summary.csv",

  "figures/week2_pca_exploratory.png",
  "figures/week2_sample_correlation_exploratory.png",

  "results/week2_differential_expression_exploratory.csv",
  "results/week2_differential_expression_summary.csv"
)

missing_outputs <- expected_outputs[
  !file.exists(expected_outputs)
]

if (length(missing_outputs) > 0L) {
  stop(
    paste0(
      "Week 2 execution completed but expected outputs are missing:\n  ",
      paste(missing_outputs, collapse = "\n  ")
    ),
    call. = FALSE
  )
}


empty_outputs <- expected_outputs[
  file.info(expected_outputs)$size <= 0
]

if (length(empty_outputs) > 0L) {
  stop(
    paste0(
      "Week 2 produced empty output files:\n  ",
      paste(empty_outputs, collapse = "\n  ")
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# 5. Validate central Week 2 results
# ------------------------------------------------------------------------------

preprocessing_summary <- read.csv(
  "results/week2_preprocessing_summary.csv",
  stringsAsFactors = FALSE
)

summary_value <- function(metric) {
  preprocessing_summary$value[
    preprocessing_summary$metric == metric
  ]
}

stopifnot(
  identical(as.numeric(summary_value("input_genes")), 17002),
  identical(as.numeric(summary_value("input_samples")), 36),
  identical(
    as.numeric(summary_value("exploratory_genes_retained")),
    15097
  ),
  identical(
    as.numeric(summary_value("exploratory_genes_removed")),
    1905
  ),
  identical(as.numeric(summary_value("responders")), 22),
  identical(as.numeric(summary_value("nonresponders")), 14)
)


expr <- readRDS(
  "data_processed/GSE160638_logCPM_exploratory.rds"
)

if (!identical(dim(expr), c(15097L, 36L))) {
  stop(
    paste0(
      "Unexpected exploratory expression dimensions: ",
      nrow(expr),
      " x ",
      ncol(expr),
      ". Expected 15097 x 36."
    ),
    call. = FALSE
  )
}

if (anyNA(expr)) {
  stop(
    "Missing values detected in exploratory expression matrix.",
    call. = FALSE
  )
}

if (!all(is.finite(expr))) {
  stop(
    "Non-finite values detected in exploratory expression matrix.",
    call. = FALSE
  )
}


pca_scores <- read.csv(
  "results/week2_pca_scores_exploratory.csv",
  stringsAsFactors = FALSE
)

sample_correlation <- read.csv(
  "results/week2_sample_correlation_exploratory.csv",
  row.names = 1,
  check.names = FALSE
)

if (nrow(pca_scores) != 36L) {
  stop(
    "PCA output does not contain exactly 36 samples.",
    call. = FALSE
  )
}

if (
  nrow(sample_correlation) != 36L ||
  ncol(sample_correlation) != 36L
) {
  stop(
    "Sample-correlation matrix is not 36 x 36.",
    call. = FALSE
  )
}


de_summary <- read.csv(
  "results/week2_differential_expression_summary.csv",
  stringsAsFactors = FALSE
)

de_value <- function(metric) {
  de_summary$value[
    de_summary$metric == metric
  ]
}


# ------------------------------------------------------------------------------
# 6. Final report
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("WEEK 2 REPRODUCTION: PASS\n")
cat("============================================================\n")

cat("Validated Week 1 input samples: 36\n")
cat("Validated Week 1 input genes: 17002\n")
cat("Responders: 22\n")
cat("NonResponders: 14\n")

cat("\nExploratory preprocessing:\n")
cat("- Genes retained after filterByExpr: 15097\n")
cat("- Genes removed: 1905\n")
cat("- TMM normalization: generated\n")
cat("- Exploratory logCPM matrix: 15097 x 36\n")

cat("\nExploratory analysis:\n")
cat("- PCA samples: 36\n")
cat("- Sample-correlation matrix: 36 x 36\n")

if ("fdr_lt_0.05" %in% de_summary$metric) {
  cat(
    "- Differential-expression genes with FDR < 0.05: ",
    de_value("fdr_lt_0.05"),
    "\n",
    sep = ""
  )
}

if ("fdr_lt_0.10" %in% de_summary$metric) {
  cat(
    "- Differential-expression genes with FDR < 0.10: ",
    de_value("fdr_lt_0.10"),
    "\n",
    sep = ""
  )
}

cat("\nExpected outputs: present and non-empty\n")

cat("\nIMPORTANT:\n")
cat(
  "The Week 2 filtered/TMM/logCPM matrix and differential-expression ",
  "results are exploratory only.\n",
  sep = ""
)
cat(
  "They must not be used as globally preprocessed or preselected inputs ",
  "for predictive modelling.\n",
  sep = ""
)
cat(
  "Predictive preprocessing and feature selection must be fitted ",
  "independently inside training folds.\n"
)