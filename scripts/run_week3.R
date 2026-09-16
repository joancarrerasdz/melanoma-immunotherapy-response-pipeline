#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================================\n")
cat("WEEK 3 — STRICT NESTED-CV REPRODUCTION\n")
cat("============================================================\n")

required_scripts <- c(
  "scripts/06_nested_cv_strict.R",
  "scripts/07_gene_stability_strict.R",
  "scripts/08_compare_historical_strict_stability.R"
)

missing_scripts <- required_scripts[!file.exists(required_scripts)]

if (length(missing_scripts) > 0L) {
  stop(
    paste(
      "Missing Week 3 scripts:",
      paste(missing_scripts, collapse = "\n")
    ),
    call. = FALSE
  )
}

run_step <- function(step, script) {

  cat("\n")
  cat("============================================================\n")
  cat(step, "\n")
  cat("============================================================\n")

  status <- system2(
    command = "Rscript",
    args = c("--vanilla", script)
  )

  if (!identical(status, 0L)) {
    stop(
      paste("Week 3 step failed:", script),
      call. = FALSE
    )
  }

  cat("[PASS]", step, "\n")
}

run_step(
  "STEP 1/3 — Strict nested cross-validation",
  "scripts/06_nested_cv_strict.R"
)

run_step(
  "STEP 2/3 — Strict gene-selection stability",
  "scripts/07_gene_stability_strict.R"
)

run_step(
  "STEP 3/3 — Historical-versus-strict stability audit",
  "scripts/08_compare_historical_strict_stability.R"
)

# ------------------------------------------------------------
# Expected outputs
# ------------------------------------------------------------

expected_outputs <- c(
  "results/week3_inner_cv_metrics_strict.csv",
  "results/week3_inner_tuning_summary_strict.csv",
  "results/week3_outer_predictions_strict.csv",
  "results/week3_outer_fold_metrics_strict.csv",
  "results/week3_outer_selected_genes_strict.csv",
  "results/week3_gene_selection_frequency_strict.csv",
  "results/week3_nested_cv_summary_strict.csv",
  "figures/week3_nested_cv_roc_strict.png",

  "results/week3_gene_stability_tiers_strict.csv",
  "results/week3_consensus_repeated_2of5_strict.csv",
  "results/week3_consensus_majority_3of5_strict.csv",
  "results/week3_consensus_high_4of5_strict.csv",
  "results/week3_consensus_core_5of5_strict.csv",
  "figures/week3_gene_selection_stability_strict.png",

  "results/week3_historical_vs_strict_5of5.csv",
  "results/week3_historical_vs_strict_5of5_summary.csv",
  "results/week3_historical_13_genes_in_strict_pipeline.csv",
  "results/week3_historical_13_genes_in_strict_pipeline_summary.csv"
)

missing_outputs <- expected_outputs[
  !file.exists(expected_outputs)
]

if (length(missing_outputs) > 0L) {
  stop(
    paste(
      "Expected Week 3 outputs are missing:",
      paste(missing_outputs, collapse = "\n")
    ),
    call. = FALSE
  )
}

empty_outputs <- expected_outputs[
  file.info(expected_outputs)$size <= 0
]

if (length(empty_outputs) > 0L) {
  stop(
    paste(
      "Week 3 produced empty output files:",
      paste(empty_outputs, collapse = "\n")
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Strict nested-CV validation
# ------------------------------------------------------------

nested_summary <- read.csv(
  "results/week3_nested_cv_summary_strict.csv",
  stringsAsFactors = FALSE
)

get_metric <- function(df, metric) {

  value <- df$value[df$metric == metric]

  if (length(value) != 1L) {
    stop(
      paste("Expected exactly one value for metric:", metric),
      call. = FALSE
    )
  }

  as.numeric(value)
}

stopifnot(
  get_metric(nested_summary, "samples") == 36,
  get_metric(nested_summary, "input_genes") == 17002,
  get_metric(nested_summary, "outer_folds") == 5,
  get_metric(nested_summary, "inner_folds") == 3,
  get_metric(nested_summary, "genes_per_model") == 100
)

expected_performance <- c(
  accuracy = 0.8333,
  sensitivity = 0.9091,
  specificity = 0.7143,
  balanced_accuracy = 0.8117,
  auc = 0.8896
)

for (metric in names(expected_performance)) {

  observed <- get_metric(
    nested_summary,
    metric
  )

  expected <- expected_performance[[metric]]

  if (abs(observed - expected) > 1e-4) {
    stop(
      paste(
        "Unexpected",
        metric,
        ": observed",
        observed,
        "expected",
        expected
      ),
      call. = FALSE
    )
  }
}

outer_predictions <- read.csv(
  "results/week3_outer_predictions_strict.csv",
  stringsAsFactors = FALSE
)

if (nrow(outer_predictions) != 36L) {
  stop(
    "Outer predictions do not contain exactly 36 samples.",
    call. = FALSE
  )
}

if (anyDuplicated(outer_predictions$sample_name)) {
  stop(
    "Duplicate samples detected in outer predictions.",
    call. = FALSE
  )
}

if (anyNA(outer_predictions$prob_responder)) {
  stop(
    "Missing response probabilities detected.",
    call. = FALSE
  )
}

if (
  any(outer_predictions$prob_responder < 0) ||
  any(outer_predictions$prob_responder > 1)
) {
  stop(
    "Invalid response probabilities detected.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Gene-stability validation
# ------------------------------------------------------------

stability <- read.csv(
  "results/week3_gene_stability_tiers_strict.csv",
  stringsAsFactors = FALSE
)

stability_counts <- c(
  selected_1of5 = nrow(stability),
  selected_2of5 = sum(stability$outer_folds_selected >= 2),
  selected_3of5 = sum(stability$outer_folds_selected >= 3),
  selected_4of5 = sum(stability$outer_folds_selected >= 4),
  selected_5of5 = sum(stability$outer_folds_selected == 5)
)

expected_stability <- c(
  selected_1of5 = 312,
  selected_2of5 = 100,
  selected_3of5 = 51,
  selected_4of5 = 25,
  selected_5of5 = 12
)

if (!identical(
  unname(as.integer(stability_counts)),
  unname(as.integer(expected_stability))
)) {
  stop(
    "Strict gene-stability counts differ from the validated Week 3 results.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Historical-versus-strict audit validation
# ------------------------------------------------------------

historical_audit <- read.csv(
  "results/week3_historical_13_genes_in_strict_pipeline_summary.csv",
  stringsAsFactors = FALSE
)

historical_value <- function(metric) {
  get_metric(historical_audit, metric)
}

stopifnot(
  historical_value("historical_genes") == 13,
  historical_value("present_in_validated_input") == 13,
  historical_value("selected_at_least_1of5") == 0,
  historical_value("selected_at_least_2of5") == 0,
  historical_value("selected_at_least_3of5") == 0,
  historical_value("selected_at_least_4of5") == 0,
  historical_value("selected_5of5") == 0
)

core_comparison <- read.csv(
  "results/week3_historical_vs_strict_5of5_summary.csv",
  stringsAsFactors = FALSE
)

core_value <- function(metric) {
  get_metric(core_comparison, metric)
}

stopifnot(
  core_value("historical_5of5_genes") == 13,
  core_value("strict_5of5_genes") == 12,
  core_value("shared_genes") == 0,
  core_value("historical_only_genes") == 13,
  core_value("strict_only_genes") == 12
)

# ------------------------------------------------------------
# Final gate
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("WEEK 3 REPRODUCTION: PASS\n")
cat("============================================================\n")

cat("Validated samples: 36\n")
cat("Validated input genes: 17002\n")
cat("Outer folds: 5\n")
cat("Inner folds: 3\n")
cat("Genes selected per training partition: 100\n")

cat("\nValidation performance:\n")
cat("- Accuracy: 0.8333\n")
cat("- Sensitivity: 0.9091\n")
cat("- Specificity: 0.7143\n")
cat("- Balanced Accuracy: 0.8117\n")
cat("- AUC: 0.8896\n")

cat("\nStrict gene-selection stability:\n")
cat("- Selected >=1/5 outer folds: 312\n")
cat("- Selected >=2/5 outer folds: 100\n")
cat("- Selected >=3/5 outer folds: 51\n")
cat("- Selected >=4/5 outer folds: 25\n")
cat("- Selected 5/5 outer folds: 12\n")

cat("\nHistorical reproducibility audit:\n")
cat("- Historical 5/5 genes: 13\n")
cat("- Historical genes present in validated input: 13/13\n")
cat("- Historical genes selected >=1/5 strict folds: 0/13\n")
cat("- Historical / strict 5/5 overlap: 0\n")

cat("\nExpected Week 3 outputs: present and non-empty\n")

cat("\nIMPORTANT:\n")
cat(
  "Validation performance is calculated exclusively from pooled outer ",
  "out-of-fold predictions.\n",
  sep = ""
)

cat(
  "Feature selection and hyperparameter tuning are restricted to training ",
  "partitions.\n",
  sep = ""
)

cat(
  "Consensus gene sets are stability-derived candidate signatures and are ",
  "not independently validated final signatures.\n",
  sep = ""
)

cat(
  "Historical full-dataset Random Forest importance is exploratory and is ",
  "not used to define the strict consensus candidates.\n",
  sep = ""
)
