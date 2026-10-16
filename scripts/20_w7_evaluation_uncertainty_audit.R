options(stringsAsFactors = FALSE)

cat("W7 — Evaluation & uncertainty audit\n")
cat("Phase: frozen input manifest + outer-fold variability\n\n")

# -------------------------------------------------------------------------
# Frozen inputs
# -------------------------------------------------------------------------

inputs <- data.frame(
  path = c(
    "results/week3_nested_cv_summary_strict.csv",
    "results/week3_outer_fold_metrics_strict.csv",
    "results/week3_outer_predictions_strict.csv",
    "results/week4_external_validation_metrics_strict.csv",
    "results/week4_external_validation_predictions_strict.csv"
  ),
  role = c(
    "W3 pooled strict nested-CV summary",
    "W3 strict outer-fold metrics",
    "W3 frozen out-of-fold predictions",
    "W4 locked external-validation metrics",
    "W4 locked external-validation predictions"
  ),
  expected_rows = c(10L, 5L, 36L, 3L, 76L),
  stringsAsFactors = FALSE
)

stopifnot(all(file.exists(inputs$path)))


# -------------------------------------------------------------------------
# Portable SHA256 helper: macOS shasum / GNU sha256sum
# -------------------------------------------------------------------------

sha256_file <- function(path) {

  shasum <- Sys.which("shasum")
  sha256sum <- Sys.which("sha256sum")

  if (nzchar(shasum)) {

    out <- system2(
      shasum,
      c("-a", "256", path),
      stdout = TRUE,
      stderr = TRUE
    )

  } else if (nzchar(sha256sum)) {

    out <- system2(
      sha256sum,
      path,
      stdout = TRUE,
      stderr = TRUE
    )

  } else {

    stop(
      "Neither shasum nor sha256sum is available; ",
      "cannot generate SHA256 manifest."
    )
  }

  status <- attr(out, "status")

  if (!is.null(status) && status != 0L) {
    stop("SHA256 command failed for: ", path)
  }

  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}


# -------------------------------------------------------------------------
# Manifest
# -------------------------------------------------------------------------

observed_rows <- vapply(
  inputs$path,
  function(path) {
    nrow(
      read.csv(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )
  },
  integer(1)
)

stopifnot(
  length(observed_rows) == length(inputs$expected_rows),
  all(unname(observed_rows) == inputs$expected_rows)
)

manifest <- data.frame(
  path = inputs$path,
  role = inputs$role,
  rows = observed_rows,
  bytes = as.numeric(file.info(inputs$path)$size),
  sha256 = vapply(inputs$path, sha256_file, character(1)),
  frozen = TRUE,
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  "results/w7_frozen_input_manifest.csv",
  row.names = FALSE,
  quote = TRUE
)

cat("Frozen input manifest: PASS\n")
cat("Inputs:", nrow(manifest), "\n\n")


# -------------------------------------------------------------------------
# W3 outer-fold variability
# -------------------------------------------------------------------------

folds <- read.csv(
  "results/week3_outer_fold_metrics_strict.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(folds) == 5L,
  identical(sort(as.integer(folds$outer_fold)), 1:5)
)

w3_summary_raw <- read.csv(
  "results/week3_nested_cv_summary_strict.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

w3_summary <- setNames(
  as.numeric(w3_summary_raw$value),
  w3_summary_raw$metric
)

metrics <- c(
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "Balanced_Accuracy",
  "AUC"
)

summary_keys <- c(
  Accuracy = "accuracy",
  Sensitivity = "sensitivity",
  Specificity = "specificity",
  Balanced_Accuracy = "balanced_accuracy",
  AUC = "auc"
)

stopifnot(
  all(metrics %in% names(folds)),
  all(summary_keys %in% names(w3_summary))
)

variability <- do.call(
  rbind,
  lapply(
    metrics,
    function(metric) {

      x <- as.numeric(folds[[metric]])
      pooled <- unname(
        w3_summary[[summary_keys[[metric]]]]
      )

      data.frame(
        metric = metric,
        outer_folds = length(x),
        fold_mean = mean(x),
        fold_sd = sd(x),
        fold_median = median(x),
        fold_min = min(x),
        fold_max = max(x),
        fold_range = max(x) - min(x),
        pooled_oof_metric = pooled,
        fold_mean_minus_pooled = mean(x) - pooled,
        stringsAsFactors = FALSE
      )
    }
  )
)

numeric_cols <- setdiff(
  names(variability),
  c("metric", "outer_folds")
)

variability[numeric_cols] <- lapply(
  variability[numeric_cols],
  function(x) round(x, 12)
)

write.csv(
  variability,
  "results/w7_outer_fold_variability.csv",
  row.names = FALSE,
  quote = TRUE
)


# -------------------------------------------------------------------------
# Methodological contract
# -------------------------------------------------------------------------

doc <- c(
  "# W7 — Evaluation and uncertainty audit",
  "",
  "## Status",
  "",
  "W7 is an additive post-hoc audit of frozen Week 3 and Week 4 predictions and metrics.",
  "",
  "## Objective",
  "",
  paste(
    "Quantify predictive-performance variability and uncertainty without",
    "retraining models, reselection of genes, threshold optimization,",
    "cohort modification or retrospective alteration of frozen analytical outputs."
  ),
  "",
  "## Frozen inputs",
  "",
  "- `results/week3_nested_cv_summary_strict.csv`",
  "- `results/week3_outer_fold_metrics_strict.csv`",
  "- `results/week3_outer_predictions_strict.csv`",
  "- `results/week4_external_validation_metrics_strict.csv`",
  "- `results/week4_external_validation_predictions_strict.csv`",
  "",
  "Their SHA256 digests are recorded in `results/w7_frozen_input_manifest.csv`.",
  "",
  "## Non-negotiable constraints",
  "",
  "- No model retraining.",
  "- No feature or gene reselection.",
  "- No modification of the frozen Week 3 nested-CV design.",
  "- No modification of the locked Week 4 deployment model.",
  "- No threshold retuning.",
  "- No addition, removal or relabelling of samples.",
  "- No retrospective modification of W1-W4 analytical outputs.",
  "",
  "## W7.3 outer-fold variability",
  "",
  paste(
    "The five strict outer folds are summarized using mean, sample standard",
    "deviation, median, minimum, maximum and range."
  ),
  "",
  paste(
    "These quantities describe empirical between-fold variability only.",
    "They are not confidence intervals and are not treated as five independent",
    "estimates of a population parameter."
  ),
  "",
  paste(
    "The pooled out-of-fold metric is reported separately because a metric",
    "computed after pooling all OOF predictions is not generally identical",
    "to the unweighted mean of fold-specific metrics."
  ),
  "",
  "## Planned W7 analyses",
  "",
  "1. Frozen-input provenance and checksum manifest.",
  "2. Outer-fold performance variability.",
  "3. No-skill and prevalence-aware reference audit.",
  "4. Prediction-level metric uncertainty.",
  "5. Internal-versus-external transportability audit.",
  "6. Final W7 QA and clean-room reproduction.",
  "",
  "## Interpretation boundary",
  "",
  paste(
    "W7 characterizes uncertainty in already-frozen results.",
    "It does not attempt to improve predictive performance."
  ),
  ""
)

writeLines(
  doc,
  "docs/w7_evaluation_uncertainty_audit.md",
  useBytes = TRUE
)


# -------------------------------------------------------------------------
# Console summary
# -------------------------------------------------------------------------

cat("Outer-fold variability: PASS\n\n")

print(
  variability,
  row.names = FALSE
)

cat("\n")
cat("W7.3 OUTPUTS\n")
cat("- results/w7_frozen_input_manifest.csv\n")
cat("- results/w7_outer_fold_variability.csv\n")
cat("- docs/w7_evaluation_uncertainty_audit.md\n")

cat("\n")
cat("W7.3: PASS\n")
