#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(pROC)
})

# ============================================================
# Week 4 — locked external validation
#
# IMPORTANT
# ---------
# This script evaluates predictions that were already frozen
# before external response labels were consulted.
#
# It does NOT:
#   - refit the model;
#   - reselect genes;
#   - change the representation;
#   - retune mtry;
#   - optimize or change the probability threshold.
#
# Primary interpretation is cohort-specific.
# Pooled external metrics are reported as secondary summaries.
# ============================================================

results_dir <- "results"
processed_dir <- "data_processed"
figures_dir <- "figures"

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  figures_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 1. Helpers
# ============================================================

safe_divide <- function(
  numerator,
  denominator
) {

  if (denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

spec_value <- function(
  spec,
  parameter
) {

  value <- spec$value[
    spec$parameter == parameter
  ]

  if (length(value) != 1L) {
    stop(
      "Expected exactly one model-spec value for: ",
      parameter,
      call. = FALSE
    )
  }

  value
}

validate_external_metadata <- function(
  metadata,
  dataset,
  expected_n,
  expected_nonresponders,
  expected_responders
) {

  required <- c(
    "dataset",
    "sample_id",
    "sample_name",
    "response_raw",
    "response",
    "treatment_time",
    "keep"
  )

  missing_columns <- setdiff(
    required,
    colnames(metadata)
  )

  if (length(missing_columns) > 0L) {
    stop(
      dataset,
      ": missing metadata columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (nrow(metadata) != expected_n) {
    stop(
      dataset,
      ": expected ",
      expected_n,
      " samples; found ",
      nrow(metadata),
      ".",
      call. = FALSE
    )
  }

  if (
    anyNA(metadata$sample_id) ||
    anyDuplicated(metadata$sample_id)
  ) {
    stop(
      dataset,
      ": invalid or duplicated sample IDs.",
      call. = FALSE
    )
  }

  if (!all(metadata$dataset == dataset)) {
    stop(
      dataset,
      ": dataset label mismatch.",
      call. = FALSE
    )
  }

  if (!all(
    as.logical(metadata$keep)
  )) {
    stop(
      dataset,
      ": non-kept samples found in kept metadata.",
      call. = FALSE
    )
  }

  if (!all(
    metadata$treatment_time ==
    "PreTreatment"
  )) {
    stop(
      dataset,
      ": non-pre-treatment samples detected.",
      call. = FALSE
    )
  }

  response <- factor(
    metadata$response,
    levels = c(
      "NonResponder",
      "Responder"
    )
  )

  if (anyNA(response)) {
    stop(
      dataset,
      ": missing or unexpected response labels.",
      call. = FALSE
    )
  }

  response_counts <- table(
    response
  )

  if (
    unname(
      response_counts["NonResponder"]
    ) != expected_nonresponders ||
    unname(
      response_counts["Responder"]
    ) != expected_responders
  ) {
    stop(
      dataset,
      ": response distribution differs from validated metadata.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

join_locked_predictions <- function(
  predictions,
  metadata,
  dataset
) {

  p <- predictions[
    predictions$dataset == dataset,
    ,
    drop = FALSE
  ]

  m <- metadata[
    metadata$dataset == dataset,
    ,
    drop = FALSE
  ]

  if (
    anyDuplicated(p$sample_id) ||
    anyDuplicated(m$sample_id)
  ) {
    stop(
      dataset,
      ": duplicated sample IDs detected.",
      call. = FALSE
    )
  }

  if (!setequal(
    p$sample_id,
    m$sample_id
  )) {

    predictions_without_metadata <- setdiff(
      p$sample_id,
      m$sample_id
    )

    metadata_without_predictions <- setdiff(
      m$sample_id,
      p$sample_id
    )

    if (
      length(
        predictions_without_metadata
      ) > 0L
    ) {

      cat(
        "\nPredictions without metadata:\n"
      )

      print(
        predictions_without_metadata
      )
    }

    if (
      length(
        metadata_without_predictions
      ) > 0L
    ) {

      cat(
        "\nMetadata without predictions:\n"
      )

      print(
        metadata_without_predictions
      )
    }

    stop(
      dataset,
      ": prediction/metadata sample sets differ.",
      call. = FALSE
    )
  }

  idx <- match(
    p$sample_id,
    m$sample_id
  )

  if (anyNA(idx)) {
    stop(
      dataset,
      ": failed to align metadata.",
      call. = FALSE
    )
  }

  data.frame(
    dataset = p$dataset,
    sample_id = p$sample_id,
    sample_name = m$sample_name[idx],
    response_raw = m$response_raw[idx],
    observed = m$response[idx],
    prob_responder = p$prob_responder,
    predicted = p$predicted,
    frozen_threshold = p$frozen_threshold,
    representation = p$representation,
    n_genes = p$n_genes,
    stringsAsFactors = FALSE
  )
}

evaluate_external <- function(
  x,
  scope,
  analysis_role
) {

  observed <- factor(
    x$observed,
    levels = c(
      "NonResponder",
      "Responder"
    )
  )

  predicted <- factor(
    x$predicted,
    levels = c(
      "NonResponder",
      "Responder"
    )
  )

  if (
    length(
      unique(observed)
    ) != 2L
  ) {
    stop(
      scope,
      ": both response classes are required.",
      call. = FALSE
    )
  }

  TP <- sum(
    observed == "Responder" &
    predicted == "Responder"
  )

  FN <- sum(
    observed == "Responder" &
    predicted == "NonResponder"
  )

  TN <- sum(
    observed == "NonResponder" &
    predicted == "NonResponder"
  )

  FP <- sum(
    observed == "NonResponder" &
    predicted == "Responder"
  )

  accuracy <- safe_divide(
    TP + TN,
    nrow(x)
  )

  sensitivity <- safe_divide(
    TP,
    TP + FN
  )

  specificity <- safe_divide(
    TN,
    TN + FP
  )

  balanced_accuracy <- mean(
    c(
      sensitivity,
      specificity
    ),
    na.rm = TRUE
  )

  ppv <- safe_divide(
    TP,
    TP + FP
  )

  npv <- safe_divide(
    TN,
    TN + FN
  )

  f1 <- safe_divide(
    2 * TP,
    2 * TP + FP + FN
  )

  observed_numeric <- as.integer(
    observed == "Responder"
  )

  brier <- mean(
    (
      x$prob_responder -
      observed_numeric
    )^2
  )

  roc_object <- roc(
    response = observed,
    predictor = x$prob_responder,
    levels = c(
      "NonResponder",
      "Responder"
    ),
    direction = "<",
    quiet = TRUE
  )

  auc_value <- as.numeric(
    auc(
      roc_object
    )
  )

  auc_ci <- tryCatch(
    {
      as.numeric(
        ci.auc(
          roc_object,
          conf.level = 0.95,
          method = "delong"
        )
      )
    },
    error = function(e) {
      c(
        NA_real_,
        auc_value,
        NA_real_
      )
    }
  )

  metrics <- data.frame(
    scope = scope,
    analysis_role = analysis_role,
    samples = nrow(x),
    responders = sum(
      observed == "Responder"
    ),
    nonresponders = sum(
      observed == "NonResponder"
    ),
    threshold = unique(
      x$frozen_threshold
    ),
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Balanced_Accuracy = balanced_accuracy,
    PPV = ppv,
    NPV = npv,
    F1 = f1,
    AUC = auc_value,
    AUC_CI95_lower = auc_ci[1],
    AUC_CI95_upper = auc_ci[3],
    Brier_Score = brier,
    stringsAsFactors = FALSE
  )

  list(
    metrics = metrics,
    roc = roc_object
  )
}

# ============================================================
# 2. Verify the frozen Gate 4 artifacts BEFORE evaluation
# ============================================================

freeze_manifest_file <- file.path(
  results_dir,
  "week4_deployment_freeze_manifest_strict.csv"
)

freeze_manifest <- read.csv(
  freeze_manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_manifest_columns <- c(
  "file",
  "md5"
)

if (!all(
  required_manifest_columns %in%
  colnames(freeze_manifest)
)) {
  stop(
    "Invalid Gate 4 freeze manifest.",
    call. = FALSE
  )
}

if (!all(
  file.exists(
    freeze_manifest$file
  )
)) {

  missing_files <- freeze_manifest$file[
    !file.exists(
      freeze_manifest$file
    )
  ]

  stop(
    "Frozen files missing: ",
    paste(
      missing_files,
      collapse = ", "
    ),
    call. = FALSE
  )
}

current_md5 <- unname(
  tools::md5sum(
    freeze_manifest$file
  )
)

md5_match <- identical(
  as.character(current_md5),
  as.character(freeze_manifest$md5)
)

if (!md5_match) {

  audit <- data.frame(
    file = freeze_manifest$file,
    frozen_md5 = freeze_manifest$md5,
    current_md5 = current_md5,
    match = (
      current_md5 ==
      freeze_manifest$md5
    ),
    stringsAsFactors = FALSE
  )

  print(
    audit
  )

  stop(
    "Gate 4 freeze integrity check FAILED.",
    call. = FALSE
  )
}

freeze_commit <- tryCatch(
  {
    trimws(
      system2(
        "git",
        c(
          "log",
          "-1",
          "--format=%H",
          "--",
          "results/week4_external_predictions_blind_strict.csv"
        ),
        stdout = TRUE
      )
    )
  },
  error = function(e) {
    NA_character_
  }
)

cat("\n============================================================\n")
cat("GATE 4 FREEZE INTEGRITY\n")
cat("============================================================\n")

cat(
  "Frozen files verified:",
  nrow(freeze_manifest),
  "\n"
)

cat(
  "MD5 integrity:",
  md5_match,
  "\n"
)

cat(
  "Frozen blind-prediction commit:",
  freeze_commit,
  "\n"
)

# ============================================================
# 3. Read the already-frozen blind predictions
# ============================================================

predictions_file <- file.path(
  results_dir,
  "week4_external_predictions_blind_strict.csv"
)

predictions <- read.csv(
  predictions_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_prediction_columns <- c(
  "dataset",
  "sample_id",
  "prob_responder",
  "predicted",
  "frozen_threshold",
  "representation",
  "n_genes"
)

missing_prediction_columns <- setdiff(
  required_prediction_columns,
  colnames(predictions)
)

if (
  length(
    missing_prediction_columns
  ) > 0L
) {
  stop(
    "Missing blind-prediction columns: ",
    paste(
      missing_prediction_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  nrow(predictions) != 76L ||
  sum(
    predictions$dataset ==
    "GSE91061"
  ) != 49L ||
  sum(
    predictions$dataset ==
    "GSE78220"
  ) != 27L
) {
  stop(
    "Frozen external prediction counts differ from Gate 4.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    paste(
      predictions$dataset,
      predictions$sample_id,
      sep = "::"
    )
  )
) {
  stop(
    "Duplicated external prediction keys.",
    call. = FALSE
  )
}

if (
  anyNA(
    predictions$prob_responder
  ) ||
  any(
    predictions$prob_responder < 0 |
    predictions$prob_responder > 1
  )
) {
  stop(
    "Invalid frozen external probabilities.",
    call. = FALSE
  )
}

if (!all(
  predictions$representation ==
  "samplewise_rank"
)) {
  stop(
    "Frozen representation differs from samplewise_rank.",
    call. = FALSE
  )
}

if (!all(
  predictions$n_genes ==
  12L
)) {
  stop(
    "Frozen external prediction gene count differs from 12.",
    call. = FALSE
  )
}

# ============================================================
# 4. Verify locked model specification
# ============================================================

model_spec <- read.csv(
  file.path(
    results_dir,
    "week4_deployment_model_spec_strict.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

frozen_threshold <- as.numeric(
  spec_value(
    model_spec,
    "frozen_threshold"
  )
)

selected_mtry <- as.integer(
  spec_value(
    model_spec,
    "selected_mtry"
  )
)

ntree <- as.integer(
  spec_value(
    model_spec,
    "ntree"
  )
)

if (
  frozen_threshold != 0.48 ||
  selected_mtry != 1L ||
  ntree != 500L
) {
  stop(
    "Frozen deployment specification differs from Gate 4.",
    call. = FALSE
  )
}

if (!all(
  abs(
    predictions$frozen_threshold -
    frozen_threshold
  ) < 1e-12
)) {
  stop(
    "Prediction threshold differs from frozen model specification.",
    call. = FALSE
  )
}

predicted_from_probability <- ifelse(
  predictions$prob_responder >=
    frozen_threshold,
  "Responder",
  "NonResponder"
)

if (!identical(
  predicted_from_probability,
  predictions$predicted
)) {
  stop(
    "Frozen class labels do not match frozen probabilities and threshold.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("LOCKED MODEL SPECIFICATION\n")
cat("============================================================\n")

cat(
  "Genes: 12\n"
)

cat(
  "Representation: samplewise_rank\n"
)

cat(
  "mtry:",
  selected_mtry,
  "\n"
)

cat(
  "ntree:",
  ntree,
  "\n"
)

cat(
  "Frozen threshold:",
  frozen_threshold,
  "\n"
)

cat(
  "Model refitted during evaluation: FALSE\n"
)

cat(
  "Threshold optimized on external labels: FALSE\n"
)

# ============================================================
# 5. External labels are opened HERE for the first time
#    in the locked external-evaluation stage
# ============================================================

metadata_91061 <- read.csv(
  file.path(
    processed_dir,
    "metadata_GSE91061_harmonized_kept.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata_78220 <- read.csv(
  file.path(
    processed_dir,
    "metadata_GSE78220_harmonized_kept.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

validate_external_metadata(
  metadata_91061,
  dataset = "GSE91061",
  expected_n = 49L,
  expected_nonresponders = 39L,
  expected_responders = 10L
)

validate_external_metadata(
  metadata_78220,
  dataset = "GSE78220",
  expected_n = 27L,
  expected_nonresponders = 12L,
  expected_responders = 15L
)

cat("\n============================================================\n")
cat("EXTERNAL LABEL CONTRACT\n")
cat("============================================================\n")

cat(
  "GSE91061: 49 samples | 39 NonResponders | 10 Responders\n"
)

cat(
  "GSE78220: 27 samples | 12 NonResponders | 15 Responders\n"
)

# ============================================================
# 6. Locked join
# ============================================================

joined_91061 <- join_locked_predictions(
  predictions,
  metadata_91061,
  "GSE91061"
)

joined_78220 <- join_locked_predictions(
  predictions,
  metadata_78220,
  "GSE78220"
)

joined_external <- rbind(
  joined_91061,
  joined_78220
)

if (
  nrow(joined_external) != 76L ||
  anyNA(joined_external$observed)
) {
  stop(
    "Locked external join failed.",
    call. = FALSE
  )
}

write.csv(
  joined_external,
  file.path(
    results_dir,
    "week4_external_validation_predictions_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. External performance
#
# Primary:
#   dataset-specific results
#
# Secondary:
#   pooled results
# ============================================================

eval_91061 <- evaluate_external(
  joined_91061,
  scope = "GSE91061",
  analysis_role = "primary_external_cohort"
)

eval_78220 <- evaluate_external(
  joined_78220,
  scope = "GSE78220",
  analysis_role = "primary_external_cohort"
)

eval_pooled <- evaluate_external(
  joined_external,
  scope = "POOLED",
  analysis_role = "secondary_pooled_summary"
)

metrics <- rbind(
  eval_91061$metrics,
  eval_78220$metrics,
  eval_pooled$metrics
)

write.csv(
  metrics,
  file.path(
    results_dir,
    "week4_external_validation_metrics_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. ROC figure
# ============================================================

png(
  filename = file.path(
    figures_dir,
    "week4_external_validation_roc_strict.png"
  ),
  width = 2200,
  height = 1800,
  res = 250
)

plot(
  eval_91061$roc,
  legacy.axes = TRUE,
  lwd = 2,
  main = "Week 4 locked external validation — ROC"
)

plot(
  eval_78220$roc,
  add = TRUE,
  lwd = 2,
  lty = 2
)

plot(
  eval_pooled$roc,
  add = TRUE,
  lwd = 2,
  lty = 3
)

abline(
  a = 0,
  b = 1,
  lty = 4
)

legend(
  "bottomright",
  legend = c(
    paste0(
      "GSE91061 AUC = ",
      sprintf(
        "%.3f",
        eval_91061$metrics$AUC
      )
    ),
    paste0(
      "GSE78220 AUC = ",
      sprintf(
        "%.3f",
        eval_78220$metrics$AUC
      )
    ),
    paste0(
      "Pooled AUC = ",
      sprintf(
        "%.3f",
        eval_pooled$metrics$AUC
      )
    )
  ),
  lty = c(
    1,
    2,
    3
  ),
  lwd = 2,
  bty = "n"
)

dev.off()

# ============================================================
# 9. Evaluation audit
# ============================================================

evaluation_audit <- data.frame(
  item = c(
    "gate4_freeze_md5_verified",
    "blind_prediction_freeze_commit",
    "blind_predictions",
    "external_samples_joined",
    "genes",
    "representation",
    "selected_mtry",
    "ntree",
    "frozen_threshold",
    "model_refitted_after_labels_opened",
    "genes_reselected_after_labels_opened",
    "representation_changed_after_labels_opened",
    "mtry_changed_after_labels_opened",
    "threshold_changed_after_labels_opened",
    "threshold_optimized_on_external_labels",
    "primary_reporting",
    "pooled_reporting"
  ),
  value = c(
    "TRUE",
    freeze_commit,
    "76",
    "76",
    "12",
    "samplewise_rank",
    as.character(
      selected_mtry
    ),
    as.character(
      ntree
    ),
    as.character(
      frozen_threshold
    ),
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "cohort_specific",
    "secondary"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  evaluation_audit,
  file.path(
    results_dir,
    "week4_external_validation_audit_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10. Final gate output
# ============================================================

cat("\n============================================================\n")
cat("WEEK 4 LOCKED EXTERNAL VALIDATION: COMPLETE\n")
cat("============================================================\n")

for (
  i in seq_len(
    nrow(metrics)
  )
) {

  cat("\n")
  cat(
    metrics$scope[i],
    "\n"
  )

  cat(
    "Role:",
    metrics$analysis_role[i],
    "\n"
  )

  cat(
    "Samples:",
    metrics$samples[i],
    "\n"
  )

  cat(
    "Responders:",
    metrics$responders[i],
    "\n"
  )

  cat(
    "NonResponders:",
    metrics$nonresponders[i],
    "\n"
  )

  cat(
    "TP / TN / FP / FN:",
    metrics$TP[i],
    "/",
    metrics$TN[i],
    "/",
    metrics$FP[i],
    "/",
    metrics$FN[i],
    "\n"
  )

  cat(
    "Accuracy:",
    sprintf(
      "%.4f",
      metrics$Accuracy[i]
    ),
    "\n"
  )

  cat(
    "Sensitivity:",
    sprintf(
      "%.4f",
      metrics$Sensitivity[i]
    ),
    "\n"
  )

  cat(
    "Specificity:",
    sprintf(
      "%.4f",
      metrics$Specificity[i]
    ),
    "\n"
  )

  cat(
    "Balanced Accuracy:",
    sprintf(
      "%.4f",
      metrics$Balanced_Accuracy[i]
    ),
    "\n"
  )

  cat(
    "AUC:",
    sprintf(
      "%.4f",
      metrics$AUC[i]
    ),
    "\n"
  )

  cat(
    "AUC 95% CI:",
    sprintf(
      "%.4f",
      metrics$AUC_CI95_lower[i]
    ),
    "to",
    sprintf(
      "%.4f",
      metrics$AUC_CI95_upper[i]
    ),
    "\n"
  )

  cat(
    "Brier score:",
    sprintf(
      "%.4f",
      metrics$Brier_Score[i]
    ),
    "\n"
  )
}

cat("\nIMPORTANT:\n")

cat(
  "No model component was changed after external labels were opened.\n"
)

cat(
  "Dataset-specific external results are the primary results.\n"
)

cat(
  "The pooled result is secondary because the two external cohorts differ in prevalence and cohort composition.\n"
)

cat(
  "External results must be reported as observed, regardless of performance.\n"
)

cat("\nOutputs:\n")

cat(
  "- results/week4_external_validation_predictions_strict.csv\n"
)

cat(
  "- results/week4_external_validation_metrics_strict.csv\n"
)

cat(
  "- results/week4_external_validation_audit_strict.csv\n"
)

cat(
  "- figures/week4_external_validation_roc_strict.png\n"
)
