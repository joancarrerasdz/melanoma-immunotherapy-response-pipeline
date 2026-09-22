#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================================\n")
cat("WEEK 4 — LOCKED EXTERNAL-VALIDATION REPRODUCTION\n")
cat("============================================================\n")

# ============================================================
# IMPORTANT METHODOLOGICAL BOUNDARY
#
# scripts/16_freeze_deployment_model_strict.R is deliberately
# NOT executed here.
#
# The deployment model and blind external predictions were
# frozen before external labels were opened. Week 4 reproduction
# verifies those frozen artifacts cryptographically and then
# reproduces the locked evaluation.
# ============================================================

required_scripts <- c(
  "scripts/14_external_input_contract_strict.R",
  "scripts/15_external_representation_audit_strict.R",
  "scripts/17_external_validation_locked_strict.R"
)

missing_scripts <- required_scripts[
  !file.exists(required_scripts)
]

if (length(missing_scripts) > 0L) {
  stop(
    paste(
      "Missing Week 4 scripts:",
      paste(missing_scripts, collapse = "\n")
    ),
    call. = FALSE
  )
}

required_frozen_files <- c(
  "results/week4_deployment_freeze_manifest_strict.csv",
  "results/week4_deployment_model_spec_strict.csv",
  "results/week4_deployment_rf_model_strict.rds",
  "results/week4_deployment_oof_predictions_strict.csv",
  "results/week4_external_predictions_blind_strict.csv",
  "results/week4_external_predictions_blind_manifest_strict.csv",
  "results/week4_expression_representation_decision_strict.csv"
)

missing_frozen <- required_frozen_files[
  !file.exists(required_frozen_files)
]

if (length(missing_frozen) > 0L) {
  stop(
    paste(
      "Missing frozen Week 4 artifacts:",
      paste(missing_frozen, collapse = "\n")
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
    args = c(
      "--vanilla",
      script
    )
  )

  if (!identical(status, 0L)) {
    stop(
      paste(
        "Week 4 step failed:",
        script
      ),
      call. = FALSE
    )
  }

  cat(
    "[PASS]",
    step,
    "\n"
  )
}

parameter_value <- function(
  df,
  parameter
) {

  value <- df$value[
    df$parameter == parameter
  ]

  if (length(value) != 1L) {
    stop(
      paste(
        "Expected exactly one value for parameter:",
        parameter
      ),
      call. = FALSE
    )
  }

  as.character(value)
}

audit_value <- function(
  df,
  item
) {

  value <- df$value[
    df$item == item
  ]

  if (length(value) != 1L) {
    stop(
      paste(
        "Expected exactly one audit value for:",
        item
      ),
      call. = FALSE
    )
  }

  as.character(value)
}

assert_close <- function(
  observed,
  expected,
  tolerance = 1e-8,
  label = "metric"
) {

  if (
    length(observed) != 1L ||
    is.na(observed) ||
    abs(observed - expected) > tolerance
  ) {
    stop(
      paste(
        "Unexpected",
        label,
        ": observed",
        observed,
        "expected",
        expected
      ),
      call. = FALSE
    )
  }
}

# ============================================================
# STEP 1 — reconstruct external inputs
# ============================================================

run_step(
  "STEP 1/3 — Strict external-input contract",
  "scripts/14_external_input_contract_strict.R"
)

# ============================================================
# STEP 2 — reproduce representation audit
# ============================================================

run_step(
  "STEP 2/3 — External expression-representation audit",
  "scripts/15_external_representation_audit_strict.R"
)

# ============================================================
# Verify external-input contract
# ============================================================

input_contract <- read.csv(
  "results/week4_external_input_contract_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(input_contract) != 2L) {
  stop(
    "External-input contract must contain exactly two datasets.",
    call. = FALSE
  )
}

get_dataset_row <- function(
  df,
  dataset
) {

  x <- df[
    df$dataset == dataset,
    ,
    drop = FALSE
  ]

  if (nrow(x) != 1L) {
    stop(
      paste(
        "Expected exactly one row for",
        dataset
      ),
      call. = FALSE
    )
  }

  x
}

contract_91061 <- get_dataset_row(
  input_contract,
  "GSE91061"
)

contract_78220 <- get_dataset_row(
  input_contract,
  "GSE78220"
)

stopifnot(
  contract_91061$validated_samples == 49L,
  contract_91061$matched_samples == 49L,
  contract_91061$required_core_genes == 12L,
  contract_91061$available_core_genes == 12L,
  contract_91061$duplicated_expression_sample_keys == 0L,
  contract_91061$duplicated_metadata_sample_keys == 0L,
  contract_91061$missing_sample_matches == 0L,
  !as.logical(
    contract_91061$response_used_in_expression_preparation
  ),

  contract_78220$validated_samples == 27L,
  contract_78220$matched_samples == 27L,
  contract_78220$required_core_genes == 12L,
  contract_78220$available_core_genes == 12L,
  contract_78220$duplicated_expression_sample_keys == 0L,
  contract_78220$duplicated_metadata_sample_keys == 0L,
  contract_78220$missing_sample_matches == 0L,
  !as.logical(
    contract_78220$response_used_in_expression_preparation
  )
)

# ============================================================
# Verify frozen representation decision
# ============================================================

representation_decision <- read.csv(
  "results/week4_expression_representation_decision_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  parameter_value(
    representation_decision,
    "gene_space"
  ) == "strict_week3_core_12",

  parameter_value(
    representation_decision,
    "n_genes"
  ) == "12",

  parameter_value(
    representation_decision,
    "representation"
  ) == "samplewise_rank",

  parameter_value(
    representation_decision,
    "ranking_scope"
  ) == "within_sample_across_frozen_12_gene_core",

  parameter_value(
    representation_decision,
    "external_labels_used"
  ) == "FALSE",

  parameter_value(
    representation_decision,
    "external_performance_calculated"
  ) == "FALSE",

  parameter_value(
    representation_decision,
    "week3_auc_reused_as_deployment_performance"
  ) == "FALSE",

  parameter_value(
    representation_decision,
    "gate"
  ) == "PASS"
)

# ============================================================
# Verify Gate 4 freeze integrity
# ============================================================

freeze_manifest <- read.csv(
  "results/week4_deployment_freeze_manifest_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(
  c(
    "file",
    "md5"
  ) %in% colnames(freeze_manifest)
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

  missing <- freeze_manifest$file[
    !file.exists(
      freeze_manifest$file
    )
  ]

  stop(
    paste(
      "Frozen files missing:",
      paste(
        missing,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

current_md5 <- unname(
  tools::md5sum(
    freeze_manifest$file
  )
)

md5_match <- (
  as.character(current_md5) ==
  as.character(freeze_manifest$md5)
)

if (!all(md5_match)) {

  md5_audit <- data.frame(
    file = freeze_manifest$file,
    frozen_md5 = freeze_manifest$md5,
    current_md5 = current_md5,
    match = md5_match,
    stringsAsFactors = FALSE
  )

  print(md5_audit)

  stop(
    "Gate 4 frozen artifacts changed during Week 4 reproduction.",
    call. = FALSE
  )
}

cat("\n")
cat("[PASS] Gate 4 frozen-artifact MD5 verification\n")

# ============================================================
# Verify frozen model specification
# ============================================================

model_spec <- read.csv(
  "results/week4_deployment_model_spec_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  parameter_value(
    model_spec,
    "training_dataset"
  ) == "GSE160638",

  parameter_value(
    model_spec,
    "training_samples"
  ) == "36",

  parameter_value(
    model_spec,
    "gene_set"
  ) == "strict_week3_core_5of5",

  parameter_value(
    model_spec,
    "n_genes"
  ) == "12",

  parameter_value(
    model_spec,
    "representation"
  ) == "samplewise_rank",

  parameter_value(
    model_spec,
    "selected_mtry"
  ) == "1",

  parameter_value(
    model_spec,
    "ntree"
  ) == "500",

  parameter_value(
    model_spec,
    "frozen_threshold"
  ) == "0.48",

  parameter_value(
    model_spec,
    "external_labels_used"
  ) == "FALSE",

  parameter_value(
    model_spec,
    "external_performance_calculated"
  ) == "FALSE"
)

# ============================================================
# Verify blind external predictions
# ============================================================

blind_predictions <- read.csv(
  "results/week4_external_predictions_blind_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(blind_predictions) != 76L) {
  stop(
    "Expected exactly 76 frozen blind predictions.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    paste(
      blind_predictions$dataset,
      blind_predictions$sample_id,
      sep = "::"
    )
  )
) {
  stop(
    "Duplicated frozen external prediction keys.",
    call. = FALSE
  )
}

if (
  sum(
    blind_predictions$dataset ==
    "GSE91061"
  ) != 49L ||
  sum(
    blind_predictions$dataset ==
    "GSE78220"
  ) != 27L
) {
  stop(
    "Frozen external prediction dataset counts differ from Gate 4.",
    call. = FALSE
  )
}

if (
  anyNA(
    blind_predictions$prob_responder
  ) ||
  any(
    blind_predictions$prob_responder < 0 |
    blind_predictions$prob_responder > 1
  )
) {
  stop(
    "Invalid frozen external probabilities.",
    call. = FALSE
  )
}

blind_manifest <- read.csv(
  "results/week4_external_predictions_blind_manifest_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest_91061 <- get_dataset_row(
  blind_manifest,
  "GSE91061"
)

manifest_78220 <- get_dataset_row(
  blind_manifest,
  "GSE78220"
)

stopifnot(
  manifest_91061$samples == 49L,
  manifest_91061$genes == 12L,
  manifest_91061$representation == "samplewise_rank",
  !as.logical(
    manifest_91061$external_labels_read
  ),
  !as.logical(
    manifest_91061$external_performance_calculated
  ),

  manifest_78220$samples == 27L,
  manifest_78220$genes == 12L,
  manifest_78220$representation == "samplewise_rank",
  !as.logical(
    manifest_78220$external_labels_read
  ),
  !as.logical(
    manifest_78220$external_performance_calculated
  )
)

freeze_commit <- trimws(
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

if (
  length(freeze_commit) != 1L ||
  !nzchar(freeze_commit)
) {
  stop(
    "Could not identify blind-prediction freeze commit.",
    call. = FALSE
  )
}

cat(
  "[PASS] Frozen blind predictions verified\n"
)

cat(
  "Frozen blind-prediction commit:",
  freeze_commit,
  "\n"
)

# ============================================================
# STEP 3 — locked evaluation
# ============================================================

run_step(
  "STEP 3/3 — Locked external validation",
  "scripts/17_external_validation_locked_strict.R"
)

# ============================================================
# Expected Week 4 outputs
# ============================================================

expected_outputs <- c(
  "data_processed/week4_GSE91061_core_expression_untransformed_strict.csv",
  "data_processed/week4_GSE78220_core_expression_untransformed_strict.csv",

  "results/week4_external_input_contract_strict.csv",
  "results/week4_external_expression_scale_audit_strict.csv",

  "results/week4_expression_representation_global_audit_strict.csv",
  "results/week4_expression_representation_gene_audit_strict.csv",
  "results/week4_expression_representation_similarity_strict.csv",
  "results/week4_expression_representation_decision_strict.csv",

  "results/week4_deployment_model_spec_strict.csv",
  "results/week4_deployment_oof_predictions_strict.csv",
  "results/week4_deployment_threshold_search_strict.csv",
  "results/week4_deployment_tuning_strict.csv",
  "results/week4_deployment_cv_folds_strict.csv",
  "results/week4_deployment_environment_strict.csv",
  "results/week4_deployment_freeze_manifest_strict.csv",
  "results/week4_deployment_rf_model_strict.rds",

  "results/week4_external_predictions_blind_strict.csv",
  "results/week4_external_predictions_blind_manifest_strict.csv",

  "results/week4_external_validation_predictions_strict.csv",
  "results/week4_external_validation_metrics_strict.csv",
  "results/week4_external_validation_audit_strict.csv",
  "figures/week4_external_validation_roc_strict.png"
)

missing_outputs <- expected_outputs[
  !file.exists(expected_outputs)
]

if (length(missing_outputs) > 0L) {
  stop(
    paste(
      "Expected Week 4 outputs are missing:",
      paste(
        missing_outputs,
        collapse = "\n"
      )
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
      "Week 4 contains empty output files:",
      paste(
        empty_outputs,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

# ============================================================
# Validate external results exactly as observed
# ============================================================

metrics <- read.csv(
  "results/week4_external_validation_metrics_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metric_row <- function(scope) {

  x <- metrics[
    metrics$scope == scope,
    ,
    drop = FALSE
  ]

  if (nrow(x) != 1L) {
    stop(
      paste(
        "Expected exactly one metric row for:",
        scope
      ),
      call. = FALSE
    )
  }

  x
}

gse91061 <- metric_row(
  "GSE91061"
)

gse78220 <- metric_row(
  "GSE78220"
)

pooled <- metric_row(
  "POOLED"
)

stopifnot(
  gse91061$analysis_role ==
    "primary_external_cohort",
  gse91061$samples == 49L,
  gse91061$responders == 10L,
  gse91061$nonresponders == 39L,
  gse91061$TP == 9L,
  gse91061$TN == 7L,
  gse91061$FP == 32L,
  gse91061$FN == 1L,

  gse78220$analysis_role ==
    "primary_external_cohort",
  gse78220$samples == 27L,
  gse78220$responders == 15L,
  gse78220$nonresponders == 12L,
  gse78220$TP == 12L,
  gse78220$TN == 6L,
  gse78220$FP == 6L,
  gse78220$FN == 3L,

  pooled$analysis_role ==
    "secondary_pooled_summary",
  pooled$samples == 76L,
  pooled$responders == 25L,
  pooled$nonresponders == 51L,
  pooled$TP == 21L,
  pooled$TN == 13L,
  pooled$FP == 38L,
  pooled$FN == 4L
)

assert_close(
  gse91061$threshold,
  0.48,
  label = "GSE91061 threshold"
)

assert_close(
  gse91061$Accuracy,
  0.326530612244898,
  label = "GSE91061 Accuracy"
)

assert_close(
  gse91061$Sensitivity,
  0.9,
  label = "GSE91061 Sensitivity"
)

assert_close(
  gse91061$Specificity,
  0.179487179487179,
  label = "GSE91061 Specificity"
)

assert_close(
  gse91061$Balanced_Accuracy,
  0.53974358974359,
  label = "GSE91061 Balanced Accuracy"
)

assert_close(
  gse91061$AUC,
  0.664102564102564,
  label = "GSE91061 AUC"
)

assert_close(
  gse78220$threshold,
  0.48,
  label = "GSE78220 threshold"
)

assert_close(
  gse78220$Accuracy,
  0.666666666666667,
  label = "GSE78220 Accuracy"
)

assert_close(
  gse78220$Sensitivity,
  0.8,
  label = "GSE78220 Sensitivity"
)

assert_close(
  gse78220$Specificity,
  0.5,
  label = "GSE78220 Specificity"
)

assert_close(
  gse78220$Balanced_Accuracy,
  0.65,
  label = "GSE78220 Balanced Accuracy"
)

assert_close(
  gse78220$AUC,
  0.644444444444444,
  label = "GSE78220 AUC"
)

assert_close(
  pooled$threshold,
  0.48,
  label = "Pooled threshold"
)

assert_close(
  pooled$Accuracy,
  0.447368421052632,
  label = "Pooled Accuracy"
)

assert_close(
  pooled$Sensitivity,
  0.84,
  label = "Pooled Sensitivity"
)

assert_close(
  pooled$Specificity,
  0.254901960784314,
  label = "Pooled Specificity"
)

assert_close(
  pooled$Balanced_Accuracy,
  0.547450980392157,
  label = "Pooled Balanced Accuracy"
)

assert_close(
  pooled$AUC,
  0.571764705882353,
  label = "Pooled AUC"
)

# ============================================================
# Validate locked-evaluation audit
# ============================================================

evaluation_audit <- read.csv(
  "results/week4_external_validation_audit_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  audit_value(
    evaluation_audit,
    "gate4_freeze_md5_verified"
  ) == "TRUE",

  audit_value(
    evaluation_audit,
    "blind_prediction_freeze_commit"
  ) == freeze_commit,

  audit_value(
    evaluation_audit,
    "blind_predictions"
  ) == "76",

  audit_value(
    evaluation_audit,
    "external_samples_joined"
  ) == "76",

  audit_value(
    evaluation_audit,
    "genes"
  ) == "12",

  audit_value(
    evaluation_audit,
    "representation"
  ) == "samplewise_rank",

  audit_value(
    evaluation_audit,
    "selected_mtry"
  ) == "1",

  audit_value(
    evaluation_audit,
    "ntree"
  ) == "500",

  audit_value(
    evaluation_audit,
    "frozen_threshold"
  ) == "0.48",

  audit_value(
    evaluation_audit,
    "model_refitted_after_labels_opened"
  ) == "FALSE",

  audit_value(
    evaluation_audit,
    "genes_reselected_after_labels_opened"
  ) == "FALSE",

  audit_value(
    evaluation_audit,
    "representation_changed_after_labels_opened"
  ) == "FALSE",

  audit_value(
    evaluation_audit,
    "mtry_changed_after_labels_opened"
  ) == "FALSE",

  audit_value(
    evaluation_audit,
    "threshold_changed_after_labels_opened"
  ) == "FALSE",

  audit_value(
    evaluation_audit,
    "threshold_optimized_on_external_labels"
  ) == "FALSE",

  audit_value(
    evaluation_audit,
    "primary_reporting"
  ) == "cohort_specific",

  audit_value(
    evaluation_audit,
    "pooled_reporting"
  ) == "secondary"
)

# ============================================================
# Final gate
# ============================================================

cat("\n")
cat("============================================================\n")
cat("WEEK 4 REPRODUCTION: PASS\n")
cat("============================================================\n")

cat("Frozen Week 3 core genes: 12\n")
cat("Frozen representation: samplewise_rank\n")
cat("Frozen mtry: 1\n")
cat("Frozen ntree: 500\n")
cat("Frozen threshold: 0.48\n")
cat("Blind external predictions: 76\n")

cat("\nPrimary external results:\n")

cat(
  "- GSE91061: n=49 | AUC=0.6641 | ",
  "Accuracy=0.3265 | Sensitivity=0.9000 | ",
  "Specificity=0.1795 | Balanced Accuracy=0.5397\n",
  sep = ""
)

cat(
  "- GSE78220: n=27 | AUC=0.6444 | ",
  "Accuracy=0.6667 | Sensitivity=0.8000 | ",
  "Specificity=0.5000 | Balanced Accuracy=0.6500\n",
  sep = ""
)

cat("\nSecondary pooled summary:\n")

cat(
  "- n=76 | AUC=0.5718 | Accuracy=0.4474 | ",
  "Sensitivity=0.8400 | Specificity=0.2549 | ",
  "Balanced Accuracy=0.5475\n",
  sep = ""
)

cat("\nReproducibility safeguards:\n")
cat("- Gate 4 frozen-artifact MD5 integrity: PASS\n")
cat("- Frozen model regenerated after unblinding: FALSE\n")
cat("- Blind external predictions regenerated after unblinding: FALSE\n")
cat("- External threshold optimization: FALSE\n")
cat("- Dataset-specific external reporting: PRIMARY\n")
cat("- Pooled reporting: SECONDARY\n")

cat("\nIMPORTANT:\n")

cat(
  "scripts/16_freeze_deployment_model_strict.R was NOT executed ",
  "during final Week 4 reproduction.\n",
  sep = ""
)

cat(
  "The deployment model and blind predictions remain those frozen ",
  "before external labels were opened.\n",
  sep = ""
)

cat(
  "Observed external performance is reproduced without post hoc ",
  "model modification.\n",
  sep = ""
)

cat(
  "Week 4 external-validation results must be reported as observed.\n"
)
