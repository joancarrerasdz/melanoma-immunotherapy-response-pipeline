#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
  library(caret)
  library(randomForest)
  library(pROC)
})

# ============================================================
# Week 4 — frozen deployment model
#
# PURPOSE
# -------
# Freeze a deployment classifier using:
#
#   - GSE160638 training cohort only
#   - frozen strict Week 3 12-gene core
#   - frozen Week 4 sample-wise-rank representation
#
# External response labels are NEVER read by this script.
#
# External predictions are generated blindly and frozen before
# any independent external-performance evaluation.
# ============================================================

results_dir <- "results"
processed_dir <- "data_processed"

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# Reproducibility constants
# ------------------------------------------------------------

cv_split_seed <- 2026092201L
cv_model_seed <- 2026092202L
final_model_seed <- 2026092203L

n_internal_folds <- 5L
ntree_final <- 500L

# ============================================================
# 1. Helper functions
# ============================================================

clean_gene_ids <- function(x) {

  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)
  x <- gsub("\\s+", "", x)

  toupper(x)
}

extract_core <- function(
  x,
  core_genes,
  dataset
) {

  if (is.null(rownames(x))) {
    stop(
      dataset,
      ": expression matrix has no gene row names.",
      call. = FALSE
    )
  }

  gene_keys <- clean_gene_ids(
    rownames(x)
  )

  duplicated_core <- intersect(
    core_genes,
    unique(
      gene_keys[
        duplicated(gene_keys)
      ]
    )
  )

  if (length(duplicated_core) > 0L) {
    stop(
      dataset,
      ": duplicated strict-core genes: ",
      paste(
        duplicated_core,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  idx <- match(
    core_genes,
    gene_keys
  )

  if (anyNA(idx)) {

    missing_genes <- core_genes[
      is.na(idx)
    ]

    stop(
      dataset,
      ": missing strict-core genes: ",
      paste(
        missing_genes,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  out <- x[
    idx,
    ,
    drop = FALSE
  ]

  rownames(out) <- core_genes

  out
}

read_external_core <- function(
  path,
  core_genes,
  dataset
) {

  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x) < 2L) {
    stop(
      dataset,
      ": invalid expression file.",
      call. = FALSE
    )
  }

  genes <- clean_gene_ids(
    x[[1]]
  )

  x[[1]] <- NULL

  mat <- as.matrix(x)
  mode(mat) <- "numeric"

  rownames(mat) <- genes

  if (
    anyNA(mat) ||
    any(!is.finite(mat)) ||
    any(mat < 0)
  ) {
    stop(
      dataset,
      ": invalid external expression values.",
      call. = FALSE
    )
  }

  extract_core(
    mat,
    core_genes,
    dataset
  )
}

samplewise_rank <- function(x) {

  n_genes <- nrow(x)

  out <- vapply(
    seq_len(ncol(x)),
    function(j) {

      r <- rank(
        x[, j],
        ties.method = "average"
      )

      (r - 0.5) / n_genes
    },
    numeric(n_genes)
  )

  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)

  out
}

binary_metrics <- function(
  observed,
  probability,
  threshold
) {

  observed <- factor(
    observed,
    levels = c(
      "NonResponder",
      "Responder"
    )
  )

  predicted <- factor(
    ifelse(
      probability >= threshold,
      "Responder",
      "NonResponder"
    ),
    levels = levels(observed)
  )

  cm <- confusionMatrix(
    data = predicted,
    reference = observed,
    positive = "Responder"
  )

  sensitivity <- unname(
    cm$byClass["Sensitivity"]
  )

  specificity <- unname(
    cm$byClass["Specificity"]
  )

  data.frame(
    threshold = threshold,
    Accuracy = unname(
      cm$overall["Accuracy"]
    ),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Balanced_Accuracy = mean(
      c(
        sensitivity,
        specificity
      )
    ),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 2. Frozen Week 3 gene core
# ============================================================

core_file <- file.path(
  results_dir,
  "week3_consensus_core_5of5_strict.csv"
)

core <- read.csv(
  core_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"gene" %in% colnames(core)) {
  stop(
    "Column 'gene' missing from strict Week 3 core file.",
    call. = FALSE
  )
}

core_genes <- unique(
  clean_gene_ids(
    core$gene
  )
)

if (length(core_genes) != 12L) {
  stop(
    "Expected exactly 12 frozen strict-core genes; found ",
    length(core_genes),
    ".",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("WEEK 4 FROZEN DEPLOYMENT MODEL\n")
cat("============================================================\n")

cat(
  "Frozen genes:",
  length(core_genes),
  "\n"
)

cat(
  paste(
    core_genes,
    collapse = ", "
  ),
  "\n"
)

# ============================================================
# 3. Training cohort — GSE160638 only
# ============================================================

counts_object <- readRDS(
  file.path(
    processed_dir,
    "GSE160638_raw_counts_PD1_aligned.rds"
  )
)

if (inherits(counts_object, "DGEList")) {

  counts <- counts_object$counts

} else if (
  is.matrix(counts_object) ||
  is.data.frame(counts_object)
) {

  counts <- as.matrix(
    counts_object
  )

} else {

  stop(
    "Unsupported GSE160638 counts object.",
    call. = FALSE
  )
}

mode(counts) <- "numeric"

if (
  anyNA(counts) ||
  any(!is.finite(counts)) ||
  any(counts < 0)
) {
  stop(
    "Invalid values in GSE160638 raw counts.",
    call. = FALSE
  )
}

if (ncol(counts) != 36L) {
  stop(
    "Expected 36 GSE160638 training samples; found ",
    ncol(counts),
    ".",
    call. = FALSE
  )
}

clinical <- read.csv(
  file.path(
    processed_dir,
    "GSE160638_TableS2A_PD1_clinical.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_clinical_columns <- c(
  "patient_id",
  "response"
)

missing_clinical_columns <- setdiff(
  required_clinical_columns,
  colnames(clinical)
)

if (length(missing_clinical_columns) > 0L) {
  stop(
    "Missing clinical columns: ",
    paste(
      missing_clinical_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

clinical$patient_id <- as.character(
  clinical$patient_id
)

idx <- match(
  colnames(counts),
  clinical$patient_id
)

if (anyNA(idx)) {

  missing_samples <- colnames(counts)[
    is.na(idx)
  ]

  stop(
    "Training samples missing from clinical metadata: ",
    paste(
      missing_samples,
      collapse = ", "
    ),
    call. = FALSE
  )
}

clinical <- clinical[
  idx,
  ,
  drop = FALSE
]

if (!identical(
  colnames(counts),
  clinical$patient_id
)) {
  stop(
    "Training clinical metadata is not aligned to counts.",
    call. = FALSE
  )
}

y <- factor(
  clinical$response,
  levels = c(
    "NonResponder",
    "Responder"
  )
)

if (anyNA(y)) {
  stop(
    "Unexpected or missing training response labels.",
    call. = FALSE
  )
}

response_counts <- table(y)

if (
  unname(response_counts["NonResponder"]) != 14L ||
  unname(response_counts["Responder"]) != 22L
) {
  stop(
    "Training response distribution differs from validated Week 1 cohort.",
    call. = FALSE
  )
}

cat(
  "Training samples:",
  ncol(counts),
  "\n"
)

cat(
  "NonResponders:",
  unname(response_counts["NonResponder"]),
  "\n"
)

cat(
  "Responders:",
  unname(response_counts["Responder"]),
  "\n"
)

# ============================================================
# 4. Frozen sample-wise-rank representation
# ============================================================

dge <- DGEList(
  counts = counts
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

training_log <- cpm(
  dge,
  log = TRUE,
  prior.count = 1
)

training_log <- extract_core(
  training_log,
  core_genes,
  "GSE160638"
)

training_rank <- samplewise_rank(
  training_log
)

if (
  nrow(training_rank) != 12L ||
  ncol(training_rank) != 36L
) {
  stop(
    "Unexpected training rank-matrix dimensions.",
    call. = FALSE
  )
}

x_train <- as.data.frame(
  t(training_rank),
  check.names = FALSE
)

colnames(x_train) <- make.names(
  colnames(x_train),
  unique = TRUE
)

train_df <- x_train
train_df$response <- y

# ============================================================
# 5. Fixed internal development CV
# ============================================================

set.seed(
  cv_split_seed
)

development_index <- createFolds(
  y,
  k = n_internal_folds,
  returnTrain = TRUE
)

if (
  length(development_index) !=
  n_internal_folds
) {
  stop(
    "Unexpected number of development folds.",
    call. = FALSE
  )
}

held_out_indices <- lapply(
  development_index,
  function(train_idx) {
    setdiff(
      seq_along(y),
      train_idx
    )
  }
)

held_out_frequency <- table(
  unlist(
    held_out_indices
  )
)

if (
  length(held_out_frequency) != length(y) ||
  any(held_out_frequency != 1L)
) {
  stop(
    "Development CV does not hold out every sample exactly once.",
    call. = FALSE
  )
}

fold_assignment <- integer(
  length(y)
)

for (
  i in seq_along(
    held_out_indices
  )
) {

  fold_assignment[
    held_out_indices[[i]]
  ] <- i
}

fold_manifest <- data.frame(
  sample_id = colnames(counts),
  response = as.character(y),
  development_fold = fold_assignment,
  stringsAsFactors = FALSE
)

write.csv(
  fold_manifest,
  file.path(
    results_dir,
    "week4_deployment_cv_folds_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. Random-Forest mtry development
# ============================================================

mtry_grid <- data.frame(
  mtry = seq_len(
    ncol(x_train)
  )
)

ctrl <- trainControl(
  method = "cv",
  number = n_internal_folds,
  index = development_index,
  classProbs = TRUE,
  savePredictions = "final",
  returnResamp = "all",
  allowParallel = FALSE
)

set.seed(
  cv_model_seed
)

rf_tuned <- train(
  response ~ .,
  data = train_df,
  method = "rf",
  metric = "Accuracy",
  trControl = ctrl,
  tuneGrid = mtry_grid,
  ntree = ntree_final,
  importance = TRUE
)

best_mtry <- as.integer(
  rf_tuned$bestTune$mtry
)

cat("\n============================================================\n")
cat("DEVELOPMENT TUNING\n")
cat("============================================================\n")

cat(
  "Candidate mtry:",
  paste(
    mtry_grid$mtry,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Selected mtry:",
  best_mtry,
  "\n"
)

write.csv(
  rf_tuned$results,
  file.path(
    results_dir,
    "week4_deployment_tuning_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. Development out-of-fold probabilities
# ============================================================

oof <- rf_tuned$pred

if ("mtry" %in% colnames(oof)) {

  oof <- oof[
    oof$mtry == best_mtry,
    ,
    drop = FALSE
  ]
}

oof <- oof[
  order(oof$rowIndex),
  ,
  drop = FALSE
]

if (
  nrow(oof) != 36L ||
  !identical(
    as.integer(oof$rowIndex),
    seq_len(36L)
  )
) {
  stop(
    "Expected exactly one OOF prediction for each of 36 training samples.",
    call. = FALSE
  )
}

if (!"Responder" %in% colnames(oof)) {
  stop(
    "Responder probabilities missing from development OOF predictions.",
    call. = FALSE
  )
}

oof_probability <- as.numeric(
  oof$Responder
)

if (
  anyNA(oof_probability) ||
  any(!is.finite(oof_probability)) ||
  any(
    oof_probability < 0 |
    oof_probability > 1
  )
) {
  stop(
    "Invalid OOF probabilities.",
    call. = FALSE
  )
}

oof_observed <- factor(
  oof$obs,
  levels = c(
    "NonResponder",
    "Responder"
  )
)

oof_roc <- roc(
  response = oof_observed,
  predictor = oof_probability,
  levels = c(
    "NonResponder",
    "Responder"
  ),
  direction = "<",
  quiet = TRUE
)

development_oof_auc <- as.numeric(
  auc(oof_roc)
)

# ============================================================
# 8. Freeze probability threshold from training OOF only
# ============================================================

threshold_candidates <- seq(
  0.05,
  0.95,
  by = 0.01
)

threshold_search <- do.call(
  rbind,
  lapply(
    threshold_candidates,
    function(threshold) {

      binary_metrics(
        observed = oof_observed,
        probability = oof_probability,
        threshold = threshold
      )
    }
  )
)

selection_order <- order(
  -threshold_search$Balanced_Accuracy,
  abs(
    threshold_search$threshold -
    0.5
  ),
  threshold_search$threshold
)

threshold_search$selected <- FALSE
threshold_search$selected[
  selection_order[1]
] <- TRUE

selected_threshold <- threshold_search$threshold[
  selection_order[1]
]

selected_threshold_metrics <- threshold_search[
  selection_order[1],
  ,
  drop = FALSE
]

write.csv(
  threshold_search,
  file.path(
    results_dir,
    "week4_deployment_threshold_search_strict.csv"
  ),
  row.names = FALSE
)

cat(
  "Development OOF AUC:",
  sprintf(
    "%.4f",
    development_oof_auc
  ),
  "\n"
)

cat(
  "Selected probability threshold:",
  sprintf(
    "%.2f",
    selected_threshold
  ),
  "\n"
)

cat(
  "OOF balanced accuracy at frozen threshold:",
  sprintf(
    "%.4f",
    selected_threshold_metrics$Balanced_Accuracy
  ),
  "\n"
)

# ============================================================
# 9. Save development OOF predictions
# ============================================================

oof_predicted <- ifelse(
  oof_probability >= selected_threshold,
  "Responder",
  "NonResponder"
)

oof_output <- data.frame(
  sample_id = colnames(counts),
  observed = as.character(
    oof_observed
  ),
  prob_responder = oof_probability,
  predicted = oof_predicted,
  development_fold = fold_assignment,
  selected_mtry = best_mtry,
  frozen_threshold = selected_threshold,
  stringsAsFactors = FALSE
)

write.csv(
  oof_output,
  file.path(
    results_dir,
    "week4_deployment_oof_predictions_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10. Fit frozen final deployment model
# ============================================================

set.seed(
  final_model_seed
)

rf_final <- randomForest(
  x = x_train,
  y = y,
  ntree = ntree_final,
  mtry = best_mtry,
  importance = TRUE
)

saveRDS(
  rf_final,
  file.path(
    results_dir,
    "week4_deployment_rf_model_strict.rds"
  ),
  version = 3
)

# ============================================================
# 11. External expression — labels are NOT read
# ============================================================

gse91061_raw <- read_external_core(
  file.path(
    processed_dir,
    "week4_GSE91061_core_expression_untransformed_strict.csv"
  ),
  core_genes,
  "GSE91061"
)

gse78220_raw <- read_external_core(
  file.path(
    processed_dir,
    "week4_GSE78220_core_expression_untransformed_strict.csv"
  ),
  core_genes,
  "GSE78220"
)

gse91061_log <- log2(
  gse91061_raw + 1
)

gse78220_log <- log2(
  gse78220_raw + 1
)

gse91061_rank <- samplewise_rank(
  gse91061_log
)

gse78220_rank <- samplewise_rank(
  gse78220_log
)

x_gse91061 <- as.data.frame(
  t(gse91061_rank),
  check.names = FALSE
)

x_gse78220 <- as.data.frame(
  t(gse78220_rank),
  check.names = FALSE
)

colnames(x_gse91061) <- make.names(
  colnames(x_gse91061),
  unique = TRUE
)

colnames(x_gse78220) <- make.names(
  colnames(x_gse78220),
  unique = TRUE
)

if (
  !identical(
    colnames(x_train),
    colnames(x_gse91061)
  ) ||
  !identical(
    colnames(x_train),
    colnames(x_gse78220)
  )
) {
  stop(
    "Training/external deployment feature order differs.",
    call. = FALSE
  )
}

if (
  nrow(x_gse91061) != 49L ||
  nrow(x_gse78220) != 27L
) {
  stop(
    "Unexpected external sample counts.",
    call. = FALSE
  )
}

# ============================================================
# 12. Blind external predictions
# ============================================================

prob_gse91061 <- predict(
  rf_final,
  newdata = x_gse91061,
  type = "prob"
)[, "Responder"]

prob_gse78220 <- predict(
  rf_final,
  newdata = x_gse78220,
  type = "prob"
)[, "Responder"]

external_predictions <- rbind(
  data.frame(
    dataset = "GSE91061",
    sample_id = rownames(
      x_gse91061
    ),
    prob_responder = as.numeric(
      prob_gse91061
    ),
    predicted = ifelse(
      prob_gse91061 >=
        selected_threshold,
      "Responder",
      "NonResponder"
    ),
    frozen_threshold = selected_threshold,
    representation = "samplewise_rank",
    n_genes = 12L,
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "GSE78220",
    sample_id = rownames(
      x_gse78220
    ),
    prob_responder = as.numeric(
      prob_gse78220
    ),
    predicted = ifelse(
      prob_gse78220 >=
        selected_threshold,
      "Responder",
      "NonResponder"
    ),
    frozen_threshold = selected_threshold,
    representation = "samplewise_rank",
    n_genes = 12L,
    stringsAsFactors = FALSE
  )
)

if (
  anyNA(
    external_predictions$prob_responder
  ) ||
  any(
    external_predictions$prob_responder < 0 |
    external_predictions$prob_responder > 1
  )
) {
  stop(
    "Invalid external probabilities.",
    call. = FALSE
  )
}

forbidden_external_columns <- c(
  "response",
  "observed",
  "label",
  "outcome"
)

if (
  any(
    forbidden_external_columns %in%
    tolower(
      colnames(
        external_predictions
      )
    )
  )
) {
  stop(
    "External-label column detected in blind prediction output.",
    call. = FALSE
  )
}

write.csv(
  external_predictions,
  file.path(
    results_dir,
    "week4_external_predictions_blind_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. Blind prediction manifest
# ============================================================

external_manifest <- do.call(
  rbind,
  lapply(
    split(
      external_predictions,
      external_predictions$dataset
    ),
    function(x) {

      data.frame(
        dataset = unique(
          x$dataset
        ),
        samples = nrow(x),
        genes = unique(
          x$n_genes
        ),
        representation = unique(
          x$representation
        ),
        probability_min = min(
          x$prob_responder
        ),
        probability_max = max(
          x$prob_responder
        ),
        predicted_NonResponder = sum(
          x$predicted ==
          "NonResponder"
        ),
        predicted_Responder = sum(
          x$predicted ==
          "Responder"
        ),
        external_labels_read = FALSE,
        external_performance_calculated = FALSE,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  external_manifest,
  file.path(
    results_dir,
    "week4_external_predictions_blind_manifest_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 14. Frozen deployment specification
# ============================================================

model_spec <- data.frame(
  parameter = c(
    "training_dataset",
    "training_samples",
    "training_nonresponders",
    "training_responders",
    "gene_set",
    "n_genes",
    "representation",
    "rank_ties_method",
    "rank_formula",
    "development_cv",
    "tuning_metric",
    "mtry_candidates",
    "selected_mtry",
    "ntree",
    "threshold_selection",
    "frozen_threshold",
    "development_oof_auc",
    "development_oof_accuracy",
    "development_oof_sensitivity",
    "development_oof_specificity",
    "development_oof_balanced_accuracy",
    "cv_split_seed",
    "cv_model_seed",
    "final_model_seed",
    "external_labels_used",
    "external_performance_calculated"
  ),
  value = c(
    "GSE160638",
    "36",
    "14",
    "22",
    "strict_week3_core_5of5",
    "12",
    "samplewise_rank",
    "average",
    "(rank - 0.5) / 12",
    "5-fold_stratified",
    "Accuracy",
    paste(
      mtry_grid$mtry,
      collapse = ","
    ),
    as.character(
      best_mtry
    ),
    as.character(
      ntree_final
    ),
    "max_OOF_balanced_accuracy_ties_closest_to_0.5",
    sprintf(
      "%.2f",
      selected_threshold
    ),
    sprintf(
      "%.6f",
      development_oof_auc
    ),
    sprintf(
      "%.6f",
      selected_threshold_metrics$Accuracy
    ),
    sprintf(
      "%.6f",
      selected_threshold_metrics$Sensitivity
    ),
    sprintf(
      "%.6f",
      selected_threshold_metrics$Specificity
    ),
    sprintf(
      "%.6f",
      selected_threshold_metrics$Balanced_Accuracy
    ),
    as.character(
      cv_split_seed
    ),
    as.character(
      cv_model_seed
    ),
    as.character(
      final_model_seed
    ),
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  model_spec,
  file.path(
    results_dir,
    "week4_deployment_model_spec_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 15. Package/environment record
# ============================================================

package_versions <- data.frame(
  component = c(
    "R",
    "edgeR",
    "caret",
    "randomForest",
    "pROC"
  ),
  version = c(
    as.character(
      getRversion()
    ),
    as.character(
      packageVersion("edgeR")
    ),
    as.character(
      packageVersion("caret")
    ),
    as.character(
      packageVersion("randomForest")
    ),
    as.character(
      packageVersion("pROC")
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_versions,
  file.path(
    results_dir,
    "week4_deployment_environment_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 16. Freeze manifest
# ============================================================

freeze_files <- c(
  "results/week3_consensus_core_5of5_strict.csv",
  "data_processed/GSE160638_raw_counts_PD1_aligned.rds",
  "data_processed/week4_GSE91061_core_expression_untransformed_strict.csv",
  "data_processed/week4_GSE78220_core_expression_untransformed_strict.csv",
  "results/week4_deployment_model_spec_strict.csv",
  "results/week4_deployment_oof_predictions_strict.csv",
  "results/week4_external_predictions_blind_strict.csv",
  "results/week4_deployment_rf_model_strict.rds"
)

if (!all(
  file.exists(
    freeze_files
  )
)) {
  stop(
    "One or more freeze-manifest files are missing.",
    call. = FALSE
  )
}

freeze_manifest <- data.frame(
  file = freeze_files,
  md5 = unname(
    tools::md5sum(
      freeze_files
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  freeze_manifest,
  file.path(
    results_dir,
    "week4_deployment_freeze_manifest_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 17. Final validation assertions
# ============================================================

stopifnot(
  length(core_genes) == 12L,
  nrow(x_train) == 36L,
  nrow(x_gse91061) == 49L,
  nrow(x_gse78220) == 27L,
  nrow(oof_output) == 36L,
  nrow(external_predictions) == 76L,
  sum(
    external_predictions$dataset ==
    "GSE91061"
  ) == 49L,
  sum(
    external_predictions$dataset ==
    "GSE78220"
  ) == 27L,
  identical(
    colnames(x_train),
    colnames(x_gse91061)
  ),
  identical(
    colnames(x_train),
    colnames(x_gse78220)
  )
)

cat("\n============================================================\n")
cat("WEEK 4 FROZEN DEPLOYMENT MODEL: PASS\n")
cat("============================================================\n")

cat(
  "Training samples: 36\n"
)

cat(
  "Frozen genes: 12\n"
)

cat(
  "Representation: samplewise_rank\n"
)

cat(
  "Selected mtry:",
  best_mtry,
  "\n"
)

cat(
  "Frozen threshold:",
  sprintf(
    "%.2f",
    selected_threshold
  ),
  "\n"
)

cat(
  "Development OOF AUC:",
  sprintf(
    "%.4f",
    development_oof_auc
  ),
  "\n"
)

cat(
  "Development OOF balanced accuracy:",
  sprintf(
    "%.4f",
    selected_threshold_metrics$Balanced_Accuracy
  ),
  "\n"
)

cat(
  "Blind GSE91061 predictions:",
  sum(
    external_predictions$dataset ==
    "GSE91061"
  ),
  "\n"
)

cat(
  "Blind GSE78220 predictions:",
  sum(
    external_predictions$dataset ==
    "GSE78220"
  ),
  "\n"
)

cat(
  "External response labels used: FALSE\n"
)

cat(
  "External performance calculated: FALSE\n"
)

cat("\nIMPORTANT:\n")

cat(
  "Development OOF metrics are model-development metrics only.\n"
)

cat(
  "They are not independent external-validation estimates.\n"
)

cat(
  "The external predictions above are frozen before external-label evaluation.\n"
)

cat("\nOutputs:\n")

cat(
  "- results/week4_deployment_cv_folds_strict.csv\n"
)

cat(
  "- results/week4_deployment_tuning_strict.csv\n"
)

cat(
  "- results/week4_deployment_threshold_search_strict.csv\n"
)

cat(
  "- results/week4_deployment_oof_predictions_strict.csv\n"
)

cat(
  "- results/week4_deployment_rf_model_strict.rds\n"
)

cat(
  "- results/week4_deployment_model_spec_strict.csv\n"
)

cat(
  "- results/week4_deployment_environment_strict.csv\n"
)

cat(
  "- results/week4_external_predictions_blind_strict.csv\n"
)

cat(
  "- results/week4_external_predictions_blind_manifest_strict.csv\n"
)

cat(
  "- results/week4_deployment_freeze_manifest_strict.csv\n"
)
