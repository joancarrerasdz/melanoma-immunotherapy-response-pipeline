#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(randomForest)
  library(caret)
  library(pROC)
})

# ============================================================
# Week 3 — Strict nested cross-validation
# GSE160638 melanoma anti-PD-1 cohort
#
# Methodological rules:
#   1. Raw validated Week 1 counts are the only expression input.
#   2. Outer-test samples are never used for feature selection,
#      tuning, or model fitting.
#   3. Inner-validation samples are never used for feature
#      selection or model fitting.
#   4. Gene filtering and differential-expression ranking are
#      re-fitted inside each training partition.
#   5. Predictive expression transformation is sample-wise and
#      does not estimate parameters from held-out samples.
#   6. Final validation performance is calculated exclusively
#      from pooled outer out-of-fold predictions.
# ============================================================

set.seed(20260915)

base_dir <- normalizePath(".", mustWork = TRUE)

counts_file <- file.path(
  base_dir,
  "data_processed",
  "GSE160638_raw_counts_PD1_aligned.rds"
)

clinical_file <- file.path(
  base_dir,
  "data_processed",
  "GSE160638_TableS2A_PD1_clinical.csv"
)

results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

outer_k <- 5L
inner_k <- 3L
top_n_genes <- 100L
ntree <- 500L

# ------------------------------------------------------------
# 1. Load and validate Week 1 foundation
# ------------------------------------------------------------

if (!file.exists(counts_file)) {
  stop("Validated Week 1 count matrix is missing.", call. = FALSE)
}

if (!file.exists(clinical_file)) {
  stop("Validated Week 1 clinical table is missing.", call. = FALSE)
}

counts <- readRDS(counts_file)

clinical <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(counts)

if (!is.numeric(counts)) {
  stop("Count matrix is not numeric.", call. = FALSE)
}

if (!identical(dim(counts), c(17002L, 36L))) {
  stop(
    paste0(
      "Unexpected count dimensions: ",
      nrow(counts), " x ", ncol(counts),
      ". Expected 17002 x 36."
    ),
    call. = FALSE
  )
}

if (anyNA(counts)) {
  stop("Missing values detected in raw counts.", call. = FALSE)
}

if (any(counts < 0)) {
  stop("Negative raw counts detected.", call. = FALSE)
}

if (any(abs(counts - round(counts)) > 1e-8)) {
  stop("Non-integer-like raw counts detected.", call. = FALSE)
}

if (anyDuplicated(rownames(counts))) {
  stop("Duplicated gene identifiers detected.", call. = FALSE)
}

required_clinical <- c("patient_id", "response")

if (!all(required_clinical %in% colnames(clinical))) {
  stop(
    "Clinical table does not contain patient_id and response.",
    call. = FALSE
  )
}

if (!identical(colnames(counts), clinical$patient_id)) {
  stop(
    "Count matrix and clinical table sample order differ.",
    call. = FALSE
  )
}

y <- factor(
  clinical$response,
  levels = c("NonResponder", "Responder")
)

if (anyNA(y)) {
  stop("Unexpected response labels detected.", call. = FALSE)
}

expected_response_counts <- c(
  NonResponder = 14L,
  Responder = 22L
)

observed_response_counts <- table(y)

observed_response_counts <- setNames(
  as.integer(observed_response_counts),
  names(observed_response_counts)
)

if (
  !all(
    names(expected_response_counts) %in%
      names(observed_response_counts)
  ) ||
  !all(
    observed_response_counts[
      names(expected_response_counts)
    ] == expected_response_counts
  )
) {
  stop(
    "Response distribution differs from validated Week 1 cohort.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

samplewise_logcpm <- function(x, genes = NULL) {

  x <- as.matrix(x)

  if (!is.null(genes)) {
    x <- x[genes, , drop = FALSE]
  }

  lib_size <- colSums(x)

  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop(
      "Invalid library size encountered.",
      call. = FALSE
    )
  }

  cpm_matrix <- sweep(
    x,
    2,
    lib_size / 1e6,
    FUN = "/"
  )

  log2(cpm_matrix + 1)
}


select_genes_from_training <- function(
  counts_train,
  y_train,
  top_n = 100L
) {

  y_train <- factor(
    y_train,
    levels = c("NonResponder", "Responder")
  )

  if (length(unique(y_train)) != 2L) {
    stop(
      "Training partition does not contain both response classes.",
      call. = FALSE
    )
  }

  dge <- DGEList(
    counts = counts_train,
    group = y_train
  )

  keep <- filterByExpr(
    dge,
    group = y_train
  )

  if (sum(keep) < 2L) {
    stop(
      "Too few genes retained by filterByExpr.",
      call. = FALSE
    )
  }

  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- calcNormFactors(
    dge,
    method = "TMM"
  )

  design <- model.matrix(~ y_train)

  if (
    ncol(design) != 2L ||
    colnames(design)[2] != "y_trainResponder"
  ) {
    stop(
      "Unexpected differential-expression contrast.",
      call. = FALSE
    )
  }

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

  ranking <- topTags(
    qlf,
    n = Inf,
    sort.by = "PValue"
  )$table

  n_select <- min(
    as.integer(top_n),
    nrow(ranking)
  )

  selected <- rownames(ranking)[
    seq_len(n_select)
  ]

  list(
    genes = selected,
    n_after_filter = nrow(dge),
    ranking = ranking
  )
}


classification_metrics <- function(
  observed,
  predicted,
  prob_responder
) {

  observed <- factor(
    observed,
    levels = c("NonResponder", "Responder")
  )

  predicted <- factor(
    predicted,
    levels = c("NonResponder", "Responder")
  )

  cm <- table(
    predicted = predicted,
    observed = observed
  )

  tn <- cm["NonResponder", "NonResponder"]
  fn <- cm["NonResponder", "Responder"]
  fp <- cm["Responder", "NonResponder"]
  tp <- cm["Responder", "Responder"]

  accuracy <- (tp + tn) / sum(cm)

  sensitivity <- if ((tp + fn) > 0) {
    tp / (tp + fn)
  } else {
    NA_real_
  }

  specificity <- if ((tn + fp) > 0) {
    tn / (tn + fp)
  } else {
    NA_real_
  }

  balanced_accuracy <- mean(
    c(sensitivity, specificity),
    na.rm = TRUE
  )

  roc_obj <- roc(
    response = observed,
    predictor = prob_responder,
    levels = c("NonResponder", "Responder"),
    direction = "<",
    quiet = TRUE
  )

  data.frame(
    Accuracy = as.numeric(accuracy),
    Sensitivity = as.numeric(sensitivity),
    Specificity = as.numeric(specificity),
    Balanced_Accuracy = as.numeric(balanced_accuracy),
    AUC = as.numeric(auc(roc_obj)),
    stringsAsFactors = FALSE
  )
}


fit_predict_partition <- function(
  counts_train,
  y_train,
  counts_eval,
  mtry,
  top_n = 100L,
  seed = 1L
) {

  selection <- select_genes_from_training(
    counts_train = counts_train,
    y_train = y_train,
    top_n = top_n
  )

  selected_genes <- selection$genes

  train_expr <- samplewise_logcpm(
    counts_train,
    selected_genes
  )

  eval_expr <- samplewise_logcpm(
    counts_eval,
    selected_genes
  )

  x_train <- t(train_expr)
  x_eval <- t(eval_expr)

  if (!identical(
    colnames(x_train),
    colnames(x_eval)
  )) {
    stop(
      "Train/evaluation feature order differs.",
      call. = FALSE
    )
  }

  y_train <- factor(
    y_train,
    levels = c("NonResponder", "Responder")
  )

  fitted_mtry <- min(
    as.integer(mtry),
    ncol(x_train)
  )

  set.seed(seed)

  rf_fit <- randomForest(
    x = x_train,
    y = y_train,
    ntree = ntree,
    mtry = fitted_mtry,
    importance = TRUE
  )

  predicted <- predict(
    rf_fit,
    newdata = x_eval,
    type = "response"
  )

  probabilities <- predict(
    rf_fit,
    newdata = x_eval,
    type = "prob"
  )[, "Responder"]

  list(
    predicted = predicted,
    prob_responder = as.numeric(probabilities),
    genes = selected_genes,
    n_after_filter = selection$n_after_filter,
    fitted_mtry = fitted_mtry
  )
}

# ------------------------------------------------------------
# 3. Outer folds
# ------------------------------------------------------------

set.seed(20260915)

outer_folds <- createFolds(
  y,
  k = outer_k,
  list = TRUE,
  returnTrain = FALSE
)

if (length(outer_folds) != outer_k) {
  stop(
    "Unexpected number of outer folds.",
    call. = FALSE
  )
}

outer_fold_id <- integer(length(y))

for (i in seq_along(outer_folds)) {
  outer_fold_id[outer_folds[[i]]] <- i
}

if (any(outer_fold_id == 0L)) {
  stop(
    "At least one sample was not assigned to an outer fold.",
    call. = FALSE
  )
}

if (any(table(seq_along(y)) != 1L)) {
  stop(
    "Unexpected sample indexing state.",
    call. = FALSE
  )
}

# Fixed a priori.
# For a 100-gene model this corresponds approximately to:
# very small / sqrt(p) / p/4.
mtry_grid <- unique(
  pmax(
    1L,
    c(
      2L,
      floor(sqrt(top_n_genes)),
      floor(top_n_genes / 4)
    )
  )
)

cat("\n")
cat("============================================================\n")
cat("WEEK 3 STRICT NESTED-CV FOUNDATION\n")
cat("============================================================\n")
cat("Samples:", ncol(counts), "\n")
cat("Genes:", nrow(counts), "\n")
cat("Responders:", sum(y == "Responder"), "\n")
cat("NonResponders:", sum(y == "NonResponder"), "\n")
cat("Outer folds:", outer_k, "\n")
cat("Inner folds:", inner_k, "\n")
cat("Genes selected per training partition:", top_n_genes, "\n")
cat("mtry candidates:", paste(mtry_grid, collapse = ", "), "\n")
cat("\n")
cat("Input validation: PASS\n")
cat("Helper functions: READY\n")
cat("Outer-fold assignment: READY\n")
cat("\n")
cat("IMPORTANT:\n")
cat(
  "No Week 2 exploratory expression matrix is used as predictive input.\n"
)
cat(
  "All feature selection will be re-fitted inside training partitions.\n"
)
cat(
  "Final performance will be calculated from outer out-of-fold predictions only.\n"
)

# ------------------------------------------------------------
# 4. Strict nested cross-validation
# ------------------------------------------------------------

inner_cv_metrics <- data.frame()
inner_tuning_summary <- data.frame()
outer_predictions <- data.frame()
outer_fold_metrics <- data.frame()
outer_selected_genes <- data.frame()

for (outer_i in seq_along(outer_folds)) {

  cat("\n")
  cat("============================================================\n")
  cat("OUTER FOLD", outer_i, "/", outer_k, "\n")
  cat("============================================================\n")

  outer_test_idx <- outer_folds[[outer_i]]

  outer_train_idx <- setdiff(
    seq_len(ncol(counts)),
    outer_test_idx
  )

  counts_outer_train <- counts[
    ,
    outer_train_idx,
    drop = FALSE
  ]

  counts_outer_test <- counts[
    ,
    outer_test_idx,
    drop = FALSE
  ]

  y_outer_train <- factor(
    y[outer_train_idx],
    levels = c("NonResponder", "Responder")
  )

  y_outer_test <- factor(
    y[outer_test_idx],
    levels = c("NonResponder", "Responder")
  )

  if (length(unique(y_outer_train)) != 2L) {
    stop(
      paste0(
        "Outer fold ",
        outer_i,
        " training partition does not contain both classes."
      ),
      call. = FALSE
    )
  }

  if (length(unique(y_outer_test)) != 2L) {
    stop(
      paste0(
        "Outer fold ",
        outer_i,
        " test partition does not contain both classes."
      ),
      call. = FALSE
    )
  }

  cat(
    "Outer train:",
    length(outer_train_idx),
    "samples\n"
  )

  cat(
    "Outer test:",
    length(outer_test_idx),
    "samples\n"
  )

  cat(
    "Outer train distribution:",
    paste(
      names(table(y_outer_train)),
      as.integer(table(y_outer_train)),
      collapse = " | "
    ),
    "\n"
  )

  cat(
    "Outer test distribution:",
    paste(
      names(table(y_outer_test)),
      as.integer(table(y_outer_test)),
      collapse = " | "
    ),
    "\n"
  )

  # ----------------------------------------------------------
  # 4A. Inner CV: tune mtry
  # ----------------------------------------------------------

  set.seed(20260915 + outer_i)

  inner_folds <- createFolds(
    y_outer_train,
    k = inner_k,
    list = TRUE,
    returnTrain = FALSE
  )

  if (length(inner_folds) != inner_k) {
    stop(
      "Unexpected number of inner folds.",
      call. = FALSE
    )
  }

  outer_inner_metrics <- data.frame()

  for (inner_i in seq_along(inner_folds)) {

    cat(
      "\nInner fold",
      inner_i,
      "/",
      inner_k,
      "\n"
    )

    inner_val_local <- inner_folds[[inner_i]]

    inner_train_local <- setdiff(
      seq_len(ncol(counts_outer_train)),
      inner_val_local
    )

    counts_inner_train <- counts_outer_train[
      ,
      inner_train_local,
      drop = FALSE
    ]

    counts_inner_val <- counts_outer_train[
      ,
      inner_val_local,
      drop = FALSE
    ]

    y_inner_train <- factor(
      y_outer_train[inner_train_local],
      levels = c("NonResponder", "Responder")
    )

    y_inner_val <- factor(
      y_outer_train[inner_val_local],
      levels = c("NonResponder", "Responder")
    )

    if (length(unique(y_inner_train)) != 2L) {
      stop(
        "Inner training partition does not contain both classes.",
        call. = FALSE
      )
    }

    if (length(unique(y_inner_val)) != 2L) {
      stop(
        "Inner validation partition does not contain both classes.",
        call. = FALSE
      )
    }

    # Feature selection is fitted ONCE using only this
    # inner-training partition.
    inner_selection <- select_genes_from_training(
      counts_train = counts_inner_train,
      y_train = y_inner_train,
      top_n = top_n_genes
    )

    inner_genes <- inner_selection$genes

    inner_train_expr <- samplewise_logcpm(
      counts_inner_train,
      inner_genes
    )

    inner_val_expr <- samplewise_logcpm(
      counts_inner_val,
      inner_genes
    )

    X_inner_train <- t(inner_train_expr)
    X_inner_val <- t(inner_val_expr)

    if (!identical(
      colnames(X_inner_train),
      colnames(X_inner_val)
    )) {
      stop(
        "Inner train/validation feature order differs.",
        call. = FALSE
      )
    }

    for (mtry_i in seq_along(mtry_grid)) {

      candidate_mtry <- min(
        as.integer(mtry_grid[mtry_i]),
        ncol(X_inner_train)
      )

      set.seed(
        20261000 +
          outer_i * 1000 +
          inner_i * 100 +
          mtry_i
      )

      rf_inner <- randomForest(
        x = X_inner_train,
        y = y_inner_train,
        ntree = ntree,
        mtry = candidate_mtry,
        importance = FALSE
      )

      inner_predicted <- predict(
        rf_inner,
        newdata = X_inner_val,
        type = "response"
      )

      inner_probability <- predict(
        rf_inner,
        newdata = X_inner_val,
        type = "prob"
      )[, "Responder"]

      inner_metrics <- classification_metrics(
        observed = y_inner_val,
        predicted = inner_predicted,
        prob_responder = inner_probability
      )

      current_inner_record <- data.frame(
        outer_fold = outer_i,
        inner_fold = inner_i,
        mtry = candidate_mtry,
        training_samples = length(inner_train_local),
        validation_samples = length(inner_val_local),
        genes_after_filterByExpr =
          inner_selection$n_after_filter,
        Accuracy = inner_metrics$Accuracy,
        Sensitivity = inner_metrics$Sensitivity,
        Specificity = inner_metrics$Specificity,
        Balanced_Accuracy =
          inner_metrics$Balanced_Accuracy,
        AUC = inner_metrics$AUC,
        stringsAsFactors = FALSE
      )

      outer_inner_metrics <- rbind(
        outer_inner_metrics,
        current_inner_record
      )

      inner_cv_metrics <- rbind(
        inner_cv_metrics,
        current_inner_record
      )
    }
  }

  # ----------------------------------------------------------
  # 4B. Summarise inner-CV tuning
  # ----------------------------------------------------------

  tuning_rows <- lapply(
    mtry_grid,
    function(current_mtry) {

      candidate_rows <- outer_inner_metrics[
        outer_inner_metrics$mtry == current_mtry,
        ,
        drop = FALSE
      ]

      data.frame(
        outer_fold = outer_i,
        mtry = current_mtry,
        mean_Accuracy =
          mean(candidate_rows$Accuracy, na.rm = TRUE),
        mean_Sensitivity =
          mean(candidate_rows$Sensitivity, na.rm = TRUE),
        mean_Specificity =
          mean(candidate_rows$Specificity, na.rm = TRUE),
        mean_Balanced_Accuracy =
          mean(
            candidate_rows$Balanced_Accuracy,
            na.rm = TRUE
          ),
        mean_AUC =
          mean(candidate_rows$AUC, na.rm = TRUE),
        sd_Balanced_Accuracy =
          sd(
            candidate_rows$Balanced_Accuracy,
            na.rm = TRUE
          ),
        sd_AUC =
          sd(candidate_rows$AUC, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )

  tuning_table <- do.call(
    rbind,
    tuning_rows
  )

  tuning_table <- tuning_table[
    order(
      -tuning_table$mean_Balanced_Accuracy,
      -tuning_table$mean_AUC,
      tuning_table$mtry
    ),
    ,
    drop = FALSE
  ]

  best_mtry <- tuning_table$mtry[1]

  tuning_table$selected <- (
    tuning_table$mtry == best_mtry
  )

  inner_tuning_summary <- rbind(
    inner_tuning_summary,
    tuning_table
  )

  cat("\nInner-CV tuning summary:\n")
  print(tuning_table)

  cat(
    "\nSelected mtry:",
    best_mtry,
    "\n"
  )

  # ----------------------------------------------------------
  # 4C. Refit preprocessing and feature selection using only
  #     the complete OUTER-TRAINING partition
  # ----------------------------------------------------------

  outer_selection <- select_genes_from_training(
    counts_train = counts_outer_train,
    y_train = y_outer_train,
    top_n = top_n_genes
  )

  selected_outer_genes <- outer_selection$genes

  outer_selected_genes <- rbind(
    outer_selected_genes,
    data.frame(
      outer_fold = outer_i,
      rank = seq_along(selected_outer_genes),
      gene = selected_outer_genes,
      stringsAsFactors = FALSE
    )
  )

  X_outer_train <- t(
    samplewise_logcpm(
      counts_outer_train,
      selected_outer_genes
    )
  )

  X_outer_test <- t(
    samplewise_logcpm(
      counts_outer_test,
      selected_outer_genes
    )
  )

  if (!identical(
    colnames(X_outer_train),
    colnames(X_outer_test)
  )) {
    stop(
      "Outer train/test feature order differs.",
      call. = FALSE
    )
  }

  fitted_outer_mtry <- min(
    as.integer(best_mtry),
    ncol(X_outer_train)
  )

  # ----------------------------------------------------------
  # 4D. Fit OUTER model and predict held-out samples once
  # ----------------------------------------------------------

  set.seed(
    20262000 + outer_i
  )

  rf_outer <- randomForest(
    x = X_outer_train,
    y = y_outer_train,
    ntree = ntree,
    mtry = fitted_outer_mtry,
    importance = TRUE
  )

  outer_predicted <- predict(
    rf_outer,
    newdata = X_outer_test,
    type = "response"
  )

  outer_probability <- predict(
    rf_outer,
    newdata = X_outer_test,
    type = "prob"
  )[, "Responder"]

  current_predictions <- data.frame(
    sample_name = colnames(counts_outer_test),
    observed = as.character(y_outer_test),
    predicted = as.character(outer_predicted),
    prob_responder = as.numeric(outer_probability),
    outer_fold = outer_i,
    selected_mtry = fitted_outer_mtry,
    stringsAsFactors = FALSE
  )

  outer_predictions <- rbind(
    outer_predictions,
    current_predictions
  )

  current_outer_metrics <- classification_metrics(
    observed = y_outer_test,
    predicted = outer_predicted,
    prob_responder = outer_probability
  )

  current_outer_metrics <- data.frame(
    outer_fold = outer_i,
    training_samples = length(outer_train_idx),
    test_samples = length(outer_test_idx),
    genes_after_filterByExpr =
      outer_selection$n_after_filter,
    selected_genes =
      length(selected_outer_genes),
    selected_mtry =
      fitted_outer_mtry,
    current_outer_metrics,
    stringsAsFactors = FALSE
  )

  outer_fold_metrics <- rbind(
    outer_fold_metrics,
    current_outer_metrics
  )

  cat("\nOuter-fold performance:\n")
  print(current_outer_metrics)

  cat(
    "\n[PASS] OUTER FOLD",
    outer_i,
    "\n"
  )
}

# ------------------------------------------------------------
# 5. Validate pooled out-of-fold predictions
# ------------------------------------------------------------

if (nrow(outer_predictions) != ncol(counts)) {
  stop(
    paste0(
      "Expected 36 pooled outer predictions; found ",
      nrow(outer_predictions),
      "."
    ),
    call. = FALSE
  )
}

if (anyDuplicated(outer_predictions$sample_name)) {
  stop(
    "At least one sample received more than one outer prediction.",
    call. = FALSE
  )
}

if (!setequal(
  outer_predictions$sample_name,
  clinical$patient_id
)) {
  stop(
    "Pooled outer predictions do not cover the validated cohort exactly.",
    call. = FALSE
  )
}

outer_predictions <- outer_predictions[
  match(
    clinical$patient_id,
    outer_predictions$sample_name
  ),
  ,
  drop = FALSE
]

if (!identical(
  outer_predictions$sample_name,
  clinical$patient_id
)) {
  stop(
    "Failed to restore validated Week 1 sample order.",
    call. = FALSE
  )
}

if (!identical(
  outer_predictions$observed,
  as.character(y)
)) {
  stop(
    "Observed labels in pooled predictions differ from validated labels.",
    call. = FALSE
  )
}

if (
  any(outer_predictions$prob_responder < 0) ||
  any(outer_predictions$prob_responder > 1) ||
  anyNA(outer_predictions$prob_responder)
) {
  stop(
    "Invalid responder probabilities detected.",
    call. = FALSE
  )
}

pooled_metrics <- classification_metrics(
  observed = factor(
    outer_predictions$observed,
    levels = c("NonResponder", "Responder")
  ),
  predicted = factor(
    outer_predictions$predicted,
    levels = c("NonResponder", "Responder")
  ),
  prob_responder =
    outer_predictions$prob_responder
)

# ------------------------------------------------------------
# 6. Gene-selection stability across outer folds
# ------------------------------------------------------------

gene_frequency <- sort(
  table(outer_selected_genes$gene),
  decreasing = TRUE
)

gene_frequency_df <- data.frame(
  gene = names(gene_frequency),
  outer_folds_selected =
    as.integer(gene_frequency),
  selection_fraction =
    as.integer(gene_frequency) / outer_k,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 7. Save validation outputs
# ------------------------------------------------------------

write.csv(
  inner_cv_metrics,
  file.path(
    results_dir,
    "week3_inner_cv_metrics_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  inner_tuning_summary,
  file.path(
    results_dir,
    "week3_inner_tuning_summary_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  outer_predictions,
  file.path(
    results_dir,
    "week3_outer_predictions_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  outer_fold_metrics,
  file.path(
    results_dir,
    "week3_outer_fold_metrics_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  outer_selected_genes,
  file.path(
    results_dir,
    "week3_outer_selected_genes_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  gene_frequency_df,
  file.path(
    results_dir,
    "week3_gene_selection_frequency_strict.csv"
  ),
  row.names = FALSE
)

nested_cv_summary <- data.frame(
  metric = c(
    "samples",
    "input_genes",
    "outer_folds",
    "inner_folds",
    "genes_per_model",
    "accuracy",
    "sensitivity",
    "specificity",
    "balanced_accuracy",
    "auc"
  ),
  value = c(
    ncol(counts),
    nrow(counts),
    outer_k,
    inner_k,
    top_n_genes,
    pooled_metrics$Accuracy,
    pooled_metrics$Sensitivity,
    pooled_metrics$Specificity,
    pooled_metrics$Balanced_Accuracy,
    pooled_metrics$AUC
  ),
  stringsAsFactors = FALSE
)

write.csv(
  nested_cv_summary,
  file.path(
    results_dir,
    "week3_nested_cv_summary_strict.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. ROC curve from pooled OUT-OF-FOLD predictions only
# ------------------------------------------------------------

pooled_roc <- roc(
  response = factor(
    outer_predictions$observed,
    levels = c("NonResponder", "Responder")
  ),
  predictor =
    outer_predictions$prob_responder,
  levels = c("NonResponder", "Responder"),
  direction = "<",
  quiet = TRUE
)

png(
  filename = file.path(
    figures_dir,
    "week3_nested_cv_roc_strict.png"
  ),
  width = 2200,
  height = 1800,
  res = 250
)

plot(
  pooled_roc,
  legacy.axes = TRUE,
  main = paste0(
    "GSE160638 — strict nested CV ROC (AUC = ",
    sprintf("%.3f", as.numeric(auc(pooled_roc))),
    ")"
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

dev.off()

# ------------------------------------------------------------
# 9. Final Week 3 validation summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("WEEK 3 STRICT NESTED CROSS-VALIDATION: PASS\n")
cat("============================================================\n")

cat(
  "Validated samples:",
  nrow(outer_predictions),
  "\n"
)

cat(
  "Each sample predicted out-of-fold exactly once: TRUE\n"
)

cat(
  "Feature selection isolated inside training partitions: TRUE\n"
)

cat(
  "Hyperparameter tuning isolated inside outer training sets: TRUE\n"
)

cat(
  "Outer-test samples excluded from fitting and tuning: TRUE\n"
)

cat(
  "Accuracy:",
  sprintf("%.4f", pooled_metrics$Accuracy),
  "\n"
)

cat(
  "Sensitivity:",
  sprintf("%.4f", pooled_metrics$Sensitivity),
  "\n"
)

cat(
  "Specificity:",
  sprintf("%.4f", pooled_metrics$Specificity),
  "\n"
)

cat(
  "Balanced Accuracy:",
  sprintf("%.4f", pooled_metrics$Balanced_Accuracy),
  "\n"
)

cat(
  "AUC:",
  sprintf("%.4f", pooled_metrics$AUC),
  "\n"
)

cat(
  "Unique genes selected across outer folds:",
  nrow(gene_frequency_df),
  "\n"
)

cat("\n")
cat("Validation performance above is calculated exclusively\n")
cat("from pooled outer out-of-fold predictions.\n")
cat("No Week 2 exploratory expression matrix is used as\n")
cat("predictive model input.\n")
