options(stringsAsFactors = FALSE)

cat("W7.4 — Reference baselines + metric uncertainty\n\n")

set.seed(20261016)


# ============================================================
# Helpers
# ============================================================

read_frozen <- function(path) {
  read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


auc_rank <- function(y, p) {

  y <- as.integer(y)

  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)

  stopifnot(n1 > 0L, n0 > 0L)

  ranks <- rank(
    p,
    ties.method = "average"
  )

  (
    sum(ranks[y == 1L]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}


metric_vector <- function(
  dat,
  threshold = NULL,
  use_frozen_prediction = FALSE
) {

  y <- as.integer(
    dat$observed == "Responder"
  )

  p <- as.numeric(
    dat$prob_responder
  )

  if (use_frozen_prediction) {

    pred <- as.integer(
      dat$predicted == "Responder"
    )

  } else {

    stopifnot(!is.null(threshold))

    pred <- as.integer(
      p >= threshold
    )
  }

  tp <- sum(y == 1L & pred == 1L)
  tn <- sum(y == 0L & pred == 0L)
  fp <- sum(y == 0L & pred == 1L)
  fn <- sum(y == 1L & pred == 0L)

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)

  c(
    Accuracy = mean(pred == y),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Balanced_Accuracy =
      (sensitivity + specificity) / 2,
    AUC = auc_rank(y, p),
    Brier = mean((p - y)^2)
  )
}


stratified_indices <- function(dat, strata) {

  key <- interaction(
    dat[strata],
    drop = TRUE,
    lex.order = TRUE
  )

  groups <- split(
    seq_len(nrow(dat)),
    key
  )

  unlist(
    lapply(
      groups,
      function(idx) {
        sample(
          idx,
          size = length(idx),
          replace = TRUE
        )
      }
    ),
    use.names = FALSE
  )
}


bootstrap_ci <- function(
  dat,
  threshold = NULL,
  use_frozen_prediction = FALSE,
  strata,
  B = 5000L
) {

  point <- metric_vector(
    dat,
    threshold = threshold,
    use_frozen_prediction =
      use_frozen_prediction
  )

  boot <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(point)
  )

  colnames(boot) <- names(point)

  for (b in seq_len(B)) {

    idx <- stratified_indices(
      dat,
      strata
    )

    boot[b, ] <- metric_vector(
      dat[idx, , drop = FALSE],
      threshold = threshold,
      use_frozen_prediction =
        use_frozen_prediction
    )
  }

  lower <- apply(
    boot,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE,
    names = FALSE
  )

  upper <- apply(
    boot,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE,
    names = FALSE
  )

  list(
    estimate = point,
    lower = lower,
    upper = upper
  )
}


reference_values <- function(dat) {

  y <- as.integer(
    dat$observed == "Responder"
  )

  prevalence <- mean(y)

  c(
    Accuracy =
      max(prevalence, 1 - prevalence),

    Balanced_Accuracy = 0.5,

    AUC = 0.5,

    Brier =
      prevalence * (1 - prevalence)
  )
}


# ============================================================
# Frozen inputs
# ============================================================

w3 <- read_frozen(
  "results/week3_outer_predictions_strict.csv"
)

w4 <- read_frozen(
  "results/week4_external_validation_predictions_strict.csv"
)

stopifnot(
  nrow(w3) == 36L,
  nrow(w4) == 76L,
  length(unique(w3$sample_name)) == 36L
)


# ============================================================
# Scopes
# ============================================================

scopes <- list(

  W3_OOF = list(
    data = w3,
    threshold = NULL,
    frozen_prediction = TRUE,
    strata = c("observed"),
    bootstrap_scheme =
      "outcome-stratified"
  ),

  W4_GSE91061 = list(
    data = w4[
      w4$dataset == "GSE91061",
      ,
      drop = FALSE
    ],
    threshold = 0.48,
    frozen_prediction = FALSE,
    strata = c("observed"),
    bootstrap_scheme =
      "outcome-stratified"
  ),

  W4_GSE78220 = list(
    data = w4[
      w4$dataset == "GSE78220",
      ,
      drop = FALSE
    ],
    threshold = 0.48,
    frozen_prediction = FALSE,
    strata = c("observed"),
    bootstrap_scheme =
      "outcome-stratified"
  ),

  W4_POOLED = list(
    data = w4,
    threshold = 0.48,
    frozen_prediction = FALSE,
    strata = c(
      "dataset",
      "observed"
    ),
    bootstrap_scheme =
      "dataset-and-outcome-stratified"
  )
)


# ============================================================
# Reference baseline audit
# ============================================================

reference_rows <- list()

for (scope_name in names(scopes)) {

  cfg <- scopes[[scope_name]]
  dat <- cfg$data

  actual <- metric_vector(
    dat,
    threshold = cfg$threshold,
    use_frozen_prediction =
      cfg$frozen_prediction
  )

  ref <- reference_values(dat)

  prevalence <- mean(
    dat$observed == "Responder"
  )

  for (
    metric in
    c(
      "Accuracy",
      "Balanced_Accuracy",
      "AUC",
      "Brier"
    )
  ) {

    actual_value <- actual[[metric]]
    reference_value <- ref[[metric]]

    if (metric == "Brier") {

      direction <- "lower_is_better"

      margin <-
        reference_value -
        actual_value

    } else {

      direction <- "higher_is_better"

      margin <-
        actual_value -
        reference_value
    }

    reference_rows[[length(reference_rows) + 1L]] <- data.frame(

      scope = scope_name,

      n = nrow(dat),

      responders =
        sum(
          dat$observed ==
            "Responder"
        ),

      nonresponders =
        sum(
          dat$observed ==
            "NonResponder"
        ),

      prevalence_responder =
        prevalence,

      metric = metric,

      actual =
        actual_value,

      reference =
        reference_value,

      direction =
        direction,

      margin_vs_reference =
        margin,

      stringsAsFactors = FALSE
    )
  }
}

reference_audit <- do.call(
  rbind,
  reference_rows
)

numeric_cols <- c(
  "prevalence_responder",
  "actual",
  "reference",
  "margin_vs_reference"
)

reference_audit[numeric_cols] <-
  lapply(
    reference_audit[numeric_cols],
    function(x) round(x, 12)
  )

write.csv(
  reference_audit,
  "results/w7_reference_baseline_audit.csv",
  row.names = FALSE,
  quote = TRUE
)


# ============================================================
# Prediction-level bootstrap uncertainty
# ============================================================

uncertainty_rows <- list()

B <- 5000L

for (scope_name in names(scopes)) {

  cfg <- scopes[[scope_name]]

  cat(
    "Bootstrapping",
    scope_name,
    "... "
  )

  res <- bootstrap_ci(
    cfg$data,
    threshold =
      cfg$threshold,
    use_frozen_prediction =
      cfg$frozen_prediction,
    strata =
      cfg$strata,
    B = B
  )

  cat("done\n")

  for (
    metric in names(res$estimate)
  ) {

    uncertainty_rows[[length(uncertainty_rows) + 1L]] <- data.frame(

      scope =
        scope_name,

      metric =
        metric,

      estimate =
        res$estimate[[metric]],

      ci95_lower =
        res$lower[[metric]],

      ci95_upper =
        res$upper[[metric]],

      bootstrap_replicates =
        B,

      bootstrap_scheme =
        cfg$bootstrap_scheme,

      stringsAsFactors = FALSE
    )
  }
}

uncertainty <- do.call(
  rbind,
  uncertainty_rows
)

numeric_cols <- c(
  "estimate",
  "ci95_lower",
  "ci95_upper"
)

uncertainty[numeric_cols] <-
  lapply(
    uncertainty[numeric_cols],
    function(x) round(x, 12)
  )

write.csv(
  uncertainty,
  "results/w7_metric_uncertainty.csv",
  row.names = FALSE,
  quote = TRUE
)


# ============================================================
# Documentation extension
# ============================================================

doc_path <-
  "docs/w7_evaluation_uncertainty_audit.md"

doc <- c(
  "",
  "## W7.4 reference baselines",
  "",
  paste(
    "Frozen predictive performance is compared with",
    "simple non-informative or prevalence-aware references."
  ),
  "",
  "- Accuracy: majority-class classifier.",
  "- Balanced accuracy: 0.50.",
  "- ROC AUC: 0.50.",
  "- Brier score: constant predictor equal to observed responder prevalence.",
  "",
  paste(
    "For Brier score, lower values are better;",
    "for the remaining metrics, higher values are better."
  ),
  "",
  "Results are recorded in `results/w7_reference_baseline_audit.csv`.",
  "",
  "## W7.4 prediction-level uncertainty",
  "",
  paste(
    "Non-parametric bootstrap confidence intervals are calculated",
    "from already-frozen prediction records."
  ),
  "",
  paste(
    "Week 3 bootstrap samples are stratified by outcome.",
    "Week 4 cohort-specific samples are stratified by outcome,",
    "and the pooled Week 4 analysis is stratified jointly by dataset and outcome."
  ),
  "",
  paste(
    "The bootstrap uses 5,000 deterministic replicates",
    "with seed 20261016."
  ),
  "",
  paste(
    "These intervals characterize evaluation uncertainty conditional",
    "on the frozen predictions. They do not include model-training,",
    "feature-selection or hyperparameter-selection uncertainty."
  ),
  "",
  "Results are recorded in `results/w7_metric_uncertainty.csv`.",
  ""
)

existing_doc <- if (file.exists(doc_path)) {
  readLines(doc_path, warn = FALSE)
} else {
  character()
}

if (!any(grepl(
  "^## W7\\.4 reference baselines$",
  existing_doc
))) {

  cat(
    paste(doc, collapse = "\n"),
    "\n",
    file = doc_path,
    append = TRUE,
    sep = ""
  )
} else {

  cat(
    "Documentation section W7.4 already present; ",
    "no duplicate append performed.\n"
  )
}


# ============================================================
# Console summary
# ============================================================

cat("\nREFERENCE BASELINE AUDIT\n\n")

print(
  reference_audit[
    ,
    c(
      "scope",
      "metric",
      "actual",
      "reference",
      "margin_vs_reference"
    )
  ],
  row.names = FALSE
)

cat("\nMETRIC UNCERTAINTY\n\n")

print(
  uncertainty[
    uncertainty$metric %in%
      c(
        "Accuracy",
        "Balanced_Accuracy",
        "AUC",
        "Brier"
      ),
    ,
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\n")
cat("W7.4 REFERENCE BASELINE AUDIT: PASS\n")
cat("W7.4 METRIC UNCERTAINTY: PASS\n")
cat("W7.4: PASS\n")
