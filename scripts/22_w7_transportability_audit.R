options(stringsAsFactors = FALSE)

cat("W7.5 — Internal vs external transportability audit\n\n")

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

  stopifnot(
    n1 > 0L,
    n0 > 0L
  )

  ranks <- rank(
    p,
    ties.method = "average"
  )

  (
    sum(ranks[y == 1L]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}


metric_vector <- function(dat) {

  y <- as.integer(
    dat$observed == "Responder"
  )

  pred <- as.integer(
    dat$predicted == "Responder"
  )

  p <- as.numeric(
    dat$prob_responder
  )

  tp <- sum(y == 1L & pred == 1L)
  tn <- sum(y == 0L & pred == 0L)
  fp <- sum(y == 0L & pred == 1L)
  fn <- sum(y == 1L & pred == 0L)

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)

  c(
    Accuracy =
      mean(pred == y),

    Balanced_Accuracy =
      (sensitivity + specificity) / 2,

    AUC =
      auc_rank(y, p),

    Brier =
      mean((p - y)^2)
  )
}


stratified_indices <- function(
  dat,
  strata
) {

  strata_data <- dat[
    ,
    strata,
    drop = FALSE
  ]

  key <- do.call(
    interaction,
    c(
      strata_data,
      list(
        drop = TRUE,
        lex.order = TRUE
      )
    )
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


transport_gap <- function(
  internal,
  external
) {

  stopifnot(
    identical(
      names(internal),
      names(external)
    )
  )

  out <- internal - external

  # Brier is lower-is-better.
  # Positive gap must still mean poorer external performance.
  out["Brier"] <-
    external["Brier"] -
    internal["Brier"]

  out
}


bootstrap_transport <- function(
  internal,
  external,
  external_strata,
  B = 5000L
) {

  internal_point <-
    metric_vector(internal)

  external_point <-
    metric_vector(external)

  point <- transport_gap(
    internal_point,
    external_point
  )

  boot <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(point)
  )

  colnames(boot) <-
    names(point)

  for (b in seq_len(B)) {

    idx_internal <-
      stratified_indices(
        internal,
        "observed"
      )

    idx_external <-
      stratified_indices(
        external,
        external_strata
      )

    int_metrics <-
      metric_vector(
        internal[
          idx_internal,
          ,
          drop = FALSE
        ]
      )

    ext_metrics <-
      metric_vector(
        external[
          idx_external,
          ,
          drop = FALSE
        ]
      )

    boot[b, ] <-
      transport_gap(
        int_metrics,
        ext_metrics
      )
  }

  lower <- apply(
    boot,
    2,
    quantile,
    probs = 0.025,
    names = FALSE,
    na.rm = TRUE
  )

  upper <- apply(
    boot,
    2,
    quantile,
    probs = 0.975,
    names = FALSE,
    na.rm = TRUE
  )

  list(
    internal = internal_point,
    external = external_point,
    gap = point,
    lower = lower,
    upper = upper
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

cat("Frozen inputs: PASS\n")
cat("W3 OOF samples :", nrow(w3), "\n")
cat("W4 samples     :", nrow(w4), "\n\n")


# ============================================================
# External scopes
# ============================================================

external_scopes <- list(

  W4_GSE91061 = list(
    data = w4[
      w4$dataset == "GSE91061",
      ,
      drop = FALSE
    ],
    strata = "observed",
    scheme =
      "independent outcome-stratified bootstrap"
  ),

  W4_GSE78220 = list(
    data = w4[
      w4$dataset == "GSE78220",
      ,
      drop = FALSE
    ],
    strata = "observed",
    scheme =
      "independent outcome-stratified bootstrap"
  ),

  W4_POOLED = list(
    data = w4,
    strata = c(
      "dataset",
      "observed"
    ),
    scheme =
      "independent dataset-and-outcome-stratified bootstrap"
  )
)


# ============================================================
# Transportability audit
# ============================================================

B <- 5000L

rows <- list()

for (scope_name in names(external_scopes)) {

  cfg <- external_scopes[[scope_name]]

  cat(
    "Comparing W3_OOF vs",
    scope_name,
    "... "
  )

  res <- bootstrap_transport(
    internal = w3,
    external = cfg$data,
    external_strata = cfg$strata,
    B = B
  )

  cat("done\n")

  internal_prev <-
    mean(
      w3$observed ==
        "Responder"
    )

  external_prev <-
    mean(
      cfg$data$observed ==
        "Responder"
    )

  for (metric in names(res$gap)) {

    gap <- unname(
      res$gap[[metric]]
    )

    lower <- unname(
      res$lower[[metric]]
    )

    upper <- unname(
      res$upper[[metric]]
    )

    interpretation <-
      if (lower > 0) {

        "EXTERNAL_DEGRADATION_SUPPORTED"

      } else if (upper < 0) {

        "EXTERNAL_IMPROVEMENT_SUPPORTED"

      } else if (gap > 0) {

        "DEGRADATION_DIRECTION_UNCERTAIN"

      } else {

        "NO_DEGRADATION_DIRECTION_UNCERTAIN"
      }

    rows[[length(rows) + 1L]] <- data.frame(

      internal_scope =
        "W3_OOF",

      external_scope =
        scope_name,

      metric =
        metric,

      internal_n =
        nrow(w3),

      external_n =
        nrow(cfg$data),

      internal_responder_prevalence =
        internal_prev,

      external_responder_prevalence =
        external_prev,

      internal_estimate =
        unname(
          res$internal[[metric]]
        ),

      external_estimate =
        unname(
          res$external[[metric]]
        ),

      degradation_gap =
        gap,

      gap_ci95_lower =
        lower,

      gap_ci95_upper =
        upper,

      gap_ci_excludes_zero =
        lower > 0 ||
        upper < 0,

      interpretation =
        interpretation,

      bootstrap_replicates =
        B,

      bootstrap_scheme =
        cfg$scheme,

      stringsAsFactors = FALSE
    )
  }
}

transport <- do.call(
  rbind,
  rows
)

numeric_cols <- c(
  "internal_responder_prevalence",
  "external_responder_prevalence",
  "internal_estimate",
  "external_estimate",
  "degradation_gap",
  "gap_ci95_lower",
  "gap_ci95_upper"
)

transport[numeric_cols] <-
  lapply(
    transport[numeric_cols],
    function(x) round(x, 12)
  )

write.csv(
  transport,
  "results/w7_transportability_audit.csv",
  row.names = FALSE,
  quote = TRUE
)


# ============================================================
# Internal consistency against W7.4
# ============================================================

uncertainty <- read_frozen(
  "results/w7_metric_uncertainty.csv"
)

for (
  scope_name in
  c(
    "W4_GSE91061",
    "W4_GSE78220",
    "W4_POOLED"
  )
) {

  for (
    metric in
    c(
      "Accuracy",
      "Balanced_Accuracy",
      "AUC",
      "Brier"
    )
  ) {

    expected_internal <-
      uncertainty$estimate[
        uncertainty$scope == "W3_OOF" &
          uncertainty$metric == metric
      ]

    expected_external <-
      uncertainty$estimate[
        uncertainty$scope == scope_name &
          uncertainty$metric == metric
      ]

    observed_row <-
      transport[
        transport$external_scope == scope_name &
          transport$metric == metric,
        ,
        drop = FALSE
      ]

    stopifnot(
      nrow(observed_row) == 1L,
      length(expected_internal) == 1L,
      length(expected_external) == 1L,
      isTRUE(
        all.equal(
          observed_row$internal_estimate,
          expected_internal,
          tolerance = 1e-10
        )
      ),
      isTRUE(
        all.equal(
          observed_row$external_estimate,
          expected_external,
          tolerance = 1e-10
        )
      )
    )
  }
}

cat(
  "\nW7.4 / W7.5 point-estimate consistency: PASS\n"
)


# ============================================================
# Documentation — idempotent append
# ============================================================

doc_path <-
  "docs/w7_evaluation_uncertainty_audit.md"

existing_doc <-
  readLines(
    doc_path,
    warn = FALSE
  )

marker <-
  "## W7.5 internal-versus-external transportability"

if (!any(existing_doc == marker)) {

  doc <- c(
    "",
    marker,
    "",
    paste(
      "Transportability is assessed by comparing the frozen strict Week 3",
      "out-of-fold predictions with the frozen locked Week 4 external",
      "validation predictions."
    ),
    "",
    paste(
      "For Accuracy, Balanced Accuracy and ROC AUC, the degradation gap is",
      "defined as internal minus external performance."
    ),
    "",
    paste(
      "For Brier score, where lower values are preferable, the degradation",
      "gap is defined as external minus internal."
    ),
    "",
    paste(
      "Therefore, a positive degradation gap consistently indicates poorer",
      "external performance."
    ),
    "",
    paste(
      "Gap uncertainty is estimated using 5,000 independent stratified",
      "bootstrap resamples of the frozen internal and external prediction",
      "records."
    ),
    "",
    paste(
      "These comparisons quantify empirical performance transportability.",
      "They are not interpreted as causal evidence that dataset shift alone",
      "caused the observed differences."
    ),
    "",
    paste(
      "Differences in cohort composition, responder prevalence, assay context",
      "and other study-specific factors remain part of the interpretation."
    ),
    "",
    "Results are recorded in `results/w7_transportability_audit.csv`.",
    ""
  )

  cat(
    paste(
      doc,
      collapse = "\n"
    ),
    "\n",
    file = doc_path,
    append = TRUE,
    sep = ""
  )

  cat(
    "W7.5 documentation appended.\n"
  )

} else {

  cat(
    "W7.5 documentation already present; ",
    "no duplicate append performed.\n"
  )
}


# ============================================================
# Console summary
# ============================================================

cat("\nTRANSPORTABILITY SUMMARY\n\n")

print(
  transport[
    ,
    c(
      "external_scope",
      "metric",
      "internal_estimate",
      "external_estimate",
      "degradation_gap",
      "gap_ci95_lower",
      "gap_ci95_upper",
      "interpretation"
    )
  ],
  row.names = FALSE
)

cat("\n")
cat("W7.5 TRANSPORTABILITY AUDIT: PASS\n")
cat("W7.5: PASS\n")
