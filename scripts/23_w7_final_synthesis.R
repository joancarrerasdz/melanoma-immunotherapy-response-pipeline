options(stringsAsFactors = FALSE)

cat("W7.6 — Final synthesis + portfolio figure\n\n")


# ============================================================
# Helpers
# ============================================================

read_csv <- function(path) {
  read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


# ============================================================
# Required W7 inputs
# ============================================================

required <- c(
  "results/w7_outer_fold_variability.csv",
  "results/w7_reference_baseline_audit.csv",
  "results/w7_metric_uncertainty.csv",
  "results/w7_transportability_audit.csv"
)

stopifnot(
  all(file.exists(required))
)

cat("Required W7 inputs: PASS\n\n")


# ============================================================
# Load frozen W7 outputs
# ============================================================

variability <- read_csv(
  "results/w7_outer_fold_variability.csv"
)

reference <- read_csv(
  "results/w7_reference_baseline_audit.csv"
)

uncertainty <- read_csv(
  "results/w7_metric_uncertainty.csv"
)

transport <- read_csv(
  "results/w7_transportability_audit.csv"
)


# ============================================================
# Final synthesis table
# ============================================================

metrics <- c(
  "Balanced_Accuracy",
  "AUC",
  "Brier"
)

external_scopes <- c(
  "W4_GSE91061",
  "W4_GSE78220",
  "W4_POOLED"
)

summary_rows <- list()

for (scope in external_scopes) {

  for (metric in metrics) {

    tr <- transport[
      transport$external_scope == scope &
        transport$metric == metric,
      ,
      drop = FALSE
    ]

    int_unc <- uncertainty[
      uncertainty$scope == "W3_OOF" &
        uncertainty$metric == metric,
      ,
      drop = FALSE
    ]

    ext_unc <- uncertainty[
      uncertainty$scope == scope &
        uncertainty$metric == metric,
      ,
      drop = FALSE
    ]

    ref <- reference[
      reference$scope == scope &
        reference$metric == metric,
      ,
      drop = FALSE
    ]

    stopifnot(
      nrow(tr) == 1L,
      nrow(int_unc) == 1L,
      nrow(ext_unc) == 1L,
      nrow(ref) == 1L
    )

    external_vs_reference <-
      if (metric == "Brier") {
        ext_unc$estimate < ref$reference
      } else {
        ext_unc$estimate > ref$reference
      }

    external_ci_informative <-
      if (metric == "Brier") {

        # Conservative descriptive criterion:
        # entire CI below reference Brier.
        ext_unc$ci95_upper < ref$reference

      } else {

        # Entire CI above no-skill/reference level.
        ext_unc$ci95_lower > ref$reference
      }

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(

      external_scope =
        scope,

      metric =
        metric,

      internal_estimate =
        tr$internal_estimate,

      internal_ci95_lower =
        int_unc$ci95_lower,

      internal_ci95_upper =
        int_unc$ci95_upper,

      external_estimate =
        tr$external_estimate,

      external_ci95_lower =
        ext_unc$ci95_lower,

      external_ci95_upper =
        ext_unc$ci95_upper,

      external_reference =
        ref$reference,

      external_beats_reference_point =
        external_vs_reference,

      external_ci_entirely_better_than_reference =
        external_ci_informative,

      degradation_gap =
        tr$degradation_gap,

      degradation_ci95_lower =
        tr$gap_ci95_lower,

      degradation_ci95_upper =
        tr$gap_ci95_upper,

      degradation_supported =
        tr$interpretation ==
          "EXTERNAL_DEGRADATION_SUPPORTED",

      interpretation =
        tr$interpretation,

      stringsAsFactors = FALSE
    )
  }
}

final_summary <- do.call(
  rbind,
  summary_rows
)

numeric_cols <- c(
  "internal_estimate",
  "internal_ci95_lower",
  "internal_ci95_upper",
  "external_estimate",
  "external_ci95_lower",
  "external_ci95_upper",
  "external_reference",
  "degradation_gap",
  "degradation_ci95_lower",
  "degradation_ci95_upper"
)

final_summary[numeric_cols] <-
  lapply(
    final_summary[numeric_cols],
    function(x) round(x, 12)
  )

write.csv(
  final_summary,
  "results/w7_final_synthesis.csv",
  row.names = FALSE,
  quote = TRUE
)

cat("Final synthesis table: PASS\n")


# ============================================================
# Key conclusions
# ============================================================

w3_auc <- uncertainty[
  uncertainty$scope == "W3_OOF" &
    uncertainty$metric == "AUC",
  ,
  drop = FALSE
]

w3_ba <- uncertainty[
  uncertainty$scope == "W3_OOF" &
    uncertainty$metric == "Balanced_Accuracy",
  ,
  drop = FALSE
]

pooled_auc <- final_summary[
  final_summary$external_scope == "W4_POOLED" &
    final_summary$metric == "AUC",
  ,
  drop = FALSE
]

pooled_ba <- final_summary[
  final_summary$external_scope == "W4_POOLED" &
    final_summary$metric == "Balanced_Accuracy",
  ,
  drop = FALSE
]

pooled_brier <- final_summary[
  final_summary$external_scope == "W4_POOLED" &
    final_summary$metric == "Brier",
  ,
  drop = FALSE
]

stopifnot(
  nrow(w3_auc) == 1L,
  nrow(w3_ba) == 1L,
  nrow(pooled_auc) == 1L,
  nrow(pooled_ba) == 1L,
  nrow(pooled_brier) == 1L
)

conclusions <- data.frame(
  finding = c(
    "Internal discrimination",
    "Internal balanced classification",
    "Pooled external discrimination",
    "Pooled external balanced classification",
    "Pooled external probability quality"
  ),
  value = c(
    sprintf(
      "W3 OOF AUC %.3f [%.3f, %.3f]",
      w3_auc$estimate,
      w3_auc$ci95_lower,
      w3_auc$ci95_upper
    ),
    sprintf(
      "W3 OOF balanced accuracy %.3f [%.3f, %.3f]",
      w3_ba$estimate,
      w3_ba$ci95_lower,
      w3_ba$ci95_upper
    ),
    sprintf(
      "AUC degradation %.3f [%.3f, %.3f]",
      pooled_auc$degradation_gap,
      pooled_auc$degradation_ci95_lower,
      pooled_auc$degradation_ci95_upper
    ),
    sprintf(
      "Balanced-accuracy degradation %.3f [%.3f, %.3f]",
      pooled_ba$degradation_gap,
      pooled_ba$degradation_ci95_lower,
      pooled_ba$degradation_ci95_upper
    ),
    sprintf(
      "Brier degradation %.3f [%.3f, %.3f]",
      pooled_brier$degradation_gap,
      pooled_brier$degradation_ci95_lower,
      pooled_brier$degradation_ci95_upper
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  conclusions,
  "results/w7_key_conclusions.csv",
  row.names = FALSE,
  quote = TRUE
)


# ============================================================
# Portfolio figure
# ============================================================

if (!requireNamespace(
  "ggplot2",
  quietly = TRUE
)) {
  stop("ggplot2 is required for the W7 portfolio figure.")
}

plot_df <- final_summary[
  final_summary$metric %in%
    c(
      "Balanced_Accuracy",
      "AUC",
      "Brier"
    ),
  ,
  drop = FALSE
]

plot_df$metric_label <-
  factor(
    plot_df$metric,
    levels = c(
      "Balanced_Accuracy",
      "AUC",
      "Brier"
    ),
    labels = c(
      "Balanced accuracy",
      "ROC AUC",
      "Brier score"
    )
  )

plot_df$external_label <-
  factor(
    plot_df$external_scope,
    levels = c(
      "W4_GSE91061",
      "W4_GSE78220",
      "W4_POOLED"
    ),
    labels = c(
      "GSE91061",
      "GSE78220",
      "Pooled external"
    )
  )

p <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = degradation_gap,
    y = external_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = degradation_ci95_lower,
      xmax = degradation_ci95_upper
    ),
    width = 0.18,
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2.7
  ) +
  ggplot2::facet_wrap(
    ~metric_label,
    scales = "free_x",
    nrow = 1
  ) +
  ggplot2::labs(
    title =
      "Internal-to-external performance degradation",
    subtitle =
      paste(
        "Positive values indicate poorer external performance;",
        "points are frozen W3→W4 gaps with bootstrap 95% intervals"
      ),
    x =
      "External degradation gap",
    y =
      NULL,
    caption =
      paste(
        "W7 post-hoc audit; no retraining, threshold retuning",
        "or modification of frozen W1–W4 outputs."
      )
  ) +
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    strip.text =
      ggplot2::element_text(
        face = "bold"
      )
  )

ggplot2::ggsave(
  filename =
    "figures/w7_transportability_uncertainty.png",
  plot = p,
  width = 10.5,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("Portfolio figure: PASS\n")


# ============================================================
# Final documentation
# ============================================================

doc_path <-
  "docs/w7_final_synthesis.md"

doc <- c(
  "# W7 — Final evaluation and uncertainty synthesis",
  "",
  "## Purpose",
  "",
  paste(
    "Week 7 is a post-hoc evaluation audit of the frozen Week 3",
    "strict internal validation and Week 4 locked external validation."
  ),
  "",
  paste(
    "No model was retrained, no feature set was reselected,",
    "no threshold was retuned and no W1-W4 analytical output was modified."
  ),
  "",
  "## Internal performance",
  "",
  sprintf(
    "- W3 OOF ROC AUC: %.3f (95%% bootstrap CI %.3f–%.3f).",
    w3_auc$estimate,
    w3_auc$ci95_lower,
    w3_auc$ci95_upper
  ),
  sprintf(
    "- W3 OOF balanced accuracy: %.3f (95%% bootstrap CI %.3f–%.3f).",
    w3_ba$estimate,
    w3_ba$ci95_lower,
    w3_ba$ci95_upper
  ),
  "",
  "These results support strong internal predictive discrimination.",
  "",
  "## External transportability",
  "",
  sprintf(
    paste(
      "- Pooled external AUC degradation:",
      "%.3f (95%% bootstrap CI %.3f–%.3f)."
    ),
    pooled_auc$degradation_gap,
    pooled_auc$degradation_ci95_lower,
    pooled_auc$degradation_ci95_upper
  ),
  sprintf(
    paste(
      "- Pooled external balanced-accuracy degradation:",
      "%.3f (95%% bootstrap CI %.3f–%.3f)."
    ),
    pooled_ba$degradation_gap,
    pooled_ba$degradation_ci95_lower,
    pooled_ba$degradation_ci95_upper
  ),
  sprintf(
    paste(
      "- Pooled external Brier degradation:",
      "%.3f (95%% bootstrap CI %.3f–%.3f)."
    ),
    pooled_brier$degradation_gap,
    pooled_brier$degradation_ci95_lower,
    pooled_brier$degradation_ci95_upper
  ),
  "",
  paste(
    "The pooled comparison therefore supports degradation of",
    "discrimination, balanced classification and probability quality",
    "when moving from internal to external evaluation."
  ),
  "",
  "## Cohort-specific interpretation",
  "",
  paste(
    "GSE91061 shows supported degradation for balanced accuracy,",
    "ROC AUC and Brier score."
  ),
  "",
  paste(
    "GSE78220 shows a poorer external point estimate for all three",
    "metrics, but the balanced-accuracy and ROC-AUC gap intervals",
    "include zero; Brier degradation remains supported."
  ),
  "",
  "## Interpretation boundary",
  "",
  paste(
    "These results demonstrate limited transportability of the frozen",
    "predictive performance across the external datasets evaluated here."
  ),
  "",
  paste(
    "They do not establish that dataset shift alone caused the",
    "performance differences, and they should not be interpreted as",
    "a causal comparison between cohorts."
  ),
  "",
  paste(
    "Differences in cohort composition, responder prevalence,",
    "assay context and study-specific characteristics remain plausible",
    "contributors to the observed external degradation."
  ),
  "",
  "## Portfolio conclusion",
  "",
  paste(
    "The final model shows strong internal discrimination but",
    "substantially weaker and more uncertain external performance."
  ),
  "",
  paste(
    "The principal methodological conclusion is therefore not simply",
    "that the model performs well internally, but that external",
    "transportability is the dominant limitation revealed by the",
    "reproducible evaluation."
  ),
  "",
  "## Generated artifacts",
  "",
  "- `results/w7_final_synthesis.csv`",
  "- `results/w7_key_conclusions.csv`",
  "- `figures/w7_transportability_uncertainty.png`",
  "- `docs/w7_final_synthesis.md`",
  ""
)

writeLines(
  doc,
  doc_path,
  useBytes = TRUE
)

cat("Final synthesis documentation: PASS\n")


# ============================================================
# Final console summary
# ============================================================

cat("\nW7 KEY CONCLUSIONS\n\n")

print(
  conclusions,
  row.names = FALSE
)

cat("\nW7.6 OUTPUTS\n")
cat("- results/w7_final_synthesis.csv\n")
cat("- results/w7_key_conclusions.csv\n")
cat("- figures/w7_transportability_uncertainty.png\n")
cat("- docs/w7_final_synthesis.md\n")

cat("\n")
cat("W7.6 FINAL SYNTHESIS: PASS\n")
cat("W7.6: PASS\n")
