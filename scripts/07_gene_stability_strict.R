#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

input_file <- "results/week3_gene_selection_frequency_strict.csv"
results_dir <- "results"
figures_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Load strict nested-CV selection frequencies
# ------------------------------------------------------------

x <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "gene",
  "outer_folds_selected",
  "selection_fraction"
)

missing_cols <- setdiff(required_cols, colnames(x))

if (length(missing_cols) > 0L) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (nrow(x) != 312L) {
  stop(
    paste0(
      "Unexpected number of selected genes: ",
      nrow(x),
      ". Expected 312."
    ),
    call. = FALSE
  )
}

if (
  any(x$outer_folds_selected < 1L) ||
  any(x$outer_folds_selected > 5L)
) {
  stop("Invalid outer-fold selection counts.", call. = FALSE)
}

# ------------------------------------------------------------
# 2. Define stability tiers
# ------------------------------------------------------------

x$stability_tier <- ifelse(
  x$outer_folds_selected == 5L,
  "core_consensus_5of5",
  ifelse(
    x$outer_folds_selected >= 4L,
    "high_stability_4of5",
    ifelse(
      x$outer_folds_selected >= 3L,
      "majority_stable_3of5",
      ifelse(
        x$outer_folds_selected >= 2L,
        "repeated_2of5",
        "single_fold_1of5"
      )
    )
  )
)

x <- x[
  order(
    -x$outer_folds_selected,
    -x$selection_fraction,
    x$gene
  ),
]

# ------------------------------------------------------------
# 3. Consensus candidate sets
# ------------------------------------------------------------

repeated <- x[x$outer_folds_selected >= 2L, , drop = FALSE]
majority <- x[x$outer_folds_selected >= 3L, , drop = FALSE]
high <- x[x$outer_folds_selected >= 4L, , drop = FALSE]
core <- x[x$outer_folds_selected == 5L, , drop = FALSE]

stopifnot(
  nrow(x) == 312L,
  nrow(repeated) == 100L,
  nrow(majority) == 51L,
  nrow(high) == 25L,
  nrow(core) == 12L
)

# ------------------------------------------------------------
# 4. Save tables
# ------------------------------------------------------------

write.csv(
  x,
  file.path(
    results_dir,
    "week3_gene_stability_tiers_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  repeated,
  file.path(
    results_dir,
    "week3_consensus_repeated_2of5_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  majority,
  file.path(
    results_dir,
    "week3_consensus_majority_3of5_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  high,
  file.path(
    results_dir,
    "week3_consensus_high_4of5_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  core,
  file.path(
    results_dir,
    "week3_consensus_core_5of5_strict.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Stability distribution figure
# ------------------------------------------------------------

freq <- table(
  factor(
    x$outer_folds_selected,
    levels = 1:5
  )
)

png(
  filename = file.path(
    figures_dir,
    "week3_gene_selection_stability_strict.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

barplot(
  freq,
  names.arg = paste0(1:5, "/5"),
  xlab = "Outer folds in which gene was selected",
  ylab = "Number of genes",
  main = "GSE160638 — nested-CV gene-selection stability"
)

text(
  seq_along(freq) * 1.2 - 0.5,
  as.numeric(freq),
  labels = as.numeric(freq),
  pos = 3
)

dev.off()

# ------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("WEEK 3 GENE-SELECTION STABILITY: PASS\n")
cat("============================================================\n")

cat("Unique selected genes:", nrow(x), "\n")
cat("Selected in >= 2/5 folds:", nrow(repeated), "\n")
cat("Selected in >= 3/5 folds:", nrow(majority), "\n")
cat("Selected in >= 4/5 folds:", nrow(high), "\n")
cat("Selected in 5/5 folds:", nrow(core), "\n")

cat("\nCore-consensus candidates (5/5):\n")
print(core[, c(
  "gene",
  "outer_folds_selected",
  "selection_fraction"
)])

cat("\nIMPORTANT:\n")
cat(
  "These stability tiers summarize feature-selection recurrence across\n"
)
cat(
  "the outer folds of the strict nested-CV procedure.\n"
)
cat(
  "They do not constitute an independently validated final signature.\n"
)
cat(
  "Performance of any consensus signature fitted after this step must\n"
)
cat(
  "not be reported as independent validation on these same 36 samples.\n"
)

cat("\nOutputs written successfully.\n")
