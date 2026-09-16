#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

strict_core_file <-
  "results/week3_consensus_core_5of5_strict.csv"

strict_stability_file <-
  "results/week3_gene_stability_tiers_strict.csv"

counts_file <-
  "data_processed/GSE160638_raw_counts_PD1_aligned.rds"

historical_candidates <- c(
  "results/stable_genes_5of5_final_100.csv",
  "results/stable_genes_5of5_with_importance.csv"
)

historical_file <- historical_candidates[
  file.exists(historical_candidates)
][1]

required_files <- c(
  strict_core_file,
  strict_stability_file,
  counts_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    paste(
      "Missing required files:",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}

if (is.na(historical_file)) {
  stop(
    "Historical 5/5 stability output not found.",
    call. = FALSE
  )
}

historical <- read.csv(
  historical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

strict_core <- read.csv(
  strict_core_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

strict_stability <- read.csv(
  strict_stability_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts <- readRDS(counts_file)

gene_candidates <- c(
  "gene",
  "Gene",
  "symbol",
  "SYMBOL",
  "gene_symbol"
)

historical_gene_col <- intersect(
  gene_candidates,
  colnames(historical)
)

if (length(historical_gene_col) == 0L) {
  stop(
    "Could not identify historical gene column.",
    call. = FALSE
  )
}

historical_gene_col <- historical_gene_col[1]

historical_genes <- sort(
  unique(
    trimws(
      as.character(
        historical[[historical_gene_col]]
      )
    )
  )
)

strict_core_genes <- sort(
  unique(trimws(as.character(strict_core$gene)))
)

strict_stability$gene <-
  trimws(as.character(strict_stability$gene))

# ------------------------------------------------------------
# 1. Historical 5/5 versus strict 5/5
# ------------------------------------------------------------

all_core_genes <- sort(
  union(historical_genes, strict_core_genes)
)

comparison <- data.frame(
  gene = all_core_genes,
  historical_5of5 =
    all_core_genes %in% historical_genes,
  strict_5of5 =
    all_core_genes %in% strict_core_genes,
  stringsAsFactors = FALSE
)

comparison$category <- ifelse(
  comparison$historical_5of5 &
    comparison$strict_5of5,
  "shared",
  ifelse(
    comparison$strict_5of5,
    "strict_only",
    "historical_only"
  )
)

shared <- intersect(
  historical_genes,
  strict_core_genes
)

historical_only <- setdiff(
  historical_genes,
  strict_core_genes
)

strict_only <- setdiff(
  strict_core_genes,
  historical_genes
)

comparison_summary <- data.frame(
  metric = c(
    "historical_5of5_genes",
    "strict_5of5_genes",
    "shared_genes",
    "historical_only_genes",
    "strict_only_genes"
  ),
  value = c(
    length(historical_genes),
    length(strict_core_genes),
    length(shared),
    length(historical_only),
    length(strict_only)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  comparison,
  "results/week3_historical_vs_strict_5of5.csv",
  row.names = FALSE
)

write.csv(
  comparison_summary,
  "results/week3_historical_vs_strict_5of5_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Historical genes inside the full strict stability analysis
# ------------------------------------------------------------

idx <- match(
  historical_genes,
  strict_stability$gene
)

historical_audit <- data.frame(
  gene = historical_genes,
  present_in_strict_selected_312 =
    !is.na(idx),
  strict_outer_folds_selected =
    strict_stability$outer_folds_selected[idx],
  strict_selection_fraction =
    strict_stability$selection_fraction[idx],
  strict_stability_tier =
    strict_stability$stability_tier[idx],
  present_in_validated_input =
    historical_genes %in% rownames(counts),
  stringsAsFactors = FALSE
)

historical_audit$strict_outer_folds_selected[
  is.na(historical_audit$strict_outer_folds_selected)
] <- 0L

historical_audit$strict_selection_fraction[
  is.na(historical_audit$strict_selection_fraction)
] <- 0

historical_audit$strict_stability_tier[
  is.na(historical_audit$strict_stability_tier)
] <- "not_selected"

historical_audit_summary <- data.frame(
  metric = c(
    "historical_genes",
    "present_in_validated_input",
    "selected_at_least_1of5",
    "selected_at_least_2of5",
    "selected_at_least_3of5",
    "selected_at_least_4of5",
    "selected_5of5"
  ),
  value = c(
    nrow(historical_audit),
    sum(historical_audit$present_in_validated_input),
    sum(historical_audit$strict_outer_folds_selected >= 1),
    sum(historical_audit$strict_outer_folds_selected >= 2),
    sum(historical_audit$strict_outer_folds_selected >= 3),
    sum(historical_audit$strict_outer_folds_selected >= 4),
    sum(historical_audit$strict_outer_folds_selected == 5)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  historical_audit,
  "results/week3_historical_13_genes_in_strict_pipeline.csv",
  row.names = FALSE
)

write.csv(
  historical_audit_summary,
  "results/week3_historical_13_genes_in_strict_pipeline_summary.csv",
  row.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("WEEK 3 HISTORICAL VS STRICT STABILITY AUDIT\n")
cat("============================================================\n")

cat("Historical source:", historical_file, "\n")
cat("Historical 5/5 genes:", length(historical_genes), "\n")
cat("Strict 5/5 genes:", length(strict_core_genes), "\n")
cat("Shared 5/5 genes:", length(shared), "\n")

cat(
  "Historical genes present in validated input:",
  sum(historical_audit$present_in_validated_input),
  "/",
  nrow(historical_audit),
  "\n"
)

cat(
  "Historical genes selected >=1/5 strict folds:",
  sum(historical_audit$strict_outer_folds_selected >= 1),
  "\n"
)

cat(
  "Historical genes selected 5/5 strict folds:",
  sum(historical_audit$strict_outer_folds_selected == 5),
  "\n"
)

cat("\nIMPORTANT:\n")
cat(
  "This is a reproducibility audit of two different modelling procedures.\n"
)
cat(
  "Zero overlap does not demonstrate biological irrelevance of the historical genes.\n"
)
cat(
  "It demonstrates that the historical stable set is not reproduced by the strict nested-CV selection procedure.\n"
)
cat(
  "Historical full-dataset RF importance is exploratory and is not used to define the strict consensus candidates.\n"
)

cat("\nAudit outputs written successfully.\n")
