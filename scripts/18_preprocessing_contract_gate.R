#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("P1 QA GATE — U08 / P1-04 PREPROCESSING\n")
cat("============================================\n\n")

contract_file <- "docs/preprocessing_contract.md"
week3_runner <- "scripts/run_week3.R"
week3_impl <- "scripts/06_nested_cv_strict.R"
week4_runner <- "scripts/run_week4.R"
week4_rep_doc <- "docs/week4_expression_representation_decision.md"
week4_locked_doc <- "docs/week4_locked_external_validation.md"

required_files <- c(
  contract_file,
  week3_runner,
  week3_impl,
  week4_runner,
  week4_rep_doc,
  week4_locked_doc,
  "data_processed/GSE160638_raw_counts_PD1_aligned.rds",
  "data_processed/GSE160638_TableS2A_PD1_clinical.csv",
  "results/week3_gene_selection_frequency_strict.csv",
  "results/week4_expression_representation_decision_strict.csv",
  "results/week4_deployment_freeze_manifest_strict.csv",
  "results/week4_external_validation_metrics_strict.csv"
)

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

contains_regex <- function(text, pattern) {
  grepl(pattern, text, ignore.case = TRUE, perl = TRUE)
}

checks <- data.frame(
  check_id = character(),
  criterion = character(),
  description = character(),
  status = character(),
  stringsAsFactors = FALSE
)

add_check <- function(check_id, criterion, description, condition) {
  checks <<- rbind(
    checks,
    data.frame(
      check_id = check_id,
      criterion = criterion,
      description = description,
      status = if (isTRUE(condition)) "PASS" else "FAIL",
      stringsAsFactors = FALSE
    )
  )
}

# ------------------------------------------------------------
# 1. Required artifact existence
# ------------------------------------------------------------

for (f in required_files) {
  add_check(
    paste0("FILE_", basename(f)),
    "U08/P1-04",
    paste("Required artifact exists:", f),
    file.exists(f)
  )
}

missing_required <- required_files[!file.exists(required_files)]

if (length(missing_required) > 0) {
  cat("Missing required files:\n")
  cat(paste0(" - ", missing_required, collapse = "\n"), "\n\n")
}

# Stop textual checks only if their source file is unavailable.
contract_text <- if (file.exists(contract_file)) read_text(contract_file) else ""
week3_runner_text <- if (file.exists(week3_runner)) read_text(week3_runner) else ""
week3_text <- if (file.exists(week3_impl)) read_text(week3_impl) else ""
week4_runner_text <- if (file.exists(week4_runner)) read_text(week4_runner) else ""
week4_rep_text <- if (file.exists(week4_rep_doc)) read_text(week4_rep_doc) else ""
week4_locked_text <- if (file.exists(week4_locked_doc)) read_text(week4_locked_doc) else ""

# ------------------------------------------------------------
# 2. Contract structure
# ------------------------------------------------------------

add_check(
  "CONTRACT_U08",
  "U08",
  "Contract explicitly references U08 preprocessing",
  contains_regex(contract_text, "U08")
)

add_check(
  "CONTRACT_P1_04",
  "P1-04",
  "Contract explicitly references P1-04 transcriptomic preprocessing",
  contains_regex(contract_text, "P1-04")
)

add_check(
  "CONTRACT_EXPLORATORY_BOUNDARY",
  "U08/P1-04",
  "Week 2 exploratory preprocessing is separated from predictive preprocessing",
  contains_regex(contract_text, "Week 2 exploratory preprocessing")
)

add_check(
  "CONTRACT_TRAINING_ONLY",
  "U08/P1-04",
  "Predictive preprocessing is documented as training-data restricted",
  contains_regex(contract_text, "restricted to training data|estimated from training data")
)

add_check(
  "CONTRACT_OUTER_TEST",
  "U07/U08/P1-04",
  "Outer-test samples are excluded from preprocessing estimation",
  contains_regex(contract_text, "outer-test samples.*excluded")
)

add_check(
  "CONTRACT_STABILITY_BOUNDARY",
  "P1-09",
  "Stability analysis is separated from independent validation",
  contains_regex(contract_text, "post-validation stability analysis")
)

add_check(
  "CONTRACT_SAMPLEWISE_RANK",
  "U08/P1-04",
  "Frozen external representation is documented as samplewise_rank",
  contains_regex(contract_text, "samplewise_rank")
)

add_check(
  "CONTRACT_UNBLINDING",
  "U07/U08",
  "Post-unblinding modifications are explicitly prohibited",
  contains_regex(contract_text, "After unblinding")
)

add_check(
  "CONTRACT_WEEK3_RUNNER",
  "U13",
  "Week 3 reproduction command is documented",
  contains_regex(contract_text, "scripts/run_week3\\.R")
)

add_check(
  "CONTRACT_WEEK4_RUNNER",
  "U13",
  "Week 4 reproduction command is documented",
  contains_regex(contract_text, "scripts/run_week4\\.R")
)

# ------------------------------------------------------------
# 3. Week 3 implementation evidence
# ------------------------------------------------------------

add_check(
  "W3_RUNNER_STRICT_IMPLEMENTATION",
  "U13/P1-04",
  "Week 3 runner executes the strict nested-CV implementation",
  contains_regex(
    week3_runner_text,
    "scripts/06_nested_cv_strict\\.R"
  )
)

add_check(
  "W3_FILTERING_IMPLEMENTED",
  "P1-04",
  "Week 3 strict implementation filters genes inside training partitions",
  contains_regex(week3_text, "filterByExpr") &&
    contains_regex(week3_text, "select_genes_from_training")
)

add_check(
  "W3_TMM_SELECTION_TRAINING_ONLY",
  "P1-04/P1-05",
  "TMM normalization used for feature ranking is fitted inside training selection",
  contains_regex(week3_text, "calcNormFactors") &&
    contains_regex(week3_text, 'method\\s*=\\s*"TMM"')
)

add_check(
  "W3_SAMPLEWISE_TRANSFORMATION",
  "U08/P1-04",
  "Predictive expression transformation is sample-wise logCPM",
  contains_regex(week3_text, "samplewise_logcpm") &&
    contains_regex(week3_text, "colSums") &&
    contains_regex(week3_text, "log2")
)

add_check(
  "W3_HELDOUT_TRANSFORMATION_INDEPENDENT",
  "U07/U08/P1-04",
  "Held-out expression transformation does not estimate cross-sample parameters",
  contains_regex(
    week3_text,
    "does not estimate parameters from held-out samples"
  ) &&
    contains_regex(week3_text, "samplewise_logcpm")
)

add_check(
  "W3_FEATURE_SELECTION_IMPLEMENTED",
  "P1-05",
  "Week 3 code contains partition-specific selected-gene logic",
  contains_regex(week3_text, "selected_genes")
)

add_check(
  "W3_NESTED_STRUCTURE",
  "P1-05/P1-06",
  "Week 3 implementation contains outer/inner CV structure",
  contains_regex(week3_text, "outer") &&
    contains_regex(week3_text, "inner")
)

# ------------------------------------------------------------
# 4. Week 4 representation and freeze evidence
# ------------------------------------------------------------

add_check(
  "W4_REP_SAMPLEWISE_RANK",
  "U08/P1-04",
  "Week 4 representation decision records samplewise_rank",
  contains_regex(week4_rep_text, "samplewise_rank")
)

add_check(
  "W4_REP_LABEL_FREE",
  "U07/U08",
  "Week 4 representation decision documents label-free selection",
  contains_regex(
    week4_rep_text,
    "No external response labels|response labels.*FALSE|labels.*not"
  )
)

add_check(
  "W4_LOCKED_NO_RESELECTION",
  "U07/P1-05",
  "Locked external validation prohibits gene reselection",
  contains_regex(week4_locked_text, "reselect genes|gene reselection")
)

add_check(
  "W4_LOCKED_NO_REFIT",
  "U07/P1-06",
  "Locked external validation prohibits model refitting",
  contains_regex(week4_locked_text, "model refit|refit")
)

add_check(
  "W4_LOCKED_THRESHOLD",
  "U07/P1-06",
  "Locked external validation prohibits external threshold optimization",
  contains_regex(
    week4_locked_text,
    "threshold.*external|optimize.*threshold|threshold optimization"
  )
)

add_check(
  "W4_RUNNER_PRESENT",
  "U13",
  "Week 4 reproducibility runner is present and non-empty",
  file.exists(week4_runner) &&
    file.info(week4_runner)$size > 0
)

# ------------------------------------------------------------
# 5. Persist QA evidence
# ------------------------------------------------------------

dir.create("results", showWarnings = FALSE, recursive = TRUE)

output_file <- "results/qa_preprocessing_contract_gate.csv"

write.csv(
  checks,
  output_file,
  row.names = FALSE,
  quote = TRUE
)

cat("QA checks:\n\n")
print(checks, row.names = FALSE)

cat("\n--------------------------------------------\n")
cat("PASS:", sum(checks$status == "PASS"), "\n")
cat("FAIL:", sum(checks$status == "FAIL"), "\n")
cat("Output:", output_file, "\n")
cat("--------------------------------------------\n\n")

if (any(checks$status != "PASS")) {
  cat("U08 / P1-04 PREPROCESSING CONTRACT: FAIL\n")
  quit(status = 1)
}

cat("U08 / P1-04 PREPROCESSING CONTRACT: PASS\n")
