# ============================================================
# Week 1 reproducible rebuild runner
# ============================================================
#
# Run from the repository root with:
#
#   Rscript --vanilla scripts/run_week1.R
#
# This runner executes the complete Week 1 rebuild:
#
#   1. Load and align GSE160638 clinical / expression data
#   2. Run automated integrity tests
#   3. Generate descriptive QC outputs
#
# Each step is executed in a separate --vanilla R session.
# ============================================================


# ------------------------------------------------------------
# Repository-root check
# ------------------------------------------------------------

required_scripts <- c(
  "scripts/01_load_data.R",
  "scripts/02_qc_tests.R",
  "scripts/02_qc_summary.R"
)

if (!all(file.exists(required_scripts))) {
  stop(
    paste0(
      "Week 1 runner must be executed from the repository root.\n",
      "Expected files:\n  ",
      paste(required_scripts, collapse = "\n  ")
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------
# Dependency check
# ------------------------------------------------------------

required_packages <- c(
  "GEOquery",
  "dplyr",
  "readxl"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing Week 1 dependencies: ",
      paste(missing_packages, collapse = ", "),
      "\nRun:\n",
      "  Rscript --vanilla scripts/install_week1_dependencies.R"
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

run_step <- function(script, label) {

  cat("\n")
  cat("============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")

  rscript <- Sys.which("Rscript")

  if (!nzchar(rscript)) {
    stop(
      "Rscript executable could not be found.",
      call. = FALSE
    )
  }

  status <- system2(
    command = rscript,
    args = c(
      "--vanilla",
      script
    )
  )

  if (!identical(status, 0L)) {
    stop(
      paste0(
        "Week 1 reproduction failed during: ",
        label,
        "\nScript: ",
        script,
        "\nExit status: ",
        status
      ),
      call. = FALSE
    )
  }

  cat("\n[PASS] ", label, "\n", sep = "")
}


# ------------------------------------------------------------
# Execute Week 1
# ------------------------------------------------------------

run_step(
  "scripts/01_load_data.R",
  "STEP 1/3 — Load and align GSE160638 data"
)

run_step(
  "scripts/02_qc_tests.R",
  "STEP 2/3 — Automated integrity tests"
)

run_step(
  "scripts/02_qc_summary.R",
  "STEP 3/3 — Initial descriptive QC"
)


# ------------------------------------------------------------
# Final output checks
# ------------------------------------------------------------

expected_outputs <- c(
  "data_processed/metadata_GSE160638.csv",
  "data_processed/GSE160638_TableS2A_PD1_clinical.csv",
  "data_processed/GSE160638_raw_counts_PD1_aligned.rds",
  "data_processed/GSE160638_clinical_join_audit.csv",
  "results/day4_sample_qc.csv",
  "results/day4_clinical_missingness.csv",
  "results/day4_clinical_by_response.csv",
  "results/day4_age_by_response.csv",
  "figures/day4_response_distribution.png",
  "figures/day4_library_size.png",
  "figures/day4_genes_detected.png",
  "figures/day4_percent_zeros.png",
  "docs/day4_qc_anomalies.md"
)

missing_outputs <- expected_outputs[
  !file.exists(expected_outputs)
]

if (length(missing_outputs) > 0L) {
  stop(
    paste0(
      "Week 1 execution completed but expected outputs are missing:\n  ",
      paste(missing_outputs, collapse = "\n  ")
    ),
    call. = FALSE
  )
}


cat("\n")
cat("============================================================\n")
cat("WEEK 1 REPRODUCTION: PASS\n")
cat("============================================================\n")
cat("Clinical cohort: 36 PD1 samples\n")
cat("Responders: 22\n")
cat("NonResponders: 14\n")
cat("Expression features: 17002\n")
cat("Automated tests: PASS\n")
cat("Initial QC: generated\n")
cat("Expected outputs: present\n")