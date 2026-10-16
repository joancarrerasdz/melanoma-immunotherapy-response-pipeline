#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1
)

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

required_root_files <- c(
  "renv.lock",
  ".Rprofile",
  "scripts"
)

missing_root_files <- required_root_files[
  !file.exists(required_root_files)
]

if (length(missing_root_files) > 0L) {
  stop(
    paste(
      "Run this pipeline from the repository root.",
      "Missing:",
      paste(missing_root_files, collapse = ", ")
    ),
    call. = FALSE
  )
}

stages <- c(
  week1 = "scripts/run_week1.R",
  week2 = "scripts/run_week2.R",
  week3 = "scripts/run_week3.R",
  week4_locked_validation = "scripts/run_week4.R",
  preprocessing_contract = "scripts/18_preprocessing_contract_gate.R",
  reproducible_environment = "scripts/19_reproducible_environment_gate.R",
  confounding_gate = "scripts/20_confounding_gate.R",
  public_figure = "scripts/21_public_figure_polish.R"
)

missing_scripts <- stages[!file.exists(stages)]

if (length(missing_scripts) > 0L) {
  stop(
    paste(
      "Pipeline scripts are missing:",
      paste(missing_scripts, collapse = ", ")
    ),
    call. = FALSE
  )
}

rscript <- Sys.which("Rscript")

if (!nzchar(rscript)) {
  stop("Rscript was not found on PATH.", call. = FALSE)
}

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) == 0L || identical(arguments, "all")) {
  selected_stages <- stages
} else {
  unknown_stages <- setdiff(arguments, names(stages))

  if (length(unknown_stages) > 0L) {
    stop(
      paste(
        "Unknown pipeline stages:",
        paste(unknown_stages, collapse = ", "),
        "\nAvailable stages:",
        paste(names(stages), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  selected_stages <- stages[arguments]
}

cat("\n")
cat("============================================================\n")
cat("P1 REPRODUCIBLE PIPELINE\n")
cat("============================================================\n")
cat("Project root:", project_root, "\n")
cat("R executable:", rscript, "\n")
cat("Stages:", paste(names(selected_stages), collapse = ", "), "\n")
cat("Started:", format(Sys.time(), tz = "UTC"), "UTC\n")

pipeline_started <- Sys.time()

for (stage_name in names(selected_stages)) {
  script <- unname(selected_stages[[stage_name]])
  stage_started <- Sys.time()

  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("STAGE:", stage_name, "\n")
  cat("SCRIPT:", script, "\n")
  cat("------------------------------------------------------------\n")

  status <- system2(
    command = rscript,
    args = script,
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop(
      paste(
        "Pipeline stage failed:",
        stage_name,
        "with exit status",
        status
      ),
      call. = FALSE
    )
  }

  elapsed <- as.numeric(
    difftime(Sys.time(), stage_started, units = "secs")
  )

  cat(
    sprintf(
      "PASS: %s completed in %.1f seconds\n",
      stage_name,
      elapsed
    )
  )
}

total_elapsed <- as.numeric(
  difftime(Sys.time(), pipeline_started, units = "mins")
)

cat("\n")
cat("============================================================\n")
cat("P1 REPRODUCIBLE PIPELINE: PASS\n")
cat(sprintf("Completed in %.2f minutes\n", total_elapsed))
cat("============================================================\n")