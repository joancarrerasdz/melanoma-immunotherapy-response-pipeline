# ============================================================
# Restore Week 1 dependencies from renv.lock
# ============================================================
#
# Tested environment:
#
# R        4.5.1
# GEOquery 2.76.0
# dplyr    1.2.0
# readxl   1.4.5
#
# This script restores the exact locked project environment.
# It does not update renv.lock or the analytical workflow.
# ============================================================


if (!file.exists("renv.lock")) {
  stop(
    "renv.lock was not found. Run this script from the repository root.",
    call. = FALSE
  )
}

if (!requireNamespace("renv", quietly = TRUE)) {
  stop(
    "renv is unavailable. Open the project once so renv/activate.R can bootstrap it.",
    call. = FALSE
  )
}

cat("Restoring the project library from renv.lock...\n")
renv::restore(prompt = FALSE)

activate_path <- file.path("renv", "activate.R")
if (file.exists(activate_path)) {
  activate_lines <- readLines(activate_path, warn = FALSE)
  writeLines(
    sub("[[:blank:]]+$", "", activate_lines),
    activate_path,
    useBytes = TRUE
  )
}

cat("Validating the restored project environment...\n")
renv::status()

required_packages <- c(
  "GEOquery",
  "dplyr",
  "readxl"
)

still_missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(still_missing) > 0L) {
  stop(
    "Dependencies still missing after installation: ",
    paste(still_missing, collapse = ", "),
    call. = FALSE
  )
}


cat("\nWEEK 1 DEPENDENCIES AVAILABLE\n\n")

for (pkg in required_packages) {
  cat(
    sprintf(
      "%-10s %s\n",
      pkg,
      as.character(packageVersion(pkg))
    )
  )
}

cat("\n")
cat(R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")