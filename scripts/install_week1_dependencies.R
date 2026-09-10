# ============================================================
# Install Week 1 dependencies
# ============================================================
#
# Tested environment:
#
# R        4.5.1
# GEOquery 2.76.0
# dplyr    1.2.0
# readxl   1.4.5
#
# This script installs missing dependencies.
# It does not modify the analytical workflow.
# ============================================================


cran_packages <- c(
  "dplyr",
  "readxl"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0L) {

  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}


if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


if (!requireNamespace("GEOquery", quietly = TRUE)) {

  BiocManager::install(
    "GEOquery",
    ask = FALSE,
    update = FALSE
  )
}


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