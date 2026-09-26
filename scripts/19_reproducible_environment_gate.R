options(stringsAsFactors = FALSE)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

checks <- list()

add_check <- function(id,
                      criterion,
                      description,
                      pass,
                      observed = "",
                      expected = "") {

  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = id,
    criterion = criterion,
    description = description,
    observed = as.character(observed),
    expected = as.character(expected),
    status = if (isTRUE(pass)) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}

required_files <- c(
  ".Rprofile",
  "DESCRIPTION",
  "renv.lock",
  "renv/activate.R",
  "renv/settings.json",
  "results/week4_deployment_environment_strict.csv"
)

for (f in required_files) {

  add_check(
    id = paste0(
      "FILE_",
      gsub("[^A-Za-z0-9]+", "_", f)
    ),
    criterion = "U11",
    description = paste(
      "Required reproducibility artifact exists:",
      f
    ),
    pass = file.exists(f),
    observed = if (file.exists(f)) "present" else "missing",
    expected = "present"
  )
}

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is required to evaluate the U11 gate.")
}

lock <- renv::lockfile_read("renv.lock")

snapshot_type <- renv::settings$snapshot.type()

add_check(
  id = "RENV_EXPLICIT_SNAPSHOT",
  criterion = "U11",
  description = "renv uses explicit dependency snapshot mode",
  pass = identical(snapshot_type, "explicit"),
  observed = snapshot_type,
  expected = "explicit"
)

rprofile_text <- paste(
  readLines(".Rprofile", warn = FALSE),
  collapse = "\n"
)

add_check(
  id = "RPROFILE_RENV_ACTIVATION",
  criterion = "U11",
  description = ".Rprofile activates the project renv environment",
  pass = grepl(
    "renv/activate\\.R",
    rprofile_text
  ),
  observed = if (
    grepl("renv/activate\\.R", rprofile_text)
  ) {
    "renv activation present"
  } else {
    "renv activation not found"
  },
  expected = "renv activation present"
)

expected_imports <- sort(
  c(
    "AnnotationDbi",
    "caret",
    "edgeR",
    "org.Hs.eg.db",
    "pROC",
    "randomForest",
    "readxl"
  )
)

desc <- read.dcf("DESCRIPTION")

if (!"Imports" %in% colnames(desc)) {
  stop("DESCRIPTION does not contain an Imports field.")
}

imports_raw <- desc[1, "Imports"]

imports <- trimws(
  unlist(
    strsplit(
      gsub("[\r\n]", " ", imports_raw),
      ",",
      fixed = TRUE
    )
  )
)

imports <- sub(
  "\\s*\\(.*$",
  "",
  imports
)

imports <- sort(
  imports[nzchar(imports)]
)

add_check(
  id = "DESCRIPTION_DIRECT_DEPENDENCIES",
  criterion = "U11",
  description = paste(
    "DESCRIPTION declares the audited strict",
    "W3-W4 direct dependencies"
  ),
  pass = setequal(imports, expected_imports),
  observed = paste(imports, collapse = ";"),
  expected = paste(expected_imports, collapse = ";")
)

for (pkg in expected_imports) {

  package_record <- lock$Packages[[pkg]]

  present <- !is.null(package_record)

  observed_version <- if (present) {
    as.character(package_record$Version)
  } else {
    "missing"
  }

  add_check(
    id = paste0(
      "LOCK_PACKAGE_",
      gsub("[^A-Za-z0-9]+", "_", pkg)
    ),
    criterion = "U11",
    description = paste(
      "Direct dependency is recorded in renv.lock:",
      pkg
    ),
    pass = present,
    observed = observed_version,
    expected = "recorded"
  )
}

locked_r_version <- if (
  !is.null(lock$R$Version)
) {
  as.character(lock$R$Version)
} else {
  NA_character_
}

add_check(
  id = "LOCK_R_VERSION",
  criterion = "U11",
  description = "renv.lock records the frozen R version",
  pass = identical(locked_r_version, "4.5.1"),
  observed = locked_r_version,
  expected = "4.5.1"
)

locked_bioc_version <- if (
  !is.null(lock$Bioconductor$Version)
) {
  as.character(lock$Bioconductor$Version)
} else {
  NA_character_
}

add_check(
  id = "LOCK_BIOCONDUCTOR_VERSION",
  criterion = "U11",
  description = "renv.lock records Bioconductor 3.21",
  pass = identical(locked_bioc_version, "3.21"),
  observed = locked_bioc_version,
  expected = "3.21"
)

frozen_file <- "results/week4_deployment_environment_strict.csv"

frozen <- read.csv(
  frozen_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "component",
  "version"
)

if (!all(
  required_columns %in% colnames(frozen)
)) {
  stop(
    frozen_file,
    " must contain component and version columns."
  )
}

versions_equivalent <- function(x, y) {

  if (
    length(x) != 1L ||
    length(y) != 1L ||
    is.na(x) ||
    is.na(y)
  ) {
    return(FALSE)
  }

  normalize_version <- function(z) {
    gsub(
      "-",
      ".",
      trimws(as.character(z)),
      fixed = TRUE
    )
  }

  identical(
    normalize_version(x),
    normalize_version(y)
  )
}

for (i in seq_len(nrow(frozen))) {

  component <- as.character(
    frozen$component[i]
  )

  expected_version <- as.character(
    frozen$version[i]
  )

  observed_version <- if (
    identical(component, "R")
  ) {

    locked_r_version

  } else if (
    !is.null(lock$Packages[[component]])
  ) {

    as.character(
      lock$Packages[[component]]$Version
    )

  } else {

    NA_character_
  }

  add_check(
    id = paste0(
      "FROZEN_VERSION_",
      gsub(
        "[^A-Za-z0-9]+",
        "_",
        component
      )
    ),
    criterion = "U11",
    description = paste(
      "Lockfile matches frozen Week 4 provenance for",
      component
    ),
    pass = versions_equivalent(
      observed_version,
      expected_version
    ),
    observed = observed_version,
    expected = expected_version
  )
}

tracked_library <- tryCatch(
  system2(
    "git",
    c(
      "ls-files",
      "--",
      "renv/library"
    ),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) character()
)

tracked_library <- tracked_library[
  nzchar(tracked_library)
]

add_check(
  id = "RENV_LIBRARY_NOT_TRACKED",
  criterion = "U11",
  description = paste(
    "Project-local renv package library",
    "is not tracked by Git"
  ),
  pass = length(tracked_library) == 0L,
  observed = if (
    length(tracked_library) == 0L
  ) {
    "none"
  } else {
    paste(
      tracked_library,
      collapse = ";"
    )
  },
  expected = "none"
)

out <- do.call(
  rbind,
  checks
)

out_file <- file.path(
  results_dir,
  "qa_reproducible_environment_gate.csv"
)

write.csv(
  out,
  out_file,
  row.names = FALSE,
  quote = TRUE
)

cat("\n")
cat("============================================\n")
cat("U11 REPRODUCIBLE ENVIRONMENT GATE\n")
cat("============================================\n")
cat("\n")

for (i in seq_len(nrow(out))) {

  cat(
    sprintf(
      "%-42s %s\n",
      out$check_id[i],
      out$status[i]
    )
  )
}

n_pass <- sum(
  out$status == "PASS"
)

n_fail <- sum(
  out$status == "FAIL"
)

cat("\n")
cat("--------------------------------------------\n")
cat("PASS:", n_pass, "\n")
cat("FAIL:", n_fail, "\n")
cat("Output:", out_file, "\n")
cat("--------------------------------------------\n")
cat("\n")

if (n_fail > 0L) {

  cat("FAILED CHECKS:\n")

  print(
    out[
      out$status == "FAIL",
      c(
        "check_id",
        "description",
        "observed",
        "expected"
      )
    ],
    row.names = FALSE
  )

  quit(
    status = 1L
  )
}

cat("U11 REPRODUCIBLE ENVIRONMENT: PASS\n")
