# ============================================================
# Day 4 — Automated clinical / expression integrity tests
# ============================================================

clinical_file <- file.path(
  "data_processed",
  "GSE160638_TableS2A_PD1_clinical.csv"
)

metadata_file <- file.path(
  "data_processed",
  "metadata_GSE160638.csv"
)

audit_file <- file.path(
  "data_processed",
  "GSE160638_clinical_join_audit.csv"
)

counts_file <- file.path(
  "data_processed",
  "GSE160638_raw_counts_PD1_aligned.rds"
)

required_files <- c(
  clinical_file,
  metadata_file,
  audit_file,
  counts_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing required Day 3 deliverables: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

clinical <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

audit <- read.csv(
  audit_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- readRDS(counts_file)


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

fail_test <- function(label, detail = NULL) {
  msg <- paste0("[FAIL] ", label)

  if (!is.null(detail)) {
    msg <- paste0(msg, ": ", detail)
  }

  stop(msg, call. = FALSE)
}

pass_test <- function(label) {
  message("[PASS] ", label)
}

check_test <- function(condition, label, detail = NULL) {
  if (!isTRUE(condition)) {
    fail_test(label, detail)
  }

  pass_test(label)
}


# ------------------------------------------------------------
# Test 0 — Expected Table S2A-derived schema
# ------------------------------------------------------------

expected_clinical_columns <- c(
  "patient_id_raw",
  "patient_id",
  "TIL ID",
  "Cohort",
  "Gender",
  "Age at treatment",
  "Biopsy site",
  "Biopsy type",
  "BRAF mutation",
  "LDH high",
  "Response",
  "OR",
  "PFS T",
  "PFS E",
  "OS T",
  "OS E",
  "Prior MAPKi",
  "Prior Ipi",
  "Prior treatments",
  "Post treatments",
  "Immune cluster",
  "ImmuneScore",
  "response_recist",
  "response"
)

check_test(
  identical(names(clinical), expected_clinical_columns),
  "Table S2A clinical schema unchanged",
  paste(
    "Expected:",
    paste(expected_clinical_columns, collapse = " | "),
    "; observed:",
    paste(names(clinical), collapse = " | ")
  )
)


# ------------------------------------------------------------
# Test 1 — Exactly 36 unique PD1 identifiers
# ------------------------------------------------------------

clinical_ids <- clinical$patient_id

check_test(
  nrow(clinical) == 36L &&
    length(unique(clinical_ids)) == 36L &&
    !anyNA(clinical_ids) &&
    !any(trimws(clinical_ids) == ""),
  "Exactly 36 unique clinical identifiers"
)


# ------------------------------------------------------------
# Test 2 — No missing response labels
# ------------------------------------------------------------

response_columns <- c(
  "Response",
  "response_recist",
  "response"
)

response_missing <- vapply(
  clinical[response_columns],
  function(x) {
    any(is.na(x) | trimws(as.character(x)) == "")
  },
  logical(1)
)

check_test(
  !any(response_missing),
  "No missing response labels",
  paste(
    names(response_missing)[response_missing],
    collapse = ", "
  )
)


# ------------------------------------------------------------
# Test 3 — Only CR / PR / SD / PD endpoint values
# ------------------------------------------------------------

allowed_recist <- c("CR", "PR", "SD", "PD")

observed_recist <- sort(unique(clinical$Response))

check_test(
  all(observed_recist %in% allowed_recist),
  "Only CR/PR/SD/PD response values",
  paste(
    "Observed:",
    paste(observed_recist, collapse = ", ")
  )
)

check_test(
  identical(clinical$Response, clinical$response_recist),
  "Original and normalized RECIST labels agree"
)


# ------------------------------------------------------------
# Test 4 — Exactly 22 responders and 14 non-responders
# ------------------------------------------------------------

expected_binary_response <- ifelse(
  clinical$response_recist %in% c("CR", "PR"),
  "Responder",
  "NonResponder"
)

check_test(
  identical(clinical$response, expected_binary_response),
  "Binary response mapping agrees with RECIST"
)

response_counts <- table(clinical$response)

check_test(
  identical(
    as.integer(
      response_counts[c("NonResponder", "Responder")]
    ),
    c(14L, 22L)
  ),
  "Exactly 22 responders and 14 non-responders"
)

# ------------------------------------------------------------
# Test 5 — No TIL samples in primary PD1 cohort
# ------------------------------------------------------------

check_test(
  all(grepl("^PD1_[0-9]+$", clinical$patient_id)) &&
    all(grepl("^PD1_[0-9]+$", meta$sample_name)) &&
    all(grepl("^PD1_[0-9]+$", colnames(counts))) &&
    !any(grepl("^TIL_", clinical$patient_id)) &&
    !any(grepl("^TIL_", meta$sample_name)) &&
    !any(grepl("^TIL_", colnames(counts))),
  "Primary cohort contains PD1 samples only"
)


# ------------------------------------------------------------
# Test 6 — Complete clinical / GEO / counts correspondence
# ------------------------------------------------------------

check_test(
  nrow(meta) == 36L &&
    ncol(counts) == 36L &&
    nrow(audit) == 36L,
  "Clinical, GEO and counts contain 36 PD1 samples"
)

check_test(
  identical(colnames(counts), clinical$patient_id),
  "Counts and clinical identifiers have identical order"
)

check_test(
  identical(colnames(counts), meta$sample_name),
  "Counts and GEO metadata identifiers have identical order"
)

check_test(
  identical(colnames(counts), audit$patient_id),
  "Counts and audit identifiers have identical order"
)

check_test(
  all(audit$in_counts) &&
    all(audit$in_geo_metadata) &&
    all(audit$in_table_s2a) &&
    all(audit$status == "matched"),
  "Clinical join audit contains no pending correspondence"
)


# ------------------------------------------------------------
# Final Day 4 automated-test gate
# ------------------------------------------------------------

message("")
message("============================================")
message("DAY 4 AUTOMATED TESTS: PASS")
message("============================================")
message(
  "Samples: 36 PD1 | Responders: 22 | NonResponders: 14"
)
message(
  "Clinical/GEO/counts correspondence: complete"
)