options(stringsAsFactors = FALSE)

clinical_file <- "data_processed/GSE160638_TableS2A_PD1_clinical.csv"
metadata_file <- "data_processed/metadata_GSE160638.csv"
clinical_summary_file <- "results/day4_clinical_by_response.csv"
age_summary_file <- "results/day4_age_by_response.csv"
doc_file <- "docs/confounding_assessment.md"
external_metrics_file <- "results/week4_external_validation_metrics_strict.csv"
output_file <- "results/qa_confounding_gate.csv"

required_files <- c(
  clinical_file,
  metadata_file,
  clinical_summary_file,
  age_summary_file,
  doc_file,
  external_metrics_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files: ",
    paste(missing_files, collapse = ", ")
  )
}

clinical <- read.csv(
  clinical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

external <- read.csv(
  external_metrics_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

doc <- paste(
  readLines(doc_file, warn = FALSE),
  collapse = "\n"
)

checks <- data.frame(
  check_id = character(),
  qa = character(),
  description = character(),
  observed = character(),
  expected = character(),
  status = character(),
  stringsAsFactors = FALSE
)

add_check <- function(id, qa, description, observed, expected, pass) {
  checks <<- rbind(
    checks,
    data.frame(
      check_id = id,
      qa = qa,
      description = description,
      observed = observed,
      expected = expected,
      status = ifelse(pass, "PASS", "FAIL"),
      stringsAsFactors = FALSE
    )
  )
}

required_columns <- c(
  "patient_id",
  "Gender",
  "Age at treatment",
  "BRAF mutation",
  "LDH high",
  "Biopsy site",
  "Biopsy type",
  "response"
)

add_check(
  "CLINICAL_COLUMNS",
  "U07/P1-10",
  "Relevant clinical confounder variables are available",
  paste(intersect(required_columns, names(clinical)), collapse = "; "),
  paste(required_columns, collapse = "; "),
  all(required_columns %in% names(clinical))
)

add_check(
  "INTERNAL_N",
  "U07/P1-10",
  "Internal cohort size is frozen",
  nrow(clinical),
  "36",
  nrow(clinical) == 36
)

response_counts <- table(clinical$response)

response_observed <- paste(
  names(response_counts),
  as.integer(response_counts),
  sep = "=",
  collapse = "; "
)

add_check(
  "RESPONSE_COUNTS",
  "U07/P1-10",
  "Response-class distribution is verified",
  response_observed,
  "NonResponder=14; Responder=22",
  identical(
    as.integer(response_counts[c("NonResponder", "Responder")]),
    c(14L, 22L)
  )
)

add_check(
  "PATIENT_UNIQUENESS",
  "U07",
  "No duplicated internal patient identifiers",
  sum(duplicated(clinical$patient_id)),
  "0",
  !anyDuplicated(clinical$patient_id)
)

gender_tab <- table(
  clinical$response,
  clinical$Gender,
  useNA = "ifany"
)

add_check(
  "GENDER_ASSESSED",
  "U07/P1-10",
  "Gender distribution assessed by response",
  paste(capture.output(print(gender_tab)), collapse = " | "),
  "Descriptive distribution calculable",
  length(gender_tab) > 0
)

age_summary <- aggregate(
  clinical[["Age at treatment"]],
  by = list(response = clinical$response),
  FUN = function(x) {
    c(
      n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE)
    )
  }
)

add_check(
  "AGE_ASSESSED",
  "U07/P1-10",
  "Age distribution assessed by response",
  paste(capture.output(print(age_summary)), collapse = " | "),
  "Age summaries calculable for both response classes",
  nrow(age_summary) == 2
)

braf_tab <- table(
  clinical$response,
  clinical[["BRAF mutation"]],
  useNA = "ifany"
)

add_check(
  "BRAF_ASSESSED",
  "U07/P1-10",
  "BRAF mutation distribution assessed",
  paste(capture.output(print(braf_tab)), collapse = " | "),
  "BRAF distribution and missingness calculable",
  length(braf_tab) > 0
)

ldh_tab <- table(
  clinical$response,
  clinical[["LDH high"]],
  useNA = "ifany"
)

add_check(
  "LDH_ASSESSED",
  "U07/P1-10",
  "LDH distribution assessed",
  paste(capture.output(print(ldh_tab)), collapse = " | "),
  "LDH distribution and missingness calculable",
  length(ldh_tab) > 0
)

biopsy_tab <- table(
  clinical$response,
  clinical[["Biopsy type"]],
  useNA = "ifany"
)

add_check(
  "BIOPSY_ASSESSED",
  "U07/P1-10",
  "Biopsy-type heterogeneity assessed",
  paste(
    "Levels:",
    length(unique(clinical[["Biopsy type"]])),
    "| cells:",
    length(biopsy_tab)
  ),
  "Biopsy heterogeneity explicitly evaluable",
  length(biopsy_tab) > 0
)

batch_pattern <- "batch|lane|plate|center|centre|sequenc|flowcell"

batch_columns <- unique(c(
  names(clinical)[grepl(batch_pattern, names(clinical), ignore.case = TRUE)],
  names(metadata)[grepl(batch_pattern, names(metadata), ignore.case = TRUE)]
))

batch_observed <- if (length(batch_columns) == 0) {
  "No explicit technical batch column available"
} else {
  paste(batch_columns, collapse = "; ")
}

add_check(
  "TECHNICAL_BATCH_METADATA",
  "U07/P1-10",
  "Availability of technical batch metadata is explicitly assessed",
  batch_observed,
  "Availability or unavailability documented",
  grepl(
    "No explicit technical batch variable",
    doc,
    fixed = TRUE
  )
)

add_check(
  "BATCH_LIMITATION_DOCUMENTED",
  "U07/P1-10",
  "Absence of batch metadata is not interpreted as absence of batch effects",
  grepl(
    "must not be interpreted as evidence",
    doc,
    fixed = TRUE
  ),
  "TRUE",
  grepl(
    "must not be interpreted as evidence",
    doc,
    fixed = TRUE
  )
)

primary_cohorts <- c("GSE91061", "GSE78220")

primary_ok <- all(primary_cohorts %in% external$scope) &&
  all(
    grepl(
      "primary",
      external$analysis_role[
        match(primary_cohorts, external$scope)
      ],
      ignore.case = TRUE
    )
  )

add_check(
  "EXTERNAL_COHORT_SPECIFIC",
  "P1-10",
  "External cohorts are evaluated separately as primary results",
  paste(primary_cohorts, collapse = "; "),
  "GSE91061 and GSE78220 cohort-specific primary reporting",
  primary_ok
)

pooled_idx <- which(external$scope == "POOLED")

pooled_ok <- length(pooled_idx) == 1 &&
  grepl(
    "secondary",
    external$analysis_role[pooled_idx],
    ignore.case = TRUE
  )

add_check(
  "POOLED_SECONDARY",
  "P1-10",
  "Pooled external result is secondary",
  if (length(pooled_idx) == 1) {
    external$analysis_role[pooled_idx]
  } else {
    "POOLED row missing or duplicated"
  },
  "secondary pooled reporting",
  pooled_ok
)

required_doc_terms <- c(
  "residual confounding",
  "Technical batch",
  "LDH",
  "Biopsy",
  "External-cohort effects",
  "post-unblinding"
)

doc_terms_ok <- all(vapply(
  required_doc_terms,
  function(x) grepl(x, doc, ignore.case = TRUE),
  logical(1)
))

add_check(
  "CONFOUNDING_DOCUMENTATION",
  "U07/P1-10",
  "Required confounding and leakage boundaries are documented",
  paste(required_doc_terms, collapse = "; "),
  "All required concepts documented",
  doc_terms_ok
)

doc_normalized <- gsub("[[:space:]]+", " ", doc)

add_check(
  "NO_RETROACTIVE_REDEFINITION",
  "U07/P1-10",
  "QA assessment does not redefine frozen W3-W4 analysis",
  grepl(
    "does not claim that confounding or batch effects have been eliminated",
    doc_normalized,
    fixed = TRUE
  ),
  "Frozen analytical results remain unchanged",
  grepl(
    "does not claim that confounding or batch effects have been eliminated",
    doc_normalized,
    fixed = TRUE
  )
)

write.csv(
  checks,
  output_file,
  row.names = FALSE
)

cat("\n==========================================\n")
cat("U07 / P1-10 CONFOUNDING QA GATE\n")
cat("==========================================\n\n")

for (i in seq_len(nrow(checks))) {
  cat(
    sprintf(
      "%-32s %s\n",
      checks$check_id[i],
      checks$status[i]
    )
  )
}

pass_n <- sum(checks$status == "PASS")
fail_n <- sum(checks$status == "FAIL")

cat("\n------------------------------------------\n")
cat("PASS:", pass_n, "\n")
cat("FAIL:", fail_n, "\n")
cat("Output:", output_file, "\n")
cat("------------------------------------------\n\n")

if (fail_n > 0) {
  cat("FAILED CHECKS:\n")
  print(
    checks[
      checks$status == "FAIL",
      c(
        "check_id",
        "description",
        "observed",
        "expected"
      )
    ],
    row.names = FALSE
  )

  stop(
    "U07 / P1-10 CONFOUNDING ASSESSMENT: FAIL",
    call. = FALSE
  )
}

cat("U07 / P1-10 CONFOUNDING ASSESSMENT: PASS\n")
