# ============================================================
# Day 4 — Initial descriptive QC
# ============================================================
#
# This script performs descriptive QC only.
# It does NOT:
#   - remove samples
#   - impute clinical values
#   - harmonize clinical categories
#   - perform predictive inference
#
# Inputs are validated Day 3 deliverables.
# ============================================================


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

clinical_file <- file.path(
  "data_processed",
  "GSE160638_TableS2A_PD1_clinical.csv"
)

counts_file <- file.path(
  "data_processed",
  "GSE160638_raw_counts_PD1_aligned.rds"
)

if (!file.exists(clinical_file)) {
  stop(
    "Missing clinical input: ",
    clinical_file,
    call. = FALSE
  )
}

if (!file.exists(counts_file)) {
  stop(
    "Missing count matrix input: ",
    counts_file,
    call. = FALSE
  )
}


# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

for (dir in c("results", "figures", "docs")) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

clinical <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(
  readRDS(counts_file)
)


# ------------------------------------------------------------
# Basic integrity checks
# ------------------------------------------------------------

if (!is.numeric(counts)) {
  stop("Count matrix is not numeric.", call. = FALSE)
}

if (anyNA(counts)) {
  stop("Count matrix contains missing values.", call. = FALSE)
}

if (any(counts < 0)) {
  stop("Count matrix contains negative values.", call. = FALSE)
}

if (any(abs(counts - round(counts)) > 1e-8)) {
  stop(
    "Count matrix contains non-integer-like values.",
    call. = FALSE
  )
}

if (
  nrow(clinical) != 36L ||
  ncol(counts) != 36L ||
  !identical(colnames(counts), clinical$patient_id)
) {
  stop(
    "Clinical/count alignment is invalid for Day 4 QC.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# Sample-level expression QC
# ------------------------------------------------------------

library_size <- colSums(counts)

genes_detected <- colSums(counts > 0)

percent_zeros <- 100 * colMeans(counts == 0)


iqr_limits <- function(x) {
  q <- quantile(
    x,
    probs = c(0.25, 0.75),
    na.rm = TRUE,
    names = FALSE
  )

  iqr <- q[2] - q[1]

  c(
    lower = q[1] - 1.5 * iqr,
    upper = q[2] + 1.5 * iqr
  )
}


library_limits <- iqr_limits(
  log10(library_size)
)

gene_limits <- iqr_limits(
  genes_detected
)

zero_limits <- iqr_limits(
  percent_zeros
)


sample_qc <- data.frame(
  sample_name = colnames(counts),
  response = clinical$response,
  response_recist = clinical$response_recist,
  library_size = as.numeric(library_size),
  genes_detected = as.integer(genes_detected),
  percent_zeros = as.numeric(percent_zeros),
  flag_library_size = (
    log10(library_size) < library_limits["lower"] |
      log10(library_size) > library_limits["upper"]
  ),
  flag_genes_detected = (
    genes_detected < gene_limits["lower"] |
      genes_detected > gene_limits["upper"]
  ),
  flag_percent_zeros = (
    percent_zeros < zero_limits["lower"] |
      percent_zeros > zero_limits["upper"]
  ),
  stringsAsFactors = FALSE
)

sample_qc$potential_outlier <- (
  sample_qc$flag_library_size |
    sample_qc$flag_genes_detected |
    sample_qc$flag_percent_zeros
)


write.csv(
  sample_qc,
  file.path(
    "results",
    "day4_sample_qc.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Clinical missingness
# ------------------------------------------------------------

missing_n <- function(x) {
  sum(
    is.na(x) |
      trimws(as.character(x)) == ""
  )
}


clinical_qc_variables <- c(
  "Age at treatment",
  "Gender",
  "BRAF mutation",
  "LDH high",
  "Biopsy site",
  "Biopsy type"
)


clinical_missingness <- data.frame(
  variable = clinical_qc_variables,
  missing_n = vapply(
    clinical[clinical_qc_variables],
    missing_n,
    integer(1)
  ),
  total_n = nrow(clinical),
  stringsAsFactors = FALSE
)

clinical_missingness$missing_percent <- (
  100 *
    clinical_missingness$missing_n /
    clinical_missingness$total_n
)


write.csv(
  clinical_missingness,
  file.path(
    "results",
    "day4_clinical_missingness.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Age distribution by binary response
# ------------------------------------------------------------

age_by_response <- do.call(
  rbind,
  lapply(
    split(
      clinical[["Age at treatment"]],
      clinical$response
    ),
    function(x) {
      data.frame(
        n = sum(!is.na(x)),
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE)
      )
    }
  )
)

age_by_response$response <- rownames(age_by_response)
rownames(age_by_response) <- NULL

age_by_response <- age_by_response[
  ,
  c(
    "response",
    "n",
    "mean",
    "median",
    "min",
    "max"
  )
]

write.csv(
  age_by_response,
  file.path(
    "results",
    "day4_age_by_response.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Clinical distributions by binary response
# ------------------------------------------------------------

categorical_variables <- c(
  "Gender",
  "BRAF mutation",
  "LDH high",
  "Biopsy type"
)


clinical_by_response_list <- list()

row_index <- 1L


for (variable in categorical_variables) {

  x <- clinical[[variable]]

  x_display <- as.character(x)

  x_display[
    is.na(x_display) |
      trimws(x_display) == ""
  ] <- "Missing"

  tab <- as.data.frame(
    table(
      response = clinical$response,
      level = x_display,
      useNA = "no"
    ),
    stringsAsFactors = FALSE
  )

  names(tab)[3] <- "n"

  response_totals <- table(
    clinical$response
  )

  tab$percent_within_response <- (
    100 *
      tab$n /
      as.numeric(
        response_totals[tab$response]
      )
  )

  tab$variable <- variable

  tab <- tab[
    ,
    c(
      "variable",
      "response",
      "level",
      "n",
      "percent_within_response"
    )
  ]

  clinical_by_response_list[[row_index]] <- tab

  row_index <- row_index + 1L
}


clinical_by_response <- do.call(
  rbind,
  clinical_by_response_list
)


write.csv(
  clinical_by_response,
  file.path(
    "results",
    "day4_clinical_by_response.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Response distribution figure
# ------------------------------------------------------------

response_counts <- table(
  factor(
    clinical$response,
    levels = c(
      "NonResponder",
      "Responder"
    )
  )
)

png(
  file.path(
    "figures",
    "day4_response_distribution.png"
  ),
  width = 1400,
  height = 1000,
  res = 150
)

barplot(
  response_counts,
  main = "PD1 cohort response distribution",
  ylab = "Number of samples",
  xlab = "Clinical response",
  ylim = c(
    0,
    max(response_counts) * 1.2
  )
)

text(
  x = seq_along(response_counts),
  y = response_counts,
  labels = response_counts,
  pos = 3
)

dev.off()


# ------------------------------------------------------------
# Library size figure
# ------------------------------------------------------------

png(
  file.path(
    "figures",
    "day4_library_size.png"
  ),
  width = 1800,
  height = 1000,
  res = 150
)

barplot(
  sample_qc$library_size,
  names.arg = sample_qc$sample_name,
  las = 2,
  cex.names = 0.7,
  main = "Library size by PD1 sample",
  ylab = "Total raw counts",
  xlab = "Sample"
)

dev.off()


# ------------------------------------------------------------
# Genes detected figure
# ------------------------------------------------------------

png(
  file.path(
    "figures",
    "day4_genes_detected.png"
  ),
  width = 1800,
  height = 1000,
  res = 150
)

barplot(
  sample_qc$genes_detected,
  names.arg = sample_qc$sample_name,
  las = 2,
  cex.names = 0.7,
  main = "Genes detected by PD1 sample",
  ylab = "Genes with count > 0",
  xlab = "Sample"
)

abline(
  h = gene_limits["lower"],
  lty = 2
)

abline(
  h = gene_limits["upper"],
  lty = 2
)

dev.off()


# ------------------------------------------------------------
# Percent zeros figure
# ------------------------------------------------------------

png(
  file.path(
    "figures",
    "day4_percent_zeros.png"
  ),
  width = 1800,
  height = 1000,
  res = 150
)

barplot(
  sample_qc$percent_zeros,
  names.arg = sample_qc$sample_name,
  las = 2,
  cex.names = 0.7,
  main = "Zero-count percentage by PD1 sample",
  ylab = "Zero counts (%)",
  xlab = "Sample"
)

abline(
  h = zero_limits["upper"],
  lty = 2
)

dev.off()


# ------------------------------------------------------------
# Detect case-only category inconsistencies
# ------------------------------------------------------------

find_case_variants <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]

  normalized <- tolower(
    trimws(x)
  )

  split_values <- split(
    x,
    normalized
  )

  variants <- lapply(
    split_values,
    unique
  )

  variants[
    vapply(
      variants,
      length,
      integer(1)
    ) > 1L
  ]
}


biopsy_type_variants <- find_case_variants(
  clinical[["Biopsy type"]]
)

biopsy_site_variants <- find_case_variants(
  clinical[["Biopsy site"]]
)


# ------------------------------------------------------------
# QC anomaly documentation
# ------------------------------------------------------------

outliers <- sample_qc[
  sample_qc$potential_outlier,
  ,
  drop = FALSE
]


doc_file <- file.path(
  "docs",
  "day4_qc_anomalies.md"
)


con <- file(
  doc_file,
  open = "wt"
)


writeLines(
  c(
    "# Day 4 — QC anomalies and observations",
    "",
    "This document records descriptive QC findings for the primary PD1 cohort.",
    "",
    "No sample was removed, no clinical value was imputed, and no predictive inference was performed.",
    "",
    "## Expression matrix integrity",
    "",
    paste0(
      "- Matrix dimensions: ",
      nrow(counts),
      " genes × ",
      ncol(counts),
      " PD1 samples."
    ),
    "- Missing count values: 0.",
    "- Negative count values: 0.",
    "- Non-integer-like count values: 0.",
    "",
    "## Potential sample-level QC outliers",
    ""
  ),
  con
)


if (nrow(outliers) == 0L) {

  writeLines(
    "No samples were flagged by the descriptive 1.5 × IQR rule.",
    con
  )

} else {

  for (i in seq_len(nrow(outliers))) {

    writeLines(
      c(
        paste0(
          "### ",
          outliers$sample_name[i]
        ),
        "",
        paste0(
          "- Library size: ",
          outliers$library_size[i],
          "."
        ),
        paste0(
          "- Genes detected: ",
          outliers$genes_detected[i],
          "."
        ),
        paste0(
          "- Zero-count percentage: ",
          round(
            outliers$percent_zeros[i],
            3
          ),
          "%."
        ),
        paste0(
          "- Library-size flag: ",
          outliers$flag_library_size[i],
          "."
        ),
        paste0(
          "- Genes-detected flag: ",
          outliers$flag_genes_detected[i],
          "."
        ),
        paste0(
          "- Percent-zero flag: ",
          outliers$flag_percent_zeros[i],
          "."
        ),
        "",
        "This is a descriptive QC flag only. The sample is retained pending further assessment.",
        ""
      ),
      con
    )
  }
}


writeLines(
  c(
    "Genes detected and zero-count percentage are mathematically related because both use the same count > 0 detection criterion. They should therefore not be interpreted as independent QC signals.",
    "",
    "## Clinical missingness",
    ""
  ),
  con
)


for (i in seq_len(nrow(clinical_missingness))) {

  writeLines(
    paste0(
      "- ",
      clinical_missingness$variable[i],
      ": ",
      clinical_missingness$missing_n[i],
      "/",
      clinical_missingness$total_n[i],
      " missing (",
      round(
        clinical_missingness$missing_percent[i],
        1
      ),
      "%)."
    ),
    con
  )
}


writeLines(
  c(
    "",
    "Missing BRAF and LDH values are retained as source-data missingness and are not treated as pipeline failures.",
    "",
    "## Clinical category consistency",
    ""
  ),
  con
)


write_variant_section <- function(
  title,
  variants,
  con
) {

  writeLines(
    paste0("### ", title),
    con
  )

  writeLines(
    "",
    con
  )

  if (length(variants) == 0L) {

    writeLines(
      "No case-only category variants detected.",
      con
    )

  } else {

    for (nm in names(variants)) {

      writeLines(
        paste0(
          "- ",
          paste(
            variants[[nm]],
            collapse = " / "
          )
        ),
        con
      )
    }
  }

  writeLines(
    "",
    con
  )
}


write_variant_section(
  "Biopsy type",
  biopsy_type_variants,
  con
)

write_variant_section(
  "Biopsy site",
  biopsy_site_variants,
  con
)


writeLines(
  c(
    "Raw clinical categories are preserved at this stage.",
    "",
    "Capitalization, punctuation or abbreviation differences are documented rather than silently harmonized.",
    "",
    "Any future harmonization must be explicit, reproducible and traceable.",
    "",
    "## Interpretation boundary",
    "",
    "All Day 4 analyses are descriptive quality-control analyses.",
    "",
    "No association testing, predictive modelling or causal interpretation is performed at this stage."
  ),
  con
)


close(con)


# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DAY 4 INITIAL QC GENERATED\n")
cat("============================================\n")

cat(
  "Samples:",
  ncol(counts),
  "\n"
)

cat(
  "Genes:",
  nrow(counts),
  "\n"
)

cat(
  "Potential QC outliers:",
  nrow(outliers),
  "\n"
)

if (nrow(outliers) > 0L) {
  cat(
    "Flagged samples:",
    paste(
      outliers$sample_name,
      collapse = ", "
    ),
    "\n"
  )
}

cat(
  "BRAF missing:",
  missing_n(clinical[["BRAF mutation"]]),
  "/ 36\n"
)

cat(
  "LDH missing:",
  missing_n(clinical[["LDH high"]]),
  "/ 36\n"
)

cat("\nOutputs written to results/, figures/ and docs/.\n")