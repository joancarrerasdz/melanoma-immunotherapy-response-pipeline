#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ============================================================
# Week 4 — strict external-input contract
#
# Purpose:
#   1. Reconstruct external expression matrices from raw files.
#   2. Match ONLY the validated pre-treatment samples.
#   3. Project both external datasets onto the frozen Week 3
#      5/5 core-consensus gene set.
#   4. Perform this preparation WITHOUT using response labels.
#
# IMPORTANT:
#   No external response information is used for preprocessing,
#   feature definition, model fitting, hyperparameter tuning,
#   threshold selection or any other model-development step.
# ============================================================

results_dir <- "results"
processed_dir <- "data_processed"
raw_dir <- "data_raw"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na")] <- NA_character_
  x
}

clean_gene_ids <- function(x) {
  x <- clean_text(x)
  x <- sub("\\.\\d+$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

sample_key <- function(x, dataset) {

  x <- toupper(trimws(as.character(x)))

  if (dataset == "GSE78220") {
    x <- gsub("\\.BASELINE$", "", x)
    x <- gsub("_BASELINE$", "", x)
    x <- gsub("-BASELINE$", "", x)
    x <- gsub("BASELINE$", "", x)
  }

  gsub("[^A-Z0-9]", "", x)
}

write_expression_matrix <- function(x, path) {

  out <- data.frame(
    gene = rownames(x),
    x,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(
    out,
    path,
    row.names = FALSE
  )
}

scale_summary <- function(x, dataset) {

  v <- as.numeric(x)

  stopifnot(
    length(v) > 0,
    all(is.finite(v))
  )

  q <- quantile(
    v,
    probs = c(0, 0.01, 0.25, 0.50, 0.75, 0.99, 1),
    na.rm = TRUE,
    names = FALSE
  )

  data.frame(
    dataset = dataset,
    minimum = q[1],
    p01 = q[2],
    p25 = q[3],
    median = q[4],
    p75 = q[5],
    p99 = q[6],
    maximum = q[7],
    zero_fraction = mean(v == 0),
    stringsAsFactors = FALSE
  )
}

align_external_samples <- function(
  expr,
  metadata,
  dataset,
  core_genes
) {

  required_metadata <- c(
    "sample_id",
    "sample_name"
  )

  missing_columns <- setdiff(
    required_metadata,
    colnames(metadata)
  )

  if (length(missing_columns) > 0) {
    stop(
      dataset,
      ": required metadata columns missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  expr_keys <- sample_key(
    colnames(expr),
    dataset
  )

  metadata_keys <- sample_key(
    metadata$sample_name,
    dataset
  )

  if (anyDuplicated(expr_keys)) {
    stop(
      dataset,
      ": duplicated expression sample keys.",
      call. = FALSE
    )
  }

  if (anyDuplicated(metadata_keys)) {
    stop(
      dataset,
      ": duplicated metadata sample keys.",
      call. = FALSE
    )
  }

  sample_index <- match(
    metadata_keys,
    expr_keys
  )

  if (anyNA(sample_index)) {

    missing_samples <- metadata$sample_name[
      is.na(sample_index)
    ]

    stop(
      dataset,
      ": metadata samples missing from expression: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  expr_aligned <- expr[
    ,
    sample_index,
    drop = FALSE
  ]

  colnames(expr_aligned) <- metadata$sample_id

  if (!identical(
    sample_key(metadata$sample_name, dataset),
    expr_keys[sample_index]
  )) {
    stop(
      dataset,
      ": sample alignment failed.",
      call. = FALSE
    )
  }

  missing_core <- setdiff(
    core_genes,
    rownames(expr_aligned)
  )

  if (length(missing_core) > 0) {
    stop(
      dataset,
      ": missing strict core genes: ",
      paste(missing_core, collapse = ", "),
      call. = FALSE
    )
  }

  gene_index <- match(
    core_genes,
    rownames(expr_aligned)
  )

  expr_core <- expr_aligned[
    gene_index,
    ,
    drop = FALSE
  ]

  if (!identical(
    rownames(expr_core),
    core_genes
  )) {
    stop(
      dataset,
      ": frozen gene order was not preserved.",
      call. = FALSE
    )
  }

  if (anyNA(expr_core)) {
    stop(
      dataset,
      ": NA values found in strict core matrix.",
      call. = FALSE
    )
  }

  if (any(!is.finite(expr_core))) {
    stop(
      dataset,
      ": non-finite values found in strict core matrix.",
      call. = FALSE
    )
  }

  list(
    matrix = expr_core,
    matched_samples = ncol(expr_core),
    duplicated_expression_keys = sum(duplicated(expr_keys)),
    duplicated_metadata_keys = sum(duplicated(metadata_keys)),
    missing_sample_matches = sum(is.na(sample_index))
  )
}

# ============================================================
# 1. Frozen Week 3 core-consensus candidates
# ============================================================

core <- read.csv(
  file.path(
    results_dir,
    "week3_consensus_core_5of5_strict.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"gene" %in% colnames(core)) {
  stop(
    "Column 'gene' missing from Week 3 core-consensus file.",
    call. = FALSE
  )
}

core_genes <- unique(
  clean_gene_ids(core$gene)
)

if (length(core_genes) != 12L) {
  stop(
    "Expected exactly 12 frozen Week 3 core genes; found ",
    length(core_genes),
    ".",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("WEEK 4 STRICT EXTERNAL INPUT CONTRACT\n")
cat("============================================================\n")
cat("Frozen Week 3 core genes:", length(core_genes), "\n")
cat(
  paste(core_genes, collapse = ", "),
  "\n"
)

# ============================================================
# 2. GSE91061
# ============================================================

cat("\n============================================================\n")
cat("GSE91061\n")
cat("============================================================\n")

gse91061_path <- file.path(
  raw_dir,
  "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
)

con <- gzfile(
  gse91061_path,
  open = "rt"
)

gse91061_raw <- read.csv(
  con,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

close(con)

gene_ids_91061 <- clean_text(
  gse91061_raw[[1]]
)

gse91061_raw[[1]] <- NULL

expr91061 <- as.matrix(
  gse91061_raw
)

mode(expr91061) <- "numeric"

mapped_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(gene_ids_91061),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

gene_symbols_91061 <- unname(
  mapped_symbols[gene_ids_91061]
)

keep_gene_91061 <- (
  !is.na(gene_symbols_91061) &
  nzchar(gene_symbols_91061)
)

expr91061 <- expr91061[
  keep_gene_91061,
  ,
  drop = FALSE
]

gene_symbols_91061 <- clean_gene_ids(
  gene_symbols_91061[keep_gene_91061]
)

duplicated_core_91061 <- intersect(
  core_genes,
  unique(
    gene_symbols_91061[
      duplicated(gene_symbols_91061)
    ]
  )
)

if (length(duplicated_core_91061) > 0) {
  stop(
    "GSE91061: duplicated mappings affect strict core genes: ",
    paste(duplicated_core_91061, collapse = ", "),
    call. = FALSE
  )
}

rownames(expr91061) <- gene_symbols_91061

expr91061 <- expr91061[
  !duplicated(rownames(expr91061)),
  ,
  drop = FALSE
]

metadata91061 <- read.csv(
  file.path(
    processed_dir,
    "metadata_GSE91061_harmonized_kept.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

aligned91061 <- align_external_samples(
  expr = expr91061,
  metadata = metadata91061,
  dataset = "GSE91061",
  core_genes = core_genes
)

write_expression_matrix(
  aligned91061$matrix,
  file.path(
    processed_dir,
    "week4_GSE91061_core_expression_untransformed_strict.csv"
  )
)

cat(
  "Source expression dimensions:",
  nrow(expr91061), "genes x",
  ncol(expr91061), "samples\n"
)

cat(
  "Validated external samples:",
  nrow(metadata91061),
  "\n"
)

cat(
  "Matched samples:",
  aligned91061$matched_samples,
  "\n"
)

cat(
  "Strict core genes:",
  nrow(aligned91061$matrix),
  "/ 12\n"
)

cat(
  "Final core matrix:",
  nrow(aligned91061$matrix), "x",
  ncol(aligned91061$matrix),
  "\n"
)

# ============================================================
# 3. GSE78220
# ============================================================

cat("\n============================================================\n")
cat("GSE78220\n")
cat("============================================================\n")

gse78220_raw <- readxl::read_excel(
  file.path(
    raw_dir,
    "GSE78220_PatientFPKM.xlsx"
  ),
  sheet = 1
)

gse78220_raw <- as.data.frame(
  gse78220_raw,
  check.names = FALSE
)

gene_ids_78220 <- clean_gene_ids(
  gse78220_raw[[1]]
)

gse78220_raw[[1]] <- NULL

expr78220 <- as.matrix(
  gse78220_raw
)

mode(expr78220) <- "numeric"

rownames(expr78220) <- gene_ids_78220

duplicated_core_78220 <- intersect(
  core_genes,
  unique(
    rownames(expr78220)[
      duplicated(rownames(expr78220))
    ]
  )
)

if (length(duplicated_core_78220) > 0) {
  stop(
    "GSE78220: duplicated rows affect strict core genes: ",
    paste(duplicated_core_78220, collapse = ", "),
    call. = FALSE
  )
}

expr78220 <- expr78220[
  !duplicated(rownames(expr78220)),
  ,
  drop = FALSE
]

keep_pre_treatment <- !grepl(
  "ONTX",
  toupper(colnames(expr78220))
)

expr78220 <- expr78220[
  ,
  keep_pre_treatment,
  drop = FALSE
]

metadata78220 <- read.csv(
  file.path(
    processed_dir,
    "metadata_GSE78220_harmonized_kept.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

aligned78220 <- align_external_samples(
  expr = expr78220,
  metadata = metadata78220,
  dataset = "GSE78220",
  core_genes = core_genes
)

write_expression_matrix(
  aligned78220$matrix,
  file.path(
    processed_dir,
    "week4_GSE78220_core_expression_untransformed_strict.csv"
  )
)

cat(
  "Source expression dimensions after OnTx exclusion:",
  nrow(expr78220), "genes x",
  ncol(expr78220), "samples\n"
)

cat(
  "Validated external samples:",
  nrow(metadata78220),
  "\n"
)

cat(
  "Matched samples:",
  aligned78220$matched_samples,
  "\n"
)

cat(
  "Strict core genes:",
  nrow(aligned78220$matrix),
  "/ 12\n"
)

cat(
  "Final core matrix:",
  nrow(aligned78220$matrix), "x",
  ncol(aligned78220$matrix),
  "\n"
)

# ============================================================
# 4. Audit outputs
# ============================================================

audit <- rbind(
  data.frame(
    dataset = "GSE91061",
    validated_samples = nrow(metadata91061),
    matched_samples = aligned91061$matched_samples,
    required_core_genes = length(core_genes),
    available_core_genes = nrow(aligned91061$matrix),
    duplicated_expression_sample_keys =
      aligned91061$duplicated_expression_keys,
    duplicated_metadata_sample_keys =
      aligned91061$duplicated_metadata_keys,
    missing_sample_matches =
      aligned91061$missing_sample_matches,
    response_used_in_expression_preparation = FALSE
  ),
  data.frame(
    dataset = "GSE78220",
    validated_samples = nrow(metadata78220),
    matched_samples = aligned78220$matched_samples,
    required_core_genes = length(core_genes),
    available_core_genes = nrow(aligned78220$matrix),
    duplicated_expression_sample_keys =
      aligned78220$duplicated_expression_keys,
    duplicated_metadata_sample_keys =
      aligned78220$duplicated_metadata_keys,
    missing_sample_matches =
      aligned78220$missing_sample_matches,
    response_used_in_expression_preparation = FALSE
  )
)

write.csv(
  audit,
  file.path(
    results_dir,
    "week4_external_input_contract_strict.csv"
  ),
  row.names = FALSE
)

scale_audit <- rbind(
  scale_summary(
    aligned91061$matrix,
    "GSE91061"
  ),
  scale_summary(
    aligned78220$matrix,
    "GSE78220"
  )
)

write.csv(
  scale_audit,
  file.path(
    results_dir,
    "week4_external_expression_scale_audit_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 5. Final gate
# ============================================================

stopifnot(
  aligned91061$matched_samples == 49L,
  aligned78220$matched_samples == 27L,
  nrow(aligned91061$matrix) == 12L,
  nrow(aligned78220$matrix) == 12L,
  aligned91061$missing_sample_matches == 0L,
  aligned78220$missing_sample_matches == 0L
)

cat("\n============================================================\n")
cat("WEEK 4 EXTERNAL INPUT CONTRACT: PASS\n")
cat("============================================================\n")

cat("GSE91061: 12 genes x 49 samples\n")
cat("GSE78220: 12 genes x 27 samples\n")
cat("External response labels used in expression preparation: FALSE\n")

cat("\nOutputs:\n")
cat(
  "- data_processed/week4_GSE91061_core_expression_untransformed_strict.csv\n"
)
cat(
  "- data_processed/week4_GSE78220_core_expression_untransformed_strict.csv\n"
)
cat(
  "- results/week4_external_input_contract_strict.csv\n"
)
cat(
  "- results/week4_external_expression_scale_audit_strict.csv\n"
)
