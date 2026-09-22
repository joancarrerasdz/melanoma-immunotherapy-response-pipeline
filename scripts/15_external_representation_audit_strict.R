#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
})

# ============================================================
# Week 4 — external representation audit
#
# PURPOSE
# -------
# Compare expression representations across:
#
#   - GSE160638 training cohort
#   - GSE91061 external cohort
#   - GSE78220 external cohort
#
# WITHOUT using external response labels.
#
# No classifier is fitted here.
# No threshold is selected here.
# No external performance is calculated here.
# ============================================================

results_dir <- "results"
processed_dir <- "data_processed"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

clean_gene_ids <- function(x) {

  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)
  x <- gsub("\\s+", "", x)

  toupper(x)
}

extract_core <- function(x, core_genes, dataset) {

  if (is.null(rownames(x))) {
    stop(
      dataset,
      ": expression matrix has no gene row names.",
      call. = FALSE
    )
  }

  gene_keys <- clean_gene_ids(
    rownames(x)
  )

  duplicated_core <- intersect(
    core_genes,
    unique(
      gene_keys[
        duplicated(gene_keys)
      ]
    )
  )

  if (length(duplicated_core) > 0) {
    stop(
      dataset,
      ": duplicated strict-core genes: ",
      paste(duplicated_core, collapse = ", "),
      call. = FALSE
    )
  }

  idx <- match(
    core_genes,
    gene_keys
  )

  if (anyNA(idx)) {

    missing_genes <- core_genes[
      is.na(idx)
    ]

    stop(
      dataset,
      ": missing strict-core genes: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  out <- x[
    idx,
    ,
    drop = FALSE
  ]

  rownames(out) <- core_genes

  out
}

read_core_csv <- function(path, core_genes, dataset) {

  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x) < 2) {
    stop(
      dataset,
      ": invalid external-expression file.",
      call. = FALSE
    )
  }

  genes <- clean_gene_ids(
    x[[1]]
  )

  x[[1]] <- NULL

  mat <- as.matrix(x)
  mode(mat) <- "numeric"

  rownames(mat) <- genes

  extract_core(
    mat,
    core_genes,
    dataset
  )
}

samplewise_rank <- function(x) {

  n_genes <- nrow(x)

  out <- vapply(
    seq_len(ncol(x)),
    function(j) {

      r <- rank(
        x[, j],
        ties.method = "average"
      )

      (r - 0.5) / n_genes
    },
    numeric(n_genes)
  )

  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)

  out
}

samplewise_z <- function(x) {

  out <- vapply(
    seq_len(ncol(x)),
    function(j) {

      v <- x[, j]

      s <- sd(v)

      if (!is.finite(s) || s <= 0) {
        stop(
          "Sample with zero or invalid within-sample SD.",
          call. = FALSE
        )
      }

      (v - mean(v)) / s
    },
    numeric(nrow(x))
  )

  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)

  out
}

summarize_global <- function(
  x,
  dataset,
  representation
) {

  values <- as.numeric(x)

  if (
    length(values) == 0 ||
    anyNA(values) ||
    any(!is.finite(values))
  ) {
    stop(
      dataset,
      " / ",
      representation,
      ": invalid expression values.",
      call. = FALSE
    )
  }

  q <- quantile(
    values,
    probs = c(
      0,
      0.01,
      0.25,
      0.50,
      0.75,
      0.99,
      1
    ),
    names = FALSE
  )

  data.frame(
    dataset = dataset,
    representation = representation,
    genes = nrow(x),
    samples = ncol(x),
    minimum = q[1],
    p01 = q[2],
    p25 = q[3],
    median = q[4],
    p75 = q[5],
    p99 = q[6],
    maximum = q[7],
    mean = mean(values),
    sd = sd(values),
    stringsAsFactors = FALSE
  )
}

summarize_genes <- function(
  x,
  dataset,
  representation
) {

  data.frame(
    dataset = dataset,
    representation = representation,
    gene = rownames(x),
    mean = apply(x, 1, mean),
    sd = apply(x, 1, sd),
    p25 = apply(
      x,
      1,
      quantile,
      probs = 0.25
    ),
    median = apply(
      x,
      1,
      median
    ),
    p75 = apply(
      x,
      1,
      quantile,
      probs = 0.75
    ),
    stringsAsFactors = FALSE
  )
}

profile_similarity <- function(
  training_matrix,
  external_matrix,
  external_dataset,
  representation
) {

  train_profile <- apply(
    training_matrix,
    1,
    median
  )

  external_profile <- apply(
    external_matrix,
    1,
    median
  )

  stopifnot(
    identical(
      names(train_profile),
      names(external_profile)
    )
  )

  data.frame(
    external_dataset = external_dataset,
    representation = representation,
    pearson_gene_median_correlation = cor(
      train_profile,
      external_profile,
      method = "pearson"
    ),
    spearman_gene_median_correlation = cor(
      train_profile,
      external_profile,
      method = "spearman"
    ),
    mean_absolute_gene_median_difference = mean(
      abs(
        train_profile -
        external_profile
      )
    ),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 1. Frozen Week 3 core
# ============================================================

core_file <- file.path(
  results_dir,
  "week3_consensus_core_5of5_strict.csv"
)

core <- read.csv(
  core_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"gene" %in% colnames(core)) {
  stop(
    "Column 'gene' missing from strict Week 3 core file.",
    call. = FALSE
  )
}

core_genes <- unique(
  clean_gene_ids(core$gene)
)

if (length(core_genes) != 12L) {
  stop(
    "Expected 12 strict core genes; found ",
    length(core_genes),
    ".",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("WEEK 4 EXPRESSION REPRESENTATION AUDIT\n")
cat("============================================================\n")

cat(
  "Frozen strict core:",
  length(core_genes),
  "genes\n"
)

cat(
  paste(core_genes, collapse = ", "),
  "\n"
)

# ============================================================
# 2. Training cohort
# ============================================================

counts_object <- readRDS(
  file.path(
    processed_dir,
    "GSE160638_raw_counts_PD1_aligned.rds"
  )
)

if (inherits(counts_object, "DGEList")) {

  counts <- counts_object$counts

} else if (
  is.matrix(counts_object) ||
  is.data.frame(counts_object)
) {

  counts <- as.matrix(counts_object)

} else {

  stop(
    "Unsupported GSE160638 counts object.",
    call. = FALSE
  )
}

mode(counts) <- "numeric"

if (
  anyNA(counts) ||
  any(!is.finite(counts)) ||
  any(counts < 0)
) {
  stop(
    "Invalid values in GSE160638 raw counts.",
    call. = FALSE
  )
}

dge <- DGEList(
  counts = counts
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

training_log <- cpm(
  dge,
  log = TRUE,
  prior.count = 1
)

training_log <- extract_core(
  training_log,
  core_genes,
  "GSE160638"
)

cat(
  "\nTraining representation:",
  nrow(training_log),
  "genes x",
  ncol(training_log),
  "samples\n"
)

# ============================================================
# 3. External cohorts
# ============================================================

gse91061_raw <- read_core_csv(
  file.path(
    processed_dir,
    "week4_GSE91061_core_expression_untransformed_strict.csv"
  ),
  core_genes,
  "GSE91061"
)

gse78220_raw <- read_core_csv(
  file.path(
    processed_dir,
    "week4_GSE78220_core_expression_untransformed_strict.csv"
  ),
  core_genes,
  "GSE78220"
)

if (
  any(gse91061_raw < 0) ||
  any(gse78220_raw < 0)
) {
  stop(
    "Negative FPKM values detected.",
    call. = FALSE
  )
}

gse91061_log <- log2(
  gse91061_raw + 1
)

gse78220_log <- log2(
  gse78220_raw + 1
)

cat(
  "GSE91061 representation:",
  nrow(gse91061_log),
  "genes x",
  ncol(gse91061_log),
  "samples\n"
)

cat(
  "GSE78220 representation:",
  nrow(gse78220_log),
  "genes x",
  ncol(gse78220_log),
  "samples\n"
)

# ============================================================
# 4. Common label-free representations
# ============================================================

training_rank <- samplewise_rank(
  training_log
)

gse91061_rank <- samplewise_rank(
  gse91061_log
)

gse78220_rank <- samplewise_rank(
  gse78220_log
)

training_z <- samplewise_z(
  training_log
)

gse91061_z <- samplewise_z(
  gse91061_log
)

gse78220_z <- samplewise_z(
  gse78220_log
)

# ============================================================
# 5. Global audit
# ============================================================

global_audit <- do.call(
  rbind,
  list(

    summarize_global(
      training_log,
      "GSE160638",
      "native_log_scale"
    ),

    summarize_global(
      gse91061_log,
      "GSE91061",
      "native_log_scale"
    ),

    summarize_global(
      gse78220_log,
      "GSE78220",
      "native_log_scale"
    ),

    summarize_global(
      training_rank,
      "GSE160638",
      "samplewise_rank"
    ),

    summarize_global(
      gse91061_rank,
      "GSE91061",
      "samplewise_rank"
    ),

    summarize_global(
      gse78220_rank,
      "GSE78220",
      "samplewise_rank"
    ),

    summarize_global(
      training_z,
      "GSE160638",
      "samplewise_z"
    ),

    summarize_global(
      gse91061_z,
      "GSE91061",
      "samplewise_z"
    ),

    summarize_global(
      gse78220_z,
      "GSE78220",
      "samplewise_z"
    )
  )
)

write.csv(
  global_audit,
  file.path(
    results_dir,
    "week4_expression_representation_global_audit_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. Gene-level audit
# ============================================================

gene_audit <- do.call(
  rbind,
  list(

    summarize_genes(
      training_log,
      "GSE160638",
      "native_log_scale"
    ),

    summarize_genes(
      gse91061_log,
      "GSE91061",
      "native_log_scale"
    ),

    summarize_genes(
      gse78220_log,
      "GSE78220",
      "native_log_scale"
    ),

    summarize_genes(
      training_rank,
      "GSE160638",
      "samplewise_rank"
    ),

    summarize_genes(
      gse91061_rank,
      "GSE91061",
      "samplewise_rank"
    ),

    summarize_genes(
      gse78220_rank,
      "GSE78220",
      "samplewise_rank"
    ),

    summarize_genes(
      training_z,
      "GSE160638",
      "samplewise_z"
    ),

    summarize_genes(
      gse91061_z,
      "GSE91061",
      "samplewise_z"
    ),

    summarize_genes(
      gse78220_z,
      "GSE78220",
      "samplewise_z"
    )
  )
)

write.csv(
  gene_audit,
  file.path(
    results_dir,
    "week4_expression_representation_gene_audit_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. Cross-dataset gene-profile similarity
# ============================================================

similarity <- do.call(
  rbind,
  list(

    profile_similarity(
      training_log,
      gse91061_log,
      "GSE91061",
      "native_log_scale"
    ),

    profile_similarity(
      training_log,
      gse78220_log,
      "GSE78220",
      "native_log_scale"
    ),

    profile_similarity(
      training_rank,
      gse91061_rank,
      "GSE91061",
      "samplewise_rank"
    ),

    profile_similarity(
      training_rank,
      gse78220_rank,
      "GSE78220",
      "samplewise_rank"
    ),

    profile_similarity(
      training_z,
      gse91061_z,
      "GSE91061",
      "samplewise_z"
    ),

    profile_similarity(
      training_z,
      gse78220_z,
      "GSE78220",
      "samplewise_z"
    )
  )
)

write.csv(
  similarity,
  file.path(
    results_dir,
    "week4_expression_representation_similarity_strict.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. Final assertions
# ============================================================

stopifnot(
  identical(
    rownames(training_log),
    core_genes
  ),
  identical(
    rownames(gse91061_log),
    core_genes
  ),
  identical(
    rownames(gse78220_log),
    core_genes
  ),
  ncol(training_log) == 36L,
  ncol(gse91061_log) == 49L,
  ncol(gse78220_log) == 27L
)

cat("\n============================================================\n")
cat("WEEK 4 EXPRESSION REPRESENTATION AUDIT: PASS\n")
cat("============================================================\n")

cat(
  "Training samples:",
  ncol(training_log),
  "\n"
)

cat(
  "GSE91061 samples:",
  ncol(gse91061_log),
  "\n"
)

cat(
  "GSE78220 samples:",
  ncol(gse78220_log),
  "\n"
)

cat(
  "Core genes:",
  length(core_genes),
  "\n"
)

cat(
  "External response labels used: FALSE\n"
)

cat(
  "Classifier fitted: FALSE\n"
)

cat(
  "External performance calculated: FALSE\n"
)

cat("\nOutputs:\n")

cat(
  "- results/week4_expression_representation_global_audit_strict.csv\n"
)

cat(
  "- results/week4_expression_representation_gene_audit_strict.csv\n"
)

cat(
  "- results/week4_expression_representation_similarity_strict.csv\n"
)
