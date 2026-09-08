################################
# 01_load_data.R
################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

data_raw_dir <- file.path(base_dir, "data_raw")
data_processed_dir <- file.path(base_dir, "data_processed")

dir.create(data_processed_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================
# Official clinical supplement
# ===========================

clinical_dir <- file.path(data_raw_dir, "clinical")
dir.create(clinical_dir, showWarnings = FALSE, recursive = TRUE)

table_s2_url <- paste0(
  "https://ndownloader.figshare.com/files/",
  "59867554"
)

table_s2_path <- file.path(
  clinical_dir,
  "cir-22-0184_table_s2_suppst2.xlsx"
)

is_valid_xlsx <- function(path) {
  if (!file.exists(path)) return(FALSE)

  file_size <- file.info(path)$size
  if (is.na(file_size) || file_size < 1024) return(FALSE)

  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)

  signature <- readBin(connection, what = "raw", n = 2L)
  identical(signature, charToRaw("PK"))
}

download_table_s2 <- function(url, destination) {
  if (is_valid_xlsx(destination)) {
    message("Using cached Table S2: ", destination)
    return(invisible(destination))
  }

  if (file.exists(destination)) {
    unlink(destination)
  }

  temporary_file <- tempfile(
    pattern = "table_s2_",
    tmpdir = dirname(destination),
    fileext = ".xlsx"
  )

  on.exit({
    if (file.exists(temporary_file)) unlink(temporary_file)
  }, add = TRUE)

  previous_timeout <- getOption("timeout")
  on.exit(options(timeout = previous_timeout), add = TRUE)
  options(timeout = max(300L, previous_timeout))

  status <- utils::download.file(
    url = url,
    destfile = temporary_file,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )

  if (!identical(status, 0L) || !is_valid_xlsx(temporary_file)) {
    stop(
      "Table S2 download failed or did not produce a valid XLSX file.",
      call. = FALSE
    )
  }

  if (!file.rename(temporary_file, destination)) {
    stop("Could not move Table S2 to its final location.", call. = FALSE)
  }

  message(
    "Downloaded Table S2: ",
    destination,
    " (",
    file.info(destination)$size,
    " bytes)"
  )

  invisible(destination)
}

download_table_s2(table_s2_url, table_s2_path)

# ============================
# Read official clinical S2A
# ============================

available_sheets <- readxl::excel_sheets(table_s2_path)

if (!"S2A" %in% available_sheets) {
  stop("Required worksheet 'S2A' was not found in Table S2.", call. = FALSE)
}

clinical_s2a <- readxl::read_excel(
  path = table_s2_path,
  sheet = "S2A",
  skip = 2,
  na = c("", "NA"),
  .name_repair = "unique"
)

names(clinical_s2a)[1] <- "patient_id"

required_s2a_columns <- c("patient_id", "Cohort", "Response")
missing_s2a_columns <- setdiff(required_s2a_columns, names(clinical_s2a))

if (length(missing_s2a_columns) > 0L) {
  stop(
    "Missing required S2A columns: ",
    paste(missing_s2a_columns, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(clinical_s2a) != 73L || ncol(clinical_s2a) != 21L) {
  stop(
    "Unexpected S2A dimensions: ",
    nrow(clinical_s2a),
    " rows x ",
    ncol(clinical_s2a),
    " columns.",
    call. = FALSE
  )
}

message(
  "Read Table S2A: ",
  nrow(clinical_s2a),
  " rows x ",
  ncol(clinical_s2a),
  " columns"
)

# =============================
# Normalize patient identifiers
# =============================

clinical_s2a <- clinical_s2a %>%
  mutate(
    patient_id_raw = as.character(patient_id),
    patient_id = toupper(trimws(patient_id_raw)),
    patient_id = gsub("[^A-Z0-9]+", "_", patient_id),
    patient_id = gsub("_+", "_", patient_id),
    patient_id = gsub("^_+|_+$", "", patient_id),
    .before = patient_id
  )

pd1_rows <- !is.na(clinical_s2a$Cohort) &
  clinical_s2a$Cohort == "PD1"

pd1_ids <- clinical_s2a$patient_id[pd1_rows]

if (
  anyNA(pd1_ids) ||
  any(!grepl("^PD1_[0-9]+$", pd1_ids))
) {
  stop("Invalid normalized PD1 patient identifier detected.", call. = FALSE)
}

if (anyDuplicated(pd1_ids) > 0L) {
  stop("Duplicated PD1 patient identifier detected.", call. = FALSE)
}

message(
  "Normalized PD1 identifiers: ",
  length(pd1_ids),
  " unique IDs"
)

# ===================
# Filter PD1 cohort
# ===================

clinical_pd1 <- clinical_s2a %>%
  filter(Cohort == "PD1")

if (nrow(clinical_pd1) != 36L) {
  stop(
    "Unexpected PD1 cohort size: ",
    nrow(clinical_pd1),
    " patients; expected 36.",
    call. = FALSE
  )
}

if (!all(clinical_pd1$Cohort == "PD1")) {
  stop("Non-PD1 records remained after cohort filtering.", call. = FALSE)
}

message(
  "Filtered PD1 cohort: ",
  nrow(clinical_pd1),
  " patients"
)

# ========================
# Define clinical endpoint
# ========================

clinical_pd1 <- clinical_pd1 %>%
  mutate(
    response_recist = toupper(trimws(as.character(Response))),
    response = case_when(
      response_recist %in% c("CR", "PR") ~ "Responder",
      response_recist %in% c("SD", "PD") ~ "NonResponder",
      TRUE ~ NA_character_
    )
  )

if (anyNA(clinical_pd1$response)) {
  stop("Unmapped or missing RECIST response detected.", call. = FALSE)
}

recist_counts <- table(
  factor(
    clinical_pd1$response_recist,
    levels = c("CR", "PR", "SD", "PD")
  )
)

expected_recist_counts <- c(
  CR = 6L,
  PR = 16L,
  SD = 2L,
  PD = 12L
)

if (!identical(as.integer(recist_counts), unname(expected_recist_counts))) {
  stop("Unexpected RECIST response counts.", call. = FALSE)
}

class_counts <- table(
  factor(
    clinical_pd1$response,
    levels = c("Responder", "NonResponder")
  )
)

expected_class_counts <- c(
  Responder = 22L,
  NonResponder = 14L
)

if (!identical(as.integer(class_counts), unname(expected_class_counts))) {
  stop("Unexpected binary endpoint counts.", call. = FALSE)
}

message(
  "RECIST counts: CR=6, PR=16, SD=2, PD=12"
)

message(
  "Endpoint counts: Responder=22, NonResponder=14"
)

# ===========================
# 1. Carregar counts
# ===========================

counts_path <- file.path(data_raw_dir, "GSE160638_combined_raw_counts.csv.gz")

counts <- read.csv(
  counts_path,
  row.names = 1,
  check.names = FALSE
)

colnames(counts) <- trimws(colnames(counts))
rownames(counts) <- trimws(rownames(counts))

cat("Dimensions counts:", dim(counts), "\n")
cat("Primeres files:\n")
print(head(rownames(counts)))
cat("Primeres mostres:\n")
print(head(colnames(counts)))

saveRDS(
  counts,
  file.path(data_processed_dir, "GSE160638_raw_counts.rds")
)

# ===========================
# 2. Descarregar metadata de GEO
# ===========================

gse <- getGEO("GSE160638", GSEMatrix = TRUE)
eset <- gse[[1]]

pheno <- pData(eset) |> as.data.frame()

write.csv(
  pheno,
  file.path(data_processed_dir, "GSE160638_pheno_FULL.csv"),
  row.names = FALSE
)

# ===========================
# 3. Construir metadata mínima
# ===========================

char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)

pheno$char_all <- apply(
  pheno[, char_cols, drop = FALSE],
  1,
  paste,
  collapse = " | "
)

# ==========================================
# Join GEO samples to the official endpoint
# ==========================================

sample_name_clean <- sub(
  ".*\\(([^()]*)\\)\\s*$",
  "\\1",
  trimws(pheno$title)
)

geo_meta <- transmute(
  pheno,
  sample_id = geo_accession,
  patient_id = sample_name_clean,
  title = trimws(title),
  source_name = source_name_ch1,
  characteristics = char_all
) %>%
  mutate(
    patient_id = toupper(trimws(as.character(patient_id))),
    patient_id = gsub("[^A-Z0-9]+", "_", patient_id),
    patient_id = gsub("_+", "_", patient_id),
    patient_id = gsub("^_+|_+$", "", patient_id)
  )

unmatched_pd1_ids <- setdiff(
  clinical_pd1$patient_id,
  geo_meta$patient_id
)

if (length(unmatched_pd1_ids) > 0L) {
  stop(
    "PD1 identifiers missing from GEO metadata: ",
    paste(unmatched_pd1_ids, collapse = ", "),
    call. = FALSE
  )
}

meta <- geo_meta %>%
  inner_join(
    clinical_pd1 %>%
      select(patient_id, response_recist, response),
    by = "patient_id"
  ) %>%
  mutate(
    sample_name = patient_id,
    .after = sample_id
  )

if (
  nrow(meta) != 36L ||
  anyDuplicated(meta$patient_id) > 0L ||
  anyNA(meta$response)
) {
  stop("Invalid GEO–clinical endpoint join.", call. = FALSE)
}

message(
  "Joined official endpoint to GEO metadata: ",
  nrow(meta),
  " PD1 samples"
)

write.csv(
  meta,
  file.path(data_processed_dir, "metadata_GSE160638.csv"),
  row.names = FALSE
)

cat("\nMetadata generada correctament.\n")
cat("Taula response:\n")
print(table(meta$response, useNA = "ifany"))

cat("\nPrimeres sample_name netes:\n")
print(head(meta$sample_name))

cat("\nMostres de counts no presents en metadata:\n")
print(setdiff(colnames(counts), meta$sample_name))

cat("\nMostres de metadata no presents en counts:\n")
print(setdiff(meta$sample_name, colnames(counts)))

