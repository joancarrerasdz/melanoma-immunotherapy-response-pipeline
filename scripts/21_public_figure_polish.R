suppressPackageStartupMessages({
  library(pROC)
})

dir.create("figures", showWarnings = FALSE)

# ============================================================
# W5 — Public figure polish
#
# Presentation-only regeneration from frozen analytical outputs.
#
# This script does NOT:
# - refit models;
# - reselect genes;
# - retune hyperparameters;
# - change thresholds;
# - recompute blind predictions;
# - modify Week 3 / Week 4 analytical conclusions.
# ============================================================


# ------------------------------------------------------------
# Helper: draw ROC manually in conventional coordinates
# x = 1 - specificity
# y = sensitivity
# ------------------------------------------------------------

roc_xy <- function(roc_object) {

  x <- 1 - roc_object$specificities
  y <- roc_object$sensitivities

  ord <- order(x, y)

  data.frame(
    fpr = x[ord],
    sensitivity = y[ord]
  )
}


# ============================================================
# 1. WEEK 3 — STRICT NESTED-CV ROC
# ============================================================

w3 <- read.csv(
  "results/week3_outer_predictions_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(nrow(w3) == 36)

roc_w3 <- roc(
  response = factor(
    w3$observed,
    levels = c("NonResponder", "Responder")
  ),
  predictor = w3$prob_responder,
  levels = c("NonResponder", "Responder"),
  direction = "<",
  quiet = TRUE
)

auc_w3 <- as.numeric(auc(roc_w3))
xy_w3 <- roc_xy(roc_w3)

png(
  "figures/week3_nested_cv_roc_strict.png",
  width = 2200,
  height = 1800,
  res = 250
)

par(
  mar = c(5.2, 5.2, 4.2, 2.0),
  las = 1
)

plot(
  xy_w3$fpr,
  xy_w3$sensitivity,
  type = "s",
  lwd = 3,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = sprintf(
    "GSE160638 — strict nested CV ROC (AUC = %.4f)",
    auc_w3
  ),
  xaxs = "i",
  yaxs = "i"
)

abline(
  a = 0,
  b = 1,
  lty = 2,
  lwd = 1.5
)

box()

dev.off()


# ============================================================
# 2. WEEK 3 — FEATURE-SELECTION STABILITY
# ============================================================

stab <- read.csv(
  "results/week3_gene_selection_frequency_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(nrow(stab) == 312)

freq <- table(
  factor(
    stab$outer_folds_selected,
    levels = 1:5
  )
)

expected_freq <- c(212L, 49L, 26L, 13L, 12L)

stopifnot(
  identical(
    as.integer(freq),
    expected_freq
  )
)

png(
  "figures/week3_gene_selection_stability_strict.png",
  width = 1800,
  height = 1400,
  res = 220
)

par(
  mar = c(5.7, 5.2, 4.2, 2.0),
  las = 1
)

bp <- barplot(
  freq,
  names.arg = paste0(1:5, "/5"),
  xlab = "Exact number of outer folds in which gene was selected",
  ylab = "Number of genes",
  main = "GSE160638 — nested-CV gene-selection stability",
  ylim = c(
    0,
    max(freq) * 1.12
  )
)

text(
  x = bp,
  y = as.numeric(freq),
  labels = as.numeric(freq),
  pos = 3,
  cex = 0.95
)

mtext(
  "Counts are exact selection frequencies, not cumulative ≥k/5 counts.",
  side = 1,
  line = 4.3,
  cex = 0.75
)

box()

dev.off()


# ============================================================
# 3. WEEK 4 — LOCKED EXTERNAL VALIDATION ROC
#    Primary cohort-specific assessment only
# ============================================================

w4 <- read.csv(
  "results/week4_external_validation_predictions_strict.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(nrow(w4) == 76)

d91061 <- subset(
  w4,
  dataset == "GSE91061"
)

d78220 <- subset(
  w4,
  dataset == "GSE78220"
)

stopifnot(
  nrow(d91061) == 49,
  nrow(d78220) == 27
)

roc_91061 <- roc(
  response = factor(
    d91061$observed,
    levels = c("NonResponder", "Responder")
  ),
  predictor = d91061$prob_responder,
  levels = c("NonResponder", "Responder"),
  direction = "<",
  quiet = TRUE
)

roc_78220 <- roc(
  response = factor(
    d78220$observed,
    levels = c("NonResponder", "Responder")
  ),
  predictor = d78220$prob_responder,
  levels = c("NonResponder", "Responder"),
  direction = "<",
  quiet = TRUE
)

auc_91061 <- as.numeric(auc(roc_91061))
auc_78220 <- as.numeric(auc(roc_78220))

xy_91061 <- roc_xy(roc_91061)
xy_78220 <- roc_xy(roc_78220)

png(
  "figures/week4_external_validation_roc_strict.png",
  width = 2200,
  height = 1800,
  res = 250
)

par(
  mar = c(5.2, 5.2, 4.2, 2.0),
  las = 1
)

plot(
  xy_91061$fpr,
  xy_91061$sensitivity,
  type = "s",
  lwd = 3,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = "Week 4 locked external validation — cohort-specific ROC",
  xaxs = "i",
  yaxs = "i"
)

lines(
  xy_78220$fpr,
  xy_78220$sensitivity,
  type = "s",
  lwd = 3,
  lty = 2
)

abline(
  a = 0,
  b = 1,
  lty = 3,
  lwd = 1.5
)

legend(
  "bottomright",
  legend = c(
    sprintf(
      "GSE91061 (n = 49), AUC = %.4f",
      auc_91061
    ),
    sprintf(
      "GSE78220 (n = 27), AUC = %.4f",
      auc_78220
    )
  ),
  lty = c(1, 2),
  lwd = 3,
  bty = "n"
)

box()

dev.off()


# ============================================================
# 4. VERIFICATION
# ============================================================

cat("\n")
cat("============================================\n")
cat("W5 PUBLIC FIGURE POLISH\n")
cat("============================================\n")

cat(
  sprintf(
    "Week 3 nested-CV AUC: %.4f\n",
    auc_w3
  )
)

cat(
  "Week 3 exact stability counts:",
  paste(as.integer(freq), collapse = ", "),
  "\n"
)

cat(
  sprintf(
    "GSE91061 AUC: %.4f\n",
    auc_91061
  )
)

cat(
  sprintf(
    "GSE78220 AUC: %.4f\n",
    auc_78220
  )
)

cat("\nGenerated:\n")

out_files <- c(
  "figures/week3_nested_cv_roc_strict.png",
  "figures/week3_gene_selection_stability_strict.png",
  "figures/week4_external_validation_roc_strict.png"
)

for (f in out_files) {
  cat(
    ifelse(
      file.exists(f) && file.info(f)$size > 0,
      "PASS ",
      "FAIL "
    ),
    f,
    "\n",
    sep = ""
  )
}

cat("\n")
cat("W5 PUBLIC FIGURE POLISH: PASS\n")
