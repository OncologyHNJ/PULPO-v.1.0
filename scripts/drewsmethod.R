#!/usr/bin/env Rscript
# ==================================================
# Script: drewsmethod.R
# Description: Cohort-level CNV signature extraction 
#              using the Drews et al. (2022) methodology.
# Author: Marta Portasany-Rodríguez
# Created on: 2025-09-11
# Last modified: 2025-11-25
# Pipeline: PULPO 🐙 
# Dependencies: CINSignatureQuantification, limSolve, dplyr, data.table
# ==================================================

# ── 1️⃣ Parse command-line arguments ───────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
base_dir   <- args[1]
output_dir <- args[2]

cat("🔹 Starting Drews CNV signature extraction...\n")

# ── 2️⃣ Load required libraries ─────────────────────────────────────
suppressPackageStartupMessages({
  library(CINSignatureQuantification)
  library(limSolve)
  library(dplyr)
  library(data.table)
})

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── 3️⃣ Collect all patient CNV files (with early filtering) ────────
patients <- list.files(base_dir)

all_data <- lapply(patients, function(p) {
  file_path <- file.path(base_dir, p, "drews", "drewsformatted.tsv")
  if (!file.exists(file_path)) {
    message("⚠️ File not found: ", file_path)
    return(NULL)
  }
  
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 🔴 FILTRO TEMPRANO: mantener solo muestras con ≥ 5 segmentos
  if (nrow(df) < 5) {
    message("⚠️ Sample ", p, " skipped: only ", nrow(df), " CNV segments (< 5).")
    return(NULL)
  }
  
  df
})

# Remove NULL entries
all_data <- all_data[!sapply(all_data, is.null)]

if (length(all_data) == 0) {
  cat("❌ No samples with ≥ 5 CNV segments. Drews analysis will not be performed.\n")
  
  # Crear outputs mínimos para que Snakemake no pete
  saveRDS(list(error = "No samples with ≥ 5 CNV segments"),
          file.path(output_dir, "drews_signatures_results.rds"))
  write.table(data.frame(),
              file = file.path(output_dir, "predictions.csv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  png(file.path(output_dir, "activities.png"), width = 1200, height = 800)
  plot.new(); text(0.5, 0.5, "Drews CNV signatures not computed\n(no samples with ≥ 5 segments)")
  dev.off()
  
  q(save = "no", status = 0)
}

# ── 4️⃣ Merge patient-level data ────────────────────────────────────
merged_data <- do.call(rbind, all_data)
merged_data <- merged_data[, c(setdiff(names(merged_data), "sample"), "sample")]
rownames(merged_data) <- NULL

# ── 5️⃣ Chromosome normalization ────────────────────────────────────
merged_data$chromosome <- as.character(merged_data$chromosome)
merged_data$chromosome[merged_data$chromosome == "23"] <- "X"
merged_data$chromosome[merged_data$chromosome == "24"] <- "Y"
merged_data$chromosome <- factor(merged_data$chromosome,
                                 levels = c(as.character(1:22), "X", "Y"))

# ── 6️⃣ Sorting ─────────────────────────────────────────────────────
merged_data$sample <- as.character(merged_data$sample)
setDT(merged_data)
setorder(merged_data, sample, chromosome, start)

# (Opcional: sanity check por si acaso)
counts <- merged_data[, .N, by = sample]
cat("✅ Merged CNV data from", nrow(counts), "samples.\n")
print(counts)

# ── 7️⃣ Save merged dataset ─────────────────────────────────────────
write.csv(merged_data,
          file = file.path(output_dir, "drews_all_patients.csv"),
          quote = FALSE, row.names = FALSE)
saveRDS(merged_data, file.path(output_dir, "drews_all_patients.rds"))

# ── 8️⃣ Run CINSignatureQuantification (Drews method) ───────────────
mySigs <- quantifyCNSignatures(
  merged_data,
  method = "drews",
  build = "hg38",
  experimentName = "subset"
)

saveRDS(mySigs, file.path(output_dir, "drews_signatures_results.rds"))

# ── 9️⃣ Predictions (Platinum classifier) ───────────────────────────
pred_raw <- clinPredictionPlatinum(object = mySigs)

# ── Helper function for normalization ───────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

normalize_predictions <- function(x) {
  if (is.vector(x) || is.factor(x)) {
    data.frame(
      sample_id = names(x) %||% paste0("sample_", seq_along(x)),
      predicted_label = as.character(x),
      stringsAsFactors = FALSE
    )
  } else if (is.matrix(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!is.null(rownames(df)) && !("sample_id" %in% names(df))) {
      df$sample_id <- rownames(df)
    }
    df
  } else if (is.data.frame(x)) {
    x
  } else {
    stop("Unrecognized structure for 'pred_raw': ", class(x)[1])
  }
}

pred_df <- normalize_predictions(pred_raw)

# Rename columns and harmonize names
prob_cols <- grep("prob|probab|score|resist|sens", names(pred_df),
                  ignore.case = TRUE, value = TRUE)
if (!"predicted_label" %in% names(pred_df)) {
  label_col <- intersect(tolower(names(pred_df)), c("prediction", "pred", "label", "class"))
  if (length(label_col) == 1) {
    names(pred_df)[tolower(names(pred_df)) == label_col] <- "predicted_label"
  }
}
if (!"sample_id" %in% names(pred_df)) {
  sid <- rownames(pred_df)
  if (is.null(sid) || all(sid == "")) {
    maybe_cols <- tryCatch(colnames(mySigs$exposures), error = function(e) NULL)
    if (is.null(maybe_cols)) maybe_cols <- tryCatch(colnames(mySigs), error = function(e) NULL)
    sid <- maybe_cols %||% paste0("sample_", seq_len(nrow(pred_df)))
  }
  pred_df$sample_id <- sid
}

keep <- c("sample_id", "predicted_label")
if (length(prob_cols) > 0) keep <- c(keep, prob_cols)
keep <- unique(intersect(keep, names(pred_df)))
pred_out <- pred_df[, keep, drop = FALSE]

# ── 🔍 QC summary ───────────────────────────────────────────────────
cat("\n📊 Prediction summary:\n")
print(table(pred_out$predicted_label, useNA = "ifany"))
if (any(grepl("prob", names(pred_out), ignore.case = TRUE))) {
  suppressWarnings(print(summary(pred_out[, grepl("prob", names(pred_out), ignore.case = TRUE), drop = FALSE])))
}

# ── 💾 Save prediction outputs ──────────────────────────────────────
saveRDS(pred_out, file.path(output_dir, "predictions.rds"))
write.table(pred_out, file = file.path(output_dir, "predictions.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 🔬 NA feature quality control (optional) ────────────────────────
if (!is.null(tryCatch(mySigs$features, error = function(e) NULL))) {
  feat <- mySigs$features
  na_by_sample <- colSums(is.na(feat))
  write.table(
    data.frame(sample_id = names(na_by_sample), n_na_features = as.integer(na_by_sample)),
    file = file.path(output_dir, "predictions_QC_NA.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# ── 🔹 Sample-by-Component visualization ────────────────────────────
samplebycomponent <- plotSampleByComponent(object = mySigs)
saveRDS(samplebycomponent, file.path(output_dir, "samplebycomponent.rds"))

png(file.path(output_dir, "samplebycomponent.png"), width = 1200, height = 800)
plotSampleByComponent(object = mySigs)
dev.off()

# ── 🔹 Activities visualization ─────────────────────────────────────
activities <- plotActivities(object = mySigs, type = "threshold")
write.csv(getActivities(mySigs), file = file.path(output_dir, "activitiesmatrix.csv"))

saveRDS(activities, file.path(output_dir, "activities.rds"))
png(file.path(output_dir, "activities.png"), width = 1200, height = 800)
plotActivities(object = mySigs, type = "threshold")
dev.off()

# ──  🔹 Save experiment metadata ───────────────────────────────────
experiment_info <- getExperiment(mySigs)
saveRDS(experiment_info, file.path(output_dir, "experiment_info.rds"))

cat("\n✅ Drews CNV signature analysis completed successfully!\n")

