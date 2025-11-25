#!/usr/bin/env Rscript
# ==================================================
# Script: taomethod.R
# Description: CNV signature extraction following
#              Tao et al. (2023) using sigminer.
# Author: Marta Portasany-Rodríguez
# Created on: 2025-09-11
# Last modified: 2025-11-25
# Pipeline: PULPO 🐙
# Dependencies: sigminer, parallel
# ==================================================

# ── 1️⃣ Parse command-line arguments ───────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
base_dir   <- args[1]
output_dir <- args[2]
n_runs     <- as.numeric(args[3])

cat("🔹 Starting CN signature extraction (Tao method)...\n")

# ── 2️⃣ Install and load required libraries ─────────────────────────
if (!requireNamespace("sigminer", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ShixiangWang/sigminer")
}

suppressPackageStartupMessages({
  library(sigminer)
  library(parallel)
})

set.seed(2025) # for reproducibility

# ── 3️⃣ Detect available cores for parallelization ─────────────────
total_cores <- detectCores()
safe_cores <- max(1, total_cores - 8)  # leave 8 free cores for system
cat("🧠 Using", safe_cores, "of", total_cores, "available cores.\n")

# ── 4️⃣ Prepare output directory ───────────────────────────────────
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# ── 5️⃣ Collect and merge CNV patient files ─────────────────────────
patients <- list.files(base_dir)

all_data <- lapply(patients, function(p) {
  file_path <- file.path(base_dir, p, "tao", "taoformatted.tsv")
  if (file.exists(file_path)) {
    df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Ensure sample column exists
    if (!"sample" %in% names(df)) {
      df$sample <- p
    }
    df
  } else {
    message("⚠️ File not found: ", file_path)
    NULL
  }
})

all_data <- all_data[!sapply(all_data, is.null)]

if (length(all_data) == 0) {
  stop("❌ No taoformatted.tsv files found in ", base_dir)
}

# Merge into a single data frame
merged_data <- do.call(rbind, all_data)

# Optional: drop sex chromosomes to avoid sex-related bias
#if ("chromosome" %in% colnames(merged_data)) {
#  merged_data <- merged_data[
#    !merged_data$chromosome %in% c("X", "Y", "chrX", "chrY"),
#  ]
#}

# Reorder columns so that 'sample' is last (as in Tao examples)
merged_data <- merged_data[, c(setdiff(names(merged_data), "sample"), "sample")]
rownames(merged_data) <- NULL

# Ensure segVal is integer absolute copy number
if ("segVal" %in% colnames(merged_data)) {
  merged_data$segVal <- as.integer(round(merged_data$segVal))
} else {
  stop("❌ Column 'segVal' not found in merged_data.")
}

# ── 6️⃣ Detect allele-specific CN (minor_cn) and configure flags ────
has_minor_cn <- "minor_cn" %in% colnames(merged_data)

if (has_minor_cn) {
  message("ℹ Detected 'minor_cn' column: enabling allele-specific Tao-176 mode (LOH-aware).")
  genome_measure  <- "wg"
  complement_flag <- TRUE
  add_loh_flag    <- TRUE
} else {
  message("ℹ No 'minor_cn' column detected: using total copy-number Tao-176 mode (no LOH).")
  genome_measure  <- "called"
  complement_flag <- FALSE
  add_loh_flag    <- FALSE
}

# Optionally cap very high copy numbers globally if needed:
# options(sigminer.copynumber.max = 10L)

# ── 7️⃣ Read CNV data into a CopyNumber object ──────────────────────
cn <- read_copynumber(
  merged_data,
  seg_cols       = c("chromosome", "start", "end", "segVal"),
  genome_build   = "hg38",
  genome_measure = genome_measure,   # conditional
  complement     = complement_flag,  # conditional
  add_loh        = add_loh_flag,     # conditional
  verbose        = TRUE
)

# ── 8️⃣ Tally Tao-176 features (method = "X") ────────────────────────
cn_tally_X <- sig_tally(
  cn,
  method  = "X",
  add_loh = add_loh_flag   # same flag as above
)

# ── 8a️⃣ Visualisation: Tao-176 catalogue across the cohort ─────────
# This shows the global distribution of Tao features (176 bins) across all samples
png(file.path(output_dir, "Tao_catalogue.png"),
    width = 2400, height = 800, res = 200)
print(
  show_catalogue(
    cn_tally_X,
    mode   = "copynumber",
    method = "X",
    style  = "cosmic",
    y_tr   = function(x) log10(x + 1),
    y_lab  = "log10(count + 1)"
  )
)
dev.off()

# ── 8b️⃣ Visualisation: Tao feature distributions (optional QC plot) ─
# This shows the distribution of feature categories (e.g. size, type, etc.)
#png(file.path(output_dir, "Tao_features_distribution.png"),
 #   width = 2000, height = 1200, res = 200)
#print(
 # show_cn_features(
 #   cn_tally_X$features,
 #   method = "X"
 # )
#)


# ── 9️⃣ Matrix cleaning and NA handling ─────────────────────────────
nmf_mat_raw <- cn_tally_X$nmf_matrix

anyNA(nmf_mat_raw)
mat_clean <- nmf_mat_raw[
  rowSums(nmf_mat_raw, na.rm = TRUE) > 0,
  colSums(nmf_mat_raw, na.rm = TRUE) > 0
]
mat_clean[is.na(mat_clean)] <- 0

cat("✅ Clean matrix dimensions:", dim(mat_clean), "\n")

# ── 🔟 Estimate optimal number of CNV signatures ────────────────────
estimate <- sig_estimate(
  mat_clean,
  nrun       = n_runs,
  save_plots = TRUE,   # this will output cophenetic/silhouette vs rank plots
  cores      = safe_cores
)

survey <- estimate$survey

# Use stability and separation quality criteria
n_sig_coph <- survey$rank[which.max(survey$cophenetic)]
n_sig_sil  <- survey$rank[which.max(survey$silhouette.consensus)]

# Combined heuristic: cophenetic > 0.95 & silhouette > 0.7
good_ranks <- survey$rank[
  survey$cophenetic > 0.95 & survey$silhouette.consensus > 0.7
]

if (length(good_ranks) == 0) {
  n_sig <- n_sig_coph
} else {
  n_sig <- max(good_ranks)
}

cat("🔹 Estimated number of CNV signatures:", n_sig, "\n")

# ── 1️⃣1️⃣ Extract CNV signatures (Sigminer NMF) ────────────────────
sig_res <- sig_extract(
  mat_clean,
  n_sig = n_sig,
  nrun  = n_runs,
  cores = safe_cores
)

# --- Normalized signature and exposure matrices (for downstream analyses) ---

# Signatures:
sig_mat_raw     <- sig_signature(sig_res, normalize = "raw")
sig_mat_feature <- sig_signature(sig_res, normalize = "row") 

# Exposures:
expo_abs <- sig_exposure(sig_res, type = "absolute")
expo_rel <- sig_exposure(sig_res, type = "relative")


write.csv(
  sig_mat_raw,
  file  = file.path(output_dir, "Tao_signatures_rawMatrix.csv"),
  quote = FALSE
)

write.csv(
  sig_mat_feature,
  file  = file.path(output_dir, "Tao_signatures_rowNorm.csv"),
  quote = FALSE
)

write.csv(
  expo_abs,
  file  = file.path(output_dir, "Tao_exposures_absolute.csv"),
  quote = FALSE
)

write.csv(
  expo_rel,
  file  = file.path(output_dir, "Tao_exposures_relative.csv"),
  quote = FALSE
)

# ── 1️⃣2️⃣ Compare de novo signatures with Tao reference catalogues ─
# Tao 176-feature catalogues:
#   - CNS_PCAWG176 (PCAWG-based)
#   - CNS_TCGA176  (TCGA-based, optional)

# PCAWG176 similarity
sim_PCAWG176 <- get_sig_similarity(
  sig_res,
  sig_db = "CNS_PCAWG176"
)

if (!is.null(sim_PCAWG176)) {
  write.table(
    sim_PCAWG176,
    file      = file.path(output_dir, "Tao_similarity_CNS_PCAWG176.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
}

# (Optional) TCGA176 similarity
#sim_TCGA176 <- get_sig_similarity(
#  sig_res,
#  sig_db = "CNS_TCGA176"
#)

#if (!is.null(sim_TCGA176)) {
#  write.table(
#    sim_TCGA176,
#    file      = file.path(output_dir, "Tao_similarity_CNS_TCGA176.tsv"),
#    sep       = "\t",
#    quote     = FALSE,
#    row.names = FALSE
#  )
#}

# ── 1️⃣3️⃣ Save and visualize results ───────────────────────────────
# Signatures profile plot (COSMIC-like layout)
png(file.path(output_dir, "Tao_signatures.png"),
    width = 2000, height = 1200, res = 200)
show_sig_profile(
  sig_res,
  mode      = "copynumber",
  normalize = "feature",
  method    = "T",       
  style     = "cosmic"
)
dev.off()

# Exposures plot (per-sample contribution)
png(file.path(output_dir, "Tao_exposures.png"),
    width = 2000, height = 1200, res = 200)
show_sig_exposure(sig_res)
dev.off()

# Export raw signature and exposure matrices
write.csv(
  sig_res$Signature,
  file  = file.path(output_dir, "Tao_signatures_raw.csv"),
  quote = FALSE
)

write.csv(
  sig_res$Exposure,
  file  = file.path(output_dir, "Tao_exposures_raw.csv"),
  quote = FALSE
)

cat("💾 CNV signature extraction (Tao method) completed successfully!\n")
cat("Results saved in:", output_dir, "\n")

