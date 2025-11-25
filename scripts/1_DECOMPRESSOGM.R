#!/usr/bin/env Rscript
# ==================================================
# Script: 1.0_DECOMPRESSOGM.R
# Description: Decompresses raw Optical Genome Mapping (OGM) ZIP files
#              into patient-specific directories for downstream analysis.
# Author: Marta Portasany
# Updated: 2025-10-28
# Pipeline: PULPO
# Dependencies: stringr, readr
# ==================================================

suppressPackageStartupMessages({
  library(stringr)
  library(readr)
})

# ───────────────────────────────────────────────
# 0️⃣ LOAD ARGUMENTS AND CONFIGURATION
# ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
ogm_directory <- args[1]
work_directory <- args[2]
samples_file  <- args[3]

cat("🔹 OGM Decompression started\n")
cat("🔹 Input directory:", ogm_directory, "\n")
cat("🔹 Working directory:", work_directory, "\n")
cat("🔹 Samples file:", samples_file, "\n")

# ───────────────────────────────────────────────
# 1️⃣ LOAD SAMPLE CONFIGURATION
# ───────────────────────────────────────────────
samples <- read_tsv(samples_file, col_types = cols())

# Validate required columns
if (!all(c("sample", "anonymised") %in% colnames(samples))) {
  stop("The configuration file must contain the columns 'sample' and 'anonymised'.")
}

# ───────────────────────────────────────────────
# 2️⃣ HELPER FUNCTION: CLEAN DIRECTORY
# ───────────────────────────────────────────────
filter_and_remove_files <- function(directory) {
  smap_pattern <- "\\.smap$"  # Keep .smap files
  cnv_pattern  <- "CNV"       # Keep CNV-related files

  files_in_directory <- list.files(directory, full.names = TRUE)
  cat("🧹 Cleaning directory:", directory, "\n")

  for (file_path in files_in_directory) {
    file_name <- basename(file_path)
    if (!(str_detect(file_name, smap_pattern) | str_detect(file_name, cnv_pattern))) {
      file.remove(file_path)
      cat("   Removed:", file_name, "\n")
    }
  }
}

# ───────────────────────────────────────────────
# 3️⃣ ITERATE OVER SAMPLES
# ───────────────────────────────────────────────
for (sample_id in samples$sample) {
  sample_dir <- file.path(ogm_directory, sample_id)
  cat("\n---------------------------------------------\n")
  cat("🧬 Processing sample:", sample_id, "\n")
  cat("📁 Directory:", sample_dir, "\n")

  if (!dir.exists(sample_dir)) {
    warning(paste("⚠️ Directory does not exist:", sample_dir))
    next
  }

  setwd(sample_dir)
  zip_files <- list.files(sample_dir, pattern = "\\.zip$", full.names = FALSE)

  # ── ZIP Filtering and Decompression ─────────────────────────────
  if (length(zip_files) > 0) {

    # Exclude secondary zips (BED_ALL, SV_CNV, BED, etc.)
    valid_zips <- grep("BED_ALL|SV_CNV|BED", zip_files, value = TRUE, invert = TRUE)

    if (length(valid_zips) == 0) {
      warning(paste("⚠️ No valid ZIP files found for:", sample_id, "(only BED/SV_CNV present)"))
      next
    }

    # If multiple valid zips, select the first one
    selected_zip <- valid_zips[1]
    cat("🗜️ Selected ZIP:", selected_zip, "\n")

    # Define output path
    anonymised_name <- samples$anonymised[samples$sample == sample_id]
    output_dir <- file.path(work_directory, "results/DATA/Patients", anonymised_name, "OGMdata")

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    unzip(selected_zip, exdir = output_dir)
    cat("📦 Decompressed:", selected_zip, "→", output_dir, "\n")

    # Clean directory (keep .smap + CNV)
    filter_and_remove_files(output_dir)

    # ── Handle Annotated vs Raw smaps ─────────────────────────────
    smaps <- list.files(output_dir, pattern = "\\.smap$", full.names = TRUE)
    annotated <- grep("Annotated", smaps, value = TRUE)

    if (length(annotated) > 0) {
      not_annotated <- setdiff(smaps, annotated)
      if (length(not_annotated) > 0) {
        file.remove(not_annotated)
        cat("🧾 Removed non-Annotated smap(s):", basename(not_annotated), "\n")
      }
    } else if (length(smaps) > 0) {
      cat("⚠️ No Annotated smap found for:", sample_id, "— keeping raw smap.\n")
    } else {
      warning(paste("⚠️ No smap files found at all for:", sample_id))
    }

    # ── Create .done sentinel ─────────────────────────────
    done_file <- file.path(output_dir, ".done")
    writeLines(paste(Sys.time(), "OK", anonymised_name), done_file)
    cat("✅ Created sentinel:", done_file, "\n")

  } else {
    cat("⚠️ No ZIP files found for:", sample_id, "\n")
  }
}

cat("\n✔️ OGM decompression completed for all available samples.\n")

