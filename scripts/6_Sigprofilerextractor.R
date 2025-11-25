#!/usr/bin/env Rscript
# ==================================================
# Script: 6_Sigprofilerextractor.R
# Description: Runs SigProfilerExtractorR for Structural Variants (SVs)
#              or Copy Number Variants (CNVs) at the individual level.
# Author: Marta Portasany
# Created: 2025-02-27
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: reticulate, devtools, SigProfilerExtractorR
# ==================================================

# ───────────────────────────────────────────────
# 0️⃣ LOAD ARGUMENTS
# ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
inputdata         <- args[1]
pythondirectory   <- args[2]
output            <- args[3]
minimum_signatures <- as.numeric(args[4])
maximum_signatures <- as.numeric(args[5])
nmf_replicates     <- as.numeric(args[6])

# ───────────────────────────────────────────────
# 1️⃣ LOAD LIBRARIES
# ───────────────────────────────────────────────
suppressPackageStartupMessages({
  library(reticulate)
  library(devtools)
})

# ───────────────────────────────────────────────
# 2️⃣ INSTALL / LOAD SIGPROFILEREXTRACTOR
# ───────────────────────────────────────────────
if (!requireNamespace("SigProfilerExtractorR", quietly = TRUE)) {
  cat("🔹 Installing SigProfilerExtractorR package from GitHub...\n")
  install.packages("devtools", repos = "http://cran.us.r-project.org")
  devtools::install_github("AlexandrovLab/SigProfilerExtractorR")
  library(SigProfilerExtractorR)
  SigProfilerExtractorR::install("GRCh38", rsync = FALSE, bash = TRUE)
} else {
  cat("✅ SigProfilerExtractorR already available.\n")
  library(SigProfilerExtractorR)
}

# ───────────────────────────────────────────────
# 3️⃣ PYTHON CONFIGURATION
# ───────────────────────────────────────────────
use_python(pythondirectory)
py_run_string("import sys")
py_config()

# ───────────────────────────────────────────────
# 4️⃣ CHECK INPUT FILE
# ───────────────────────────────────────────────
is_empty_or_header_only <- function(file) {
  lines <- readLines(file, warn = FALSE)
  return(length(lines) <= 1)
}

if (!file.exists(inputdata) || file.info(inputdata)$size == 0 || is_empty_or_header_only(inputdata)) {
  cat("⚠️ Input file is empty or header-only. Creating placeholder and exiting.\n")

  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  empty_output <- file.path(output, "empty_result.txt")
  file.create(empty_output)

  cat("✅ Placeholder created at:", empty_output, "\n")
  quit(status = 0)
}

# ───────────────────────────────────────────────
# 5️⃣ SEED HANDLING (REPRODUCIBILITY)
# ───────────────────────────────────────────────
seeds_path <- file.path(output, "Seeds.txt")

if (file.exists(seeds_path)) {
  cat("🔁 Reusing existing Seeds.txt for reproducibility:\n  ", seeds_path, "\n")
} else {
  cat("🧬 No Seeds.txt found — letting SigProfilerExtractor generate one automatically.\n")
}

# ───────────────────────────────────────────────
# 6️⃣ RUN SIGPROFILEREXTRACTOR
# ───────────────────────────────────────────────
cat("🔹 Running SigProfilerExtractor with parameters:\n",
    "   Input:", inputdata, "\n",
    "   Output:", output, "\n",
    "   Min signatures:", minimum_signatures, "\n",
    "   Max signatures:", maximum_signatures, "\n",
    "   NMF replicates:", nmf_replicates, "\n"
)

tryCatch({
  SigProfilerExtractorR::sigprofilerextractor(
    input_type = "matrix",
    output = output,
    input_data = inputdata,
    reference_genome = "GRCh38",
    opportunity_genome = "GRCh38",
    minimum_signatures = minimum_signatures,
    maximum_signatures = maximum_signatures,
    nmf_replicates = nmf_replicates
  )

  cat("✅ SigProfilerExtractor finished successfully!\n")

}, error = function(e) {
  cat("❌ ERROR during SigProfilerExtractor execution:\n", e$message, "\n")
  quit(status = 1)
})

cat("✔️ Completed SigProfilerExtractorR execution.\n")

