#!/usr/bin/env Rscript
# ==================================================
# Script: 8_Sigprofilerextractorcohort.R
# Description: Runs SigProfilerExtractor for cohort-level analysis
#              of CNV or SV mutational matrices generated across samples.
# Author: Marta Portasany
# Created on: 2025-10-31
# Pipeline: PULPO 🐙
# Dependencies: reticulate, devtools, SigProfilerExtractorR
# ==================================================

# ── Command-line arguments ───────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
inputdatacohort    <- args[1]
pythondirectory    <- args[2]
outputcohort       <- args[3]
minimum_signatures <- as.numeric(args[4])
maximum_signatures <- as.numeric(args[5])
nmf_replicates     <- as.numeric(args[6])

cat("\n🧬 Starting SigProfilerExtractor (Cohort Mode)\n")
cat("Input matrix:       ", inputdatacohort, "\n")
cat("Output directory:   ", outputcohort, "\n")
cat("Min signatures:     ", minimum_signatures, "\n")
cat("Max signatures:     ", maximum_signatures, "\n")
cat("NMF replicates:     ", nmf_replicates, "\n\n")

# ── Load required libraries ──────────────────────────────
suppressPackageStartupMessages({
  library("reticulate")
  library("devtools")
})

# ── Check and install SigProfilerExtractorR if needed ────
if (!requireNamespace("SigProfilerExtractorR", quietly = TRUE)) {
  cat("⚙️  Installing SigProfilerExtractorR from GitHub...\n")
  install.packages("devtools", repos = "http://cran.us.r-project.org")
  devtools::install_github("AlexandrovLab/SigProfilerExtractorR")
  library("SigProfilerExtractorR")
  SigProfilerExtractorR::install("GRCh38", rsync = FALSE, bash = TRUE)
} else {
  library("SigProfilerExtractorR")
  message("✅ SigProfilerExtractorR is already installed.")
}

# ── Configure Python environment ─────────────────────────
use_python(pythondirectory)
py_run_string("import sys")
py_config()

# ── Global reproducibility seed ──────────────────────────
global_seed_path <- file.path(dirname(dirname(outputcohort)), "Seeds_master.txt")

if (!file.exists(global_seed_path)) {
  global_seed <- as.integer(runif(1, min = 1, max = 2^31 - 1))
  write.table(data.frame(Seed = global_seed),
              file = global_seed_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("🧬 Created global seed file:", global_seed_path, "with seed", global_seed, "\n")
} else {
  global_seed <- as.numeric(read.table(global_seed_path, header = TRUE)$Seed[1])
  cat("🔁 Using global seed from:", global_seed_path, "→", global_seed, "\n")
}

# Set same seed in R and Python
set.seed(global_seed)
py_run_string(sprintf("import random, numpy as np; random.seed(%d); np.random.seed(%d)",
                      global_seed, global_seed))

# ── Local seed handling (cohort directory) ───────────────
seeds_path <- file.path(outputcohort, "Seeds.txt")

if (file.exists(seeds_path)) {
  cat("🔁 Reusing existing Seeds.txt for reproducibility:\n  ", seeds_path, "\n")
} else {
  cat("🧬 No Seeds.txt found — letting SigProfilerExtractor generate one automatically.\n")
}

# ── Run SigProfilerExtractor (Cohort mode) ───────────────
cat("\n🚀 Running SigProfilerExtractor (Cohort Mode)...\n")

SigProfilerExtractorR::sigprofilerextractor(
  input_type = "matrix",
  output = outputcohort,
  input_data = inputdatacohort,
  reference_genome = "GRCh38",
  opportunity_genome = "GRCh38",
  minimum_signatures = minimum_signatures,
  maximum_signatures = maximum_signatures,
  nmf_replicates = nmf_replicates,
  min_nmf_iterations = 5000,
  max_nmf_iterations = 100000,
  nmf_test_conv = 1000,
  nmf_tolerance = 1e-6
)

cat("\n✅ SigProfilerExtractor (Cohort Mode) finished successfully!\n")

