#!/usr/bin/env Rscript
# ==================================================
# Script: 11_Sigprofilerextractorcnvs.R
# Description: Runs SigProfilerExtractorR for individual
#              Copy Number Variant (CNV) mutational signature analysis.
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
inputdata      <- args[1]
pythondirectory <- args[2]
outputdata     <- args[3]
min_signatures <- as.numeric(args[4])
max_signatures <- as.numeric(args[5])
nmf_replicates <- as.numeric(args[6])
errorlog       <- args[7]

cat("🔹 Starting CNV signature extraction for input:", inputdata, "\n")

# ───────────────────────────────────────────────
# 1️⃣ LOAD LIBRARIES
# ───────────────────────────────────────────────
suppressPackageStartupMessages({
  library(reticulate)
  library(devtools)
})

# ───────────────────────────────────────────────
# 2️⃣ INSTALL / LOAD SigProfilerExtractorR
# ───────────────────────────────────────────────
if (!requireNamespace("SigProfilerExtractorR", quietly = TRUE)) {
  cat("⚠️ Installing SigProfilerExtractorR package from GitHub...\n")
  install.packages("devtools", repos = "http://cran.us.r-project.org")
  devtools::install_github("AlexandrovLab/SigProfilerExtractorR")
  library(SigProfilerExtractorR)
  SigProfilerExtractorR::install("GRCh38", rsync = FALSE, bash = TRUE)
} else {
  cat("✅ SigProfilerExtractorR already available.\n")
  library(SigProfilerExtractorR)
}

# ───────────────────────────────────────────────
# 3️⃣ CONFIGURE PYTHON ENVIRONMENT
# ───────────────────────────────────────────────
use_python(pythondirectory)
py_config()

# ───────────────────────────────────────────────
# 4️⃣ INITIALIZE ERROR LOG
# ───────────────────────────────────────────────
log_dir <- dirname(errorlog)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
file.create(errorlog)

# ───────────────────────────────────────────────
# 5️⃣ CHECK INPUT FILE
# ───────────────────────────────────────────────
if (!file.exists(inputdata)) {
  cat("⚠️ CNV input file not found. Writing path to error log and exiting.\n")
  dir.create(outputdata, recursive = TRUE, showWarnings = FALSE)
  writeLines(inputdata, errorlog, sep = "\n")
  quit(status = 0)
}

# ───────────────────────────────────────────────
# 6️⃣ GLOBAL REPRODUCIBILITY SEED
# ───────────────────────────────────────────────
global_seed_path <- file.path(dirname(dirname(outputdata)), "Seeds_master.txt")

if (!file.exists(global_seed_path)) {
  global_seed <- as.integer(runif(1, min = 1, max = 2^31 - 1))
  write.table(
    data.frame(Seed = global_seed),
    file = global_seed_path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  cat("🧬 Created global seed file:", global_seed_path, "→ Seed:", global_seed, "\n")
} else {
  global_seed <- as.numeric(read.table(global_seed_path, header = TRUE)$Seed[1])
  cat("🔁 Using global seed from:", global_seed_path, "→", global_seed, "\n")
}

# Apply seed to both R and Python
set.seed(global_seed)
py_run_string(sprintf("import random, numpy as np; random.seed(%d); np.random.seed(%d)", global_seed, global_seed))

# ───────────────────────────────────────────────
# 7️⃣ SEED HANDLING (LOCAL)
# ───────────────────────────────────────────────
seeds_path <- file.path(outputdata, "Seeds.txt")
if (file.exists(seeds_path)) {
  cat("🔁 Reusing existing Seeds.txt:", seeds_path, "\n")
} else {
  cat("🧬 No Seeds.txt found — letting SigProfilerExtractor generate one automatically.\n")
}

# ───────────────────────────────────────────────
# 8️⃣ RUN SigProfilerExtractorR
# ───────────────────────────────────────────────
cat("🔹 Running SigProfilerExtractorR with parameters:\n",
    "   Input:", inputdata, "\n",
    "   Output:", outputdata, "\n",
    "   Min signatures:", min_signatures, "\n",
    "   Max signatures:", max_signatures, "\n",
    "   NMF replicates:", nmf_replicates, "\n\n"
)

tryCatch({
  SigProfilerExtractorR::sigprofilerextractor(
    input_type = "matrix",
    output = outputdata,
    input_data = inputdata,
    reference_genome = "GRCh38",
    minimum_signatures = min_signatures,
    maximum_signatures = max_signatures,
    nmf_replicates = nmf_replicates,
    min_nmf_iterations = 1000,
    max_nmf_iterations = 100000,
    nmf_test_conv = 1000,
    nmf_tolerance = 1e-8
  )
  cat("✅ SigProfilerExtractorR finished successfully!\n")

}, error = function(e) {
  cat("❌ ERROR during CNV signature extraction:", e$message, "\n")
  writeLines(paste("ERROR:", e$message), con = errorlog)
  quit(status = 1)
})

cat("✔️ Completed CNV signature extraction for:", inputdata, "\n")

