#!/usr/bin/env Rscript
# ==================================================
# Script: 10_Sigprofilermatrixgeneratorcnvs.R
# Description: Executes SigProfilerMatrixGeneratorR for individual
#              Copy Number Variant (CNV) files (single-sample analysis).
# Author: Marta Portasany
# Created: 2025-02-27
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: reticulate, devtools, SigProfilerMatrixGeneratorR
# ==================================================

# ───────────────────────────────────────────────
# 0️⃣ LOAD ARGUMENTS
# ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
inputdata       <- args[1]
pythondirectory <- args[2]
output          <- args[3]
patient         <- args[4]
sigprofilerfile <- args[5]

cat("🔹 Running CNV Matrix Generator for sample:", patient, "\n")
cat("   Input:", inputdata, "\n")
cat("   Output:", output, "\n\n")

# ───────────────────────────────────────────────
# 1️⃣ LOAD LIBRARIES AND PYTHON ENVIRONMENT
# ───────────────────────────────────────────────
suppressPackageStartupMessages({
  library(reticulate)
  library(devtools)
})

use_python(pythondirectory)
py_config()
py_run_string("import sys")

# ───────────────────────────────────────────────
# 2️⃣ ENSURE SigProfilerMatrixGeneratorR IS AVAILABLE
# ───────────────────────────────────────────────
suppressPackageStartupMessages({
  ## reticulate (por si acaso)
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cat("🔹 Installing 'reticulate' from CRAN...\n")
    install.packages("reticulate", repos = "https://cloud.r-project.org")
  }
  library(reticulate)
  
  ## SigProfilerMatrixGeneratorR: intenta cargar, si falla instala desde GitHub
  tryCatch(
    {
      library(SigProfilerMatrixGeneratorR)
      cat("✅ SigProfilerMatrixGeneratorR loaded successfully.\n")
    },
    error = function(e) {
      cat("⚠️ Package not found — installing SigProfilerMatrixGeneratorR from GitHub...\n")
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "https://cloud.r-project.org")
      }
      devtools::install_github("AlexandrovLab/SigProfilerMatrixGeneratorR")
      library(SigProfilerMatrixGeneratorR)
      cat("✅ SigProfilerMatrixGeneratorR installed and loaded.\n")
    }
  )
})

# ───────────────────────────────────────────────
# 3️⃣ READ INPUT FILE AND CHECK CONTENT
# ───────────────────────────────────────────────
tryCatch({
  cnv_data <- read.table(inputdata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  if (nrow(cnv_data) == 0) {
    cat("⚠️ CNV file is empty. Creating placeholder output file...\n")

    if (!dir.exists(output)) {
      dir.create(output, recursive = TRUE)
      cat("📂 Created directory:", output, "\n")
    }

    file.create(sigprofilerfile)
    cat("✅ Placeholder CNV matrix created at:", sigprofilerfile, "\n")

  } else {
    # ───────────────────────────────────────────────
    # 4️⃣ RUN SigProfilerMatrixGeneratorR
    # ───────────────────────────────────────────────
    cat("🔹 CNV file contains", nrow(cnv_data), "rows. Running matrix generator...\n")

    tryCatch({
      SigProfilerMatrixGeneratorR::CNVMatrixGenerator(
        file_type = "PCAWG",
        input_file = inputdata,
        project = patient,
        output_path = output
      )
      cat("✅ CNV Matrix successfully generated for", patient, "\n")
    },
    error = function(e) {
      if (grepl("objeto 'sys' no encontrado", e$message)) {
        cat("⚠️ Warning: 'sys not found' error ignored.\n")
      } else {
        stop(e)
      }
    })
  }

}, error = function(e) {
  cat("❌ ERROR:", e$message, "\n")
  cat("Creating empty output for patient:", patient, "\n")

  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
    cat("📂 Created directory:", output, "\n")
  }

  file.create(sigprofilerfile)
  cat("✅ Placeholder file created at:", sigprofilerfile, "\n")
})

cat("✔️ Finished CNV Matrix Generation for sample:", patient, "\n")

