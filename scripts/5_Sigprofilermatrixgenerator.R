#!/usr/bin/env Rscript
# ==================================================
# Script: 5_Sigprofilermatrixgenerator.R
# Description: Generates SV mutational matrices using
#              SigProfilerMatrixGeneratorR for individual samples.
# Author: Marta Portasany
# Created: 2025-02-27
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: reticulate, SigProfilerMatrixGeneratorR
# ==================================================

# ───────────────────────────────────────────────
# 0️⃣ LOAD ARGUMENTS
# ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
inputdata       <- args[1]
pythondirectory <- args[2]
output          <- args[3]
patientid <- as.character(args[4])

# ───────────────────────────────────────────────
# 1️⃣ INSTALL/ LOAD LIBRARIES
# ───────────────────────────────────────────────
suppressPackageStartupMessages({
  ## 0) Get Python path from pythondirectory argument
  python_bin     <-  pythondirectory
  
  ## 1) reticulate: expected to be installed in the conda env
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("❌ 'reticulate' is not installed in this environment.\n",
         "Please install it in the PULPO2 conda environment, e.g.:\n",
         "  conda install -c conda-forge r-reticulate",
         call. = FALSE)
  }
  library(reticulate)
  
  ## 2) Force reticulate to use the conda Python (no uv env)
  cat("🔹 Using Python from:", python_bin, "\n")
  reticulate::use_python(python_bin, required = TRUE)
  
  ## 3) Check that the Python module is available
  if (!reticulate::py_module_available("SigProfilerMatrixGenerator")) {
    stop(
      "❌ Python module 'SigProfilerMatrixGenerator' is not available in this Python:\n  ",
      python_bin, "\n",
      "Please install it in the PULPO2 env, e.g.:\n",
      "  conda activate PULPO2\n",
      "  pip install SigProfilerMatrixGenerator",
      call. = FALSE
    )
  }
  
  ## 4) devtools: expected to be pre-installed as well
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("❌ 'devtools' is not installed in this environment.\n",
         "Please install it in the PULPO2 conda environment, e.g.:\n",
         "  conda install -c conda-forge r-devtools",
         call. = FALSE)
  }
  library(devtools)
  
  ## 5) SigProfilerMatrixGeneratorR: install from GitHub if missing
  tryCatch(
    {
      library(SigProfilerMatrixGeneratorR)
      cat("✅ SigProfilerMatrixGeneratorR loaded successfully.\n")
    },
    error = function(e) {
      cat("⚠️ SigProfilerMatrixGeneratorR not found — installing from GitHub...\n")
      
      devtools::install_github(
        "AlexandrovLab/SigProfilerMatrixGeneratorR",
        build_vignettes = FALSE,
        dependencies   = TRUE
      )
      
      library(SigProfilerMatrixGeneratorR)
      cat("✅ SigProfilerMatrixGeneratorR installed and loaded.\n")
    }
  )
})
# ───────────────────────────────────────────────
# 2️⃣ HELPER FUNCTION
# ───────────────────────────────────────────────
is_empty_or_header_only <- function(file) {
  lines <- readLines(file, warn = FALSE)
  return(length(lines) == 1)  # only header line → no variants
}

# ───────────────────────────────────────────────
# 3️⃣ VALIDATE INPUT FILE
# ───────────────────────────────────────────────
if (!file.exists(inputdata) || file.info(inputdata)$size == 0 || is_empty_or_header_only(inputdata)) {
  cat("⚠️  Input file is empty or header-only. Creating placeholder output and skipping analysis.\n")

  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }

  empty_output_file <- file.path(output, paste0(patientid, ".SV32.matrix.tsv"))
  file.create(empty_output_file)
  cat("✅ Empty matrix placeholder created at:", empty_output_file, "\n")

  quit(status = 0)
}

# ───────────────────────────────────────────────
# 4️⃣ RUN SIGPROFILERMATRIXGENERATOR
# ───────────────────────────────────────────────
tryCatch({
  dirpath <- dirname(inputdata)
  SigProfilerMatrixGeneratorR::SVMatrixGenerator(
    input_dir = dirpath,
    project = patientid,
    output_dir = output
  )
  cat("✅ Matrix generation completed successfully for:", patientid, "\n")
}, error = function(e) {
  message("❌ Error during matrix generation for ", patientid, ": ", e$message)
  quit(status = 1)
})

cat("✔️ Finished SigProfilerMatrixGenerator execution.\n")

