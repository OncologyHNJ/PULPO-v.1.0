#!/usr/bin/env Rscript
# ==================================================
# Script: CheckfilesSV.R
# Description: Checks the integrity of OGM .smap files by verifying
#              that each file contains variant entries (non-comment lines).
# Author: Marta Portasany (v2025-11-05)
# ==================================================

# ───────────────────────────────────────────────
# 0️⃣ LOAD ARGUMENTS AND DEFINE LOG PATHS
# ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
patient   <- args[1]
inputsv   <- args[2]
error_log <- file.path(inputsv, "check_svs_error.txt")

# ───────────────────────────────────────────────
# 1️⃣ ERROR HANDLING FUNCTION
# ───────────────────────────────────────────────
handle_error <- function(e) {
  error_message <- paste0("❌ ERROR: ", e$message, " in sample ", patient, ".\n")
  writeLines(error_message, con = error_log)
  stop(error_message)
}

# ───────────────────────────────────────────────
# 2️⃣ MAIN EXECUTION WITH ERROR CAPTURE
# ───────────────────────────────────────────────
tryCatch({

  cat("🔹 Checking OGM .smap files in:", inputsv, "\n")

  # Search for .smap files within the directory
  smap_files <- list.files(inputsv, pattern = "\\.smap$", full.names = TRUE)
  
  # Validate the presence of .smap files
  if (length(smap_files) == 0) {
    stop("No .smap files were found in the directory.")
  }
  
  # Process each .smap file found
  for (smap_file in smap_files) {
    all_lines <- readLines(smap_file)
    
    # Filter lines that do not start with '#'
    variant_lines <- grep("^[^#]", all_lines, value = TRUE)
    
    # Check for variant content
    if (length(variant_lines) == 0) {
      stop(paste0("The file ", smap_file, " for sample ", patient, " is empty."))
    } else {
      message(paste0("✅ OK: The file ", smap_file, " for sample ", patient, " contains variants."))
    }
  }

  cat("✔️ Validation completed successfully for sample:", patient, "\n")

}, error = handle_error)

