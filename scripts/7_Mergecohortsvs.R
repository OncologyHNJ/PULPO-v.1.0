#!/usr/bin/env Rscript
# ==================================================
# Script: 7_Mergecohortsvs.R
# Description: Merges per-patient SigProfiler SV matrices (SV32) 
#              into a single cohort-level matrix for downstream 
#              signature extraction.
# Author: Marta Portasany
# Created on: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: base R (read.table, write.table)
# ==================================================

# ── Command-line arguments ───────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
directorypatients <- args[1]
outputcohort <- args[2]

cat("🔹 Merging SigProfiler SV matrices from:", directorypatients, "\n")
cat("🔹 Output cohort matrix will be saved to:", outputcohort, "\n\n")

# ── List all anonymised patient directories ───────────────
patients <- list.files(directorypatients, full.names = FALSE)
directories <- character()
files <- list()

# ── Build full paths for each patient matrix ──────────────
for (p in patients) {
  path <- file.path(directorypatients, p, "SigProfiler", "results", "MatrixGenerator", paste0(p, ".SV32.matrix.tsv"))
  directories <- c(directories, path)
}

cat("📂 Found", length(directories), "potential patient matrix files.\n")

# ── Read available matrices into a list ───────────────────
for (d in directories) {
  if (file.exists(d)) {
    df <- tryCatch(read.table(d, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                   error = function(e) { warning(paste("⚠️ Could not read:", d)); NULL })
    if (!is.null(df)) files <- append(files, list(df))
  } else {
    warning(paste("⚠️ File does not exist:", d))
  }
}

# ── Check if any files were successfully loaded ───────────
if (length(files) == 0) {
  stop("❌ No valid SV32 matrix files found in:", directorypatients)
}

cat("✅ Successfully loaded", length(files), "SV matrix files.\n")

# ── Merge all patient matrices by column ──────────────────
cohortfile <- files[[1]][1]  # Start with MutationType column

for (f in seq_along(files)) {
  df <- files[[f]]
  if (nrow(df) > 0 && ncol(df) > 1) {
    cohortfile <- cbind(cohortfile, df[2])
  } else {
    message("⚠️ File", f, "was empty or malformed → skipped.")
  }
}

# ── Write merged cohort file ──────────────────────────────
dir.create(dirname(outputcohort), recursive = TRUE, showWarnings = FALSE)
write.table(cohortfile, file = outputcohort, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n🎯 Cohort-level SV32 matrix successfully written to:\n", outputcohort, "\n")
cat("🧩 Total patients merged:", ncol(cohortfile) - 1, "\n")
cat("✅ Merge completed.\n")

