#!/usr/bin/env Rscript
# ==================================================
# Script: 12_Mergecohortcnvs.R
# Description: Merges all individual CNV48 matrix files
#              into a single cohort-level matrix (TSV format)
# Author: Marta Portasany-Rodríguez
# Created on: 2025-02-27
# Last modified: 2025-11-05
# Pipeline: PULPO 🐙 
# Dependencies: Base R (no external packages required)
# ==================================================

# ── 1️⃣ Parse command-line arguments ───────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
directorypatients <- args[1]
outputcohort      <- args[2]

# Example (for testing):
# directorypatients <- "/home/user/MARTA/LABO/PULPO/results/CNVs/Patients"
# outputcohort      <- "/home/user/MARTA/LABO/PULPO/results/CNVs/Cohort/SigProfiler/MatrixGenerator/Cohort.CNV48.matrix.tsv"

cat("🔹 Merging CNV48 matrices from directory:", directorypatients, "\n")

# ── 2️⃣ List patient directories and target matrix files ───────────
patients <- list.files(directorypatients)
directories <- file.path(directorypatients, patients,
                         "SigProfiler/results/MatrixGenerator",
                         paste0(patients, ".CNV48.matrix.tsv"))

cat("📂 Found", length(directories), "potential CNV matrix files.\n")

# ── 3️⃣ Read valid matrix files (non-empty) ─────────────────────────
files <- list()

for (path in directories) {
  if (file.exists(path)) {
    if (file.info(path)$size > 0) {
      cat("✅ Reading file:", path, "\n")
      file <- read.table(path, header = TRUE)
      files <- append(files, list(file))
    } else {
      message("⚠️  The file ", path, " is empty and will be skipped.")
    }
  } else {
    message("⚠️  The file ", path, " does not exist and will be skipped.")
  }
}

if (length(files) == 0) {
  stop("❌ No valid CNV48 matrices found in directory: ", directorypatients)
}

# ── 4️⃣ Merge matrices column-wise ─────────────────────────────────
cat("🔄 Combining", length(files), "CNV matrices into a single cohort file...\n")

cohortfile <- files[[1]][1]  # Initialize with MutationType column

for (i in seq_along(files)) {
  current <- files[[i]]
  if (nrow(current) > 0 && ncol(current) > 1) {
    cohortfile <- cbind(cohortfile, current[2])
  } else {
    message("⚠️  Skipping file ", i, " — invalid or single-column structure.")
  }
}

# ── 5️⃣ Save merged cohort matrix ──────────────────────────────────
dir.create(dirname(outputcohort), recursive = TRUE, showWarnings = FALSE)
write.table(cohortfile, file = outputcohort, sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("✅ Cohort CNV48 matrix successfully written to:\n  ", outputcohort, "\n")

