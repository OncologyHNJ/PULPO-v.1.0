#!/usr/bin/env Rscript
# ==================================================
# Script: ogm_cnv_to_sigprofiler.R
# Description: Converts OGM CNV export files (.txt/.csv) into
#              SigProfiler-compatible CNV tables (TSV format).
# Author: Marta Portasany (v2025-10-31)
# Pipeline: PULPO
# ==================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ───────────────────────────────────────────────
# 0️⃣ Command-line arguments and log setup
# ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[1]
output_file <- args[2]
platform <- "OGM"

cat("🔹 Processing CNV file:", input_file, "\n")
cat("🔹 Output path:", output_file, "\n")
cat("🔹 Platform mode:", platform, "\n")

cnv_dir <- dirname(input_file)
done_log  <- file.path(cnv_dir, "check_cnvs_done.txt")
error_log <- file.path(cnv_dir, "check_cnvs_error.txt")

if (file.exists(done_log))  file.remove(done_log)
if (file.exists(error_log)) file.remove(error_log)

# ───────────────────────────────────────────────
# 1️⃣ READ CNV FILE (.txt or .csv)
# ───────────────────────────────────────────────
if (!file.exists(input_file)) {
  msg <- paste(Sys.time(), "❌ CNV file not found:", input_file)
  writeLines(msg, error_log)
  quit(status = 0)
}

if (file.info(input_file)$size == 0) {
  msg <- paste(Sys.time(), "⚠️ CNV file is empty:", input_file)
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
  quit(status = 0)
}

is_csv <- grepl("\\.csv$", input_file, ignore.case = TRUE)
is_txt <- grepl("\\.txt$", input_file, ignore.case = TRUE)

cnv_raw <- tryCatch({
  if (is_csv) {
    read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
  } else if (is_txt) {
    # detect last header line (starts with "#Id")
    lines <- readLines(input_file)
    header_line <- grep("^#Id", lines, value = TRUE)
    if (length(header_line) == 0)
      stop("No '#Id' header line found in file.")

    column_names <- strsplit(gsub("^#", "", header_line), "\t")[[1]]
    column_names[1] <- gsub("\ufeff", "", column_names[1]) # remove BOM if any

    # read table skipping comment lines
    dat <- read.table(
      input_file, header = FALSE, sep = "\t", comment.char = "#",
      stringsAsFactors = FALSE, quote = "", fill = TRUE
    )
    colnames(dat)[seq_along(column_names)] <- column_names
    dat
  } else {
    read_tsv(input_file, show_col_types = FALSE)
  }
}, error = function(e) {
  msg <- paste(Sys.time(), "❌ Error reading CNV file:", conditionMessage(e))
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
  quit(status = 0)
})

if (nrow(cnv_raw) == 0) {
  msg <- paste(Sys.time(), "⚠️ CNV file has headers but no data:", input_file)
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
  quit(status = 0)
}

cat("✅ CNV file loaded with", nrow(cnv_raw), "rows and", ncol(cnv_raw), "columns.\n")

# ───────────────────────────────────────────────
# 2️⃣ BUILD CNV TABLE
# ───────────────────────────────────────────────
cols_lower <- tolower(colnames(cnv_raw))
Chromosome <- cnv_raw[[which(cols_lower == "chromosome")[1]]]
Start       <- cnv_raw[[which(cols_lower == "start")[1]]]
End         <- cnv_raw[[which(cols_lower == "end")[1]]]
CopyNumber  <- cnv_raw[[which(cols_lower == "copynumber")[1]]]
Type        <- cnv_raw[[which(cols_lower == "type")[1]]]

cnv <- data.frame(
  chromosome       = Chromosome,
  chromosome_start = round(as.numeric(Start)),
  chromosome_end   = round(as.numeric(End)),
  copy_number      = as.integer(CopyNumber),
  mutation_type    = gsub("(loss|gain)_masked", "\\1", Type),
  stringsAsFactors = FALSE
)

# ───────────────────────────────────────────────
# 3️⃣ NORMALIZE FOR OGM
# ───────────────────────────────────────────────
cnv <- cnv %>%
  mutate(
    mutation_type = case_when(
      copy_number == 0 ~ "loss",
      copy_number == 1 ~ "loss",
      copy_number == 2 ~ "copy neutral",
      copy_number >= 3 ~ "gain",
      TRUE ~ mutation_type
    )
  )

# ───────────────────────────────────────────────
# 4️⃣ ADD SAMPLE NAME (ANONYMISED)
# ───────────────────────────────────────────────
path_parts  <- strsplit(cnv_dir, "/")[[1]]
sample_name <- tail(path_parts[grep("Patient", path_parts)], 1)

cnv <- cnv %>%
  mutate(sample = sample_name) %>%
  select(sample, chromosome, chromosome_start, chromosome_end, copy_number, mutation_type)

# ───────────────────────────────────────────────
# 5️⃣ WRITE OUTPUT
# ───────────────────────────────────────────────
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (nrow(cnv) == 0) {
  msg <- paste(Sys.time(), "⚠️ No valid CNVs after filtering:", input_file)
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
} else {
  write_tsv(cnv, output_file)
  cat("✅ CNV formatted and written to:", output_file, "\n")
  writeLines(paste(Sys.time(), "OK", sample_name), done_log)
}

cat("✔️ Finished CNV formatting for", sample_name, "\n")

