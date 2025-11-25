# ==================================================
# Script: ogm_cnv_to_sigprofiler.R
# Description: Validate and format CNV files to SigProfiler input (handles .csv/.txt and rounds Start/End)
# Author: Marta Portasany (updated 2025-10-28)
# ==================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[1]
output_file <- args[2]

cat("🔹 Processing CNV file:", input_file, "\n")
cat("🔹 Output path:", output_file, "\n")

# Define paths for sentinel logs
cnv_dir <- dirname(input_file)
done_log <- file.path(cnv_dir, "check_cnvs_done.txt")
error_log <- file.path(cnv_dir, "check_cnvs_error.txt")

# Clean up previous logs if they exist
if (file.exists(done_log)) file.remove(done_log)
if (file.exists(error_log)) file.remove(error_log)

# -------------------------------
# 1️⃣ VALIDATION
# -------------------------------
if (!file.exists(input_file)) {
  msg <- paste(Sys.time(), "❌ ERROR: CNV file not found at", input_file)
  writeLines(msg, error_log)
  file.create(output_file)
  quit(status = 0)
}

if (file.info(input_file)$size == 0) {
  msg <- paste(Sys.time(), "⚠️ WARNING: CNV file is empty at", input_file)
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
  writeLines(paste(Sys.time(), "EMPTY CNV FILE"), done_log)
  quit(status = 0)
}

# -------------------------------
# 2️⃣ LOAD CNV DATA (robust for both TXT and CSV)
# -------------------------------
cnv_raw <- tryCatch({
  if (grepl("\\.csv$", input_file, ignore.case = TRUE)) {
    # Caso 1: extensión .csv → separador coma
    read_csv(input_file, show_col_types = FALSE)
  } else if (grepl("\\.txt$", input_file, ignore.case = TRUE)) {
    # Caso 2: extensión .txt → autodetectar delimitador
    first_line <- readLines(input_file, n = 1)
    if (grepl(",", first_line)) {
      read_csv(input_file, show_col_types = FALSE)
    } else if (grepl(";", first_line)) {
      read_delim(input_file, delim = ";", show_col_types = FALSE)
    } else {
      # Por defecto, tabulado
      suppressMessages(data.table::fread(input_file, data.table = FALSE))
    }
  } else {
    # Caso general o desconocido
    suppressMessages(data.table::fread(input_file, data.table = FALSE))
  }
}, error = function(e) {
  msg <- paste(Sys.time(), "❌ ERROR reading CNV file:", conditionMessage(e))
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
  writeLines(paste(Sys.time(), "READ ERROR"), done_log)
  quit(status = 0)
})

if (nrow(cnv_raw) == 0) {
  msg <- paste(Sys.time(), "⚠️ CNV file has headers but no data:", input_file)
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
  writeLines(paste(Sys.time(), "EMPTY AFTER READ"), done_log)
  quit(status = 0)
}

cat("✅ CNV file loaded with", nrow(cnv_raw), "rows and", ncol(cnv_raw), "columns.\n")

# -------------------------------
# 3️⃣ FORMAT CNV DATA
# -------------------------------
possible_cols <- tolower(names(cnv_raw))

chrom_col <- names(cnv_raw)[which(possible_cols %in% c("chromosome", "chr", "chrom", "refcontig"))[1]]
start_col <- names(cnv_raw)[which(possible_cols %in% c("start", "pos1", "begin", "startposition"))[1]]
end_col   <- names(cnv_raw)[which(possible_cols %in% c("end", "pos2", "stop", "endposition"))[1]]
cn_col    <- names(cnv_raw)[which(possible_cols %in% c("copynumber", "cn", "totalcn", "cn_total"))[1]]
type_col  <- names(cnv_raw)[which(possible_cols %in% c("type", "eventtype", "varianttype", "cnvtype", "cnv type"))[1]]

cnv <- cnv_raw %>%
  mutate(
    Chromosome = if (!is.null(chrom_col)) .data[[chrom_col]] else NA,
    Start      = if (!is.null(start_col)) .data[[start_col]] else NA,
    End        = if (!is.null(end_col))   .data[[end_col]]   else NA,
    CNV_Type   = if (!is.null(type_col))  .data[[type_col]]  else NA,
    CopyNumber = if (!is.null(cn_col))    .data[[cn_col]]    else NA
  ) %>%
  # 🔹 Convertir Start/End a enteros
  mutate(
    Start = suppressWarnings(as.numeric(Start)),
    End   = suppressWarnings(as.numeric(End)),
    Start = round(Start),
    End   = round(End)
  ) %>%
  # 🔹 Filtrar filas válidas
  filter(!is.na(Start) & !is.na(End) & !is.na(Chromosome)) %>%
  select(Chromosome, Start, End, CNV_Type, CopyNumber)

# Normalizar nombres de cromosoma
cnv <- cnv %>%
  mutate(
    Chromosome = as.character(Chromosome),
    Chromosome = str_replace(Chromosome, "^([0-9XYM])$", "chr\\1"),
    Chromosome = ifelse(!str_detect(Chromosome, "^chr"), paste0("chr", Chromosome), Chromosome)
  )

# -------------------------------
# 4️⃣ WRITE OUTPUT + SENTINELS
# -------------------------------
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (nrow(cnv) == 0) {
  msg <- paste(Sys.time(), "⚠️ No valid CNV rows after filtering:", input_file)
  writeLines(msg, error_log)
  write_tsv(tibble(), output_file)
  writeLines(paste(Sys.time(), "NO VALID CNVS"), done_log)
} else {
  write_tsv(cnv, output_file)
  cat("✅ CNV formatted and written to:", output_file, "\n")
  writeLines(paste(Sys.time(), "OK", basename(input_file)), done_log)
}

cat("✔️ Finished CNV formatting for", basename(input_file), "\n")

