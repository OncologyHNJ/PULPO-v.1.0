#!/usr/bin/env Rscript
# ==================================================
# Script: CNV_formatting.R
# Description: Unified CNV formatter for OGM and WES data.
#              Converts raw CNV caller outputs into a normalized structure
#              compatible with different CNV signature extraction methods
#              (e.g. Tao, Drews, Steele, Sigminer).
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# ==================================================

# ───────────────────────────────────────────────
# 0️⃣ Command-line arguments
# ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
inputdata    <- args[1]
output_base  <- args[2]
patient      <- args[3]
method       <- args[4]                # "tao" | "drews"
source_type  <- ifelse(length(args) >= 5, args[5], "ogm")     # "ogm" | "wes"
caller       <- ifelse(length(args) >= 6, args[6], "ogm")     # "cnvkit" | "facets" | ...
genome_build <- ifelse(length(args) >= 7, args[7], "hg38")
purity_in    <- ifelse(length(args) >= 8, args[8], NA)
ploidy_in    <- ifelse(length(args) >= 9, args[9], NA)

suppressPackageStartupMessages({
  library(tools)
  library(dplyr)
  library(readr)
})

# ───────────────────────────────────────────────
# 1️⃣ Define output paths
# ───────────────────────────────────────────────
if (method == "tao") {
  output_dir  <- file.path(output_base, "tao")
  output_file <- "taoformatted.tsv"
} else if (method == "drews") {
  output_dir  <- file.path(output_base, "drews")
  output_file <- "drewsformatted.tsv"
} else {
  stop("❌ Method must be 'tao' or 'drews'")
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output <- file.path(output_dir, output_file)

# ───────────────────────────────────────────────
# 2️⃣ Utility functions
# ───────────────────────────────────────────────
norm_chr <- function(x) {
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("23", "X", "x")] <- "X"
  x[x %in% c("24", "Y", "y")] <- "Y"
  as.character(x)
}

as_num <- function(v) suppressWarnings(as.numeric(v))

# ───────────────────────────────────────────────
# 3️⃣ Readers for each data source
# ───────────────────────────────────────────────

# --- OGM CNV exports ---
read_ogm <- function(path_dir) {
  fps <- list.files(path = path_dir, pattern = "CNV", full.names = TRUE)
  if (length(fps) == 0) return(NULL)

  dfs <- list()
  for (f in fps) {
    ext <- file_ext(f)
    if (ext == "txt") {
      lines <- readLines(f, warn = FALSE)
      header_lines <- grep("^#", lines, value = TRUE)
      if (length(header_lines) > 0) {
        colnames_cnv <- strsplit(gsub("^#", "", tail(header_lines, 1)), "\t")[[1]]
        data_lines <- grep("^#", lines, value = TRUE, invert = TRUE)
        cnv <- read.table(text = data_lines, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        colnames(cnv) <- colnames_cnv
      } else {
        cnv <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      }
    } else if (ext == "csv") {
      cnv <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
    } else {
      next
    }
    dfs[[length(dfs) + 1]] <- cnv
  }

  if (!length(dfs)) return(NULL)
  bind_rows(dfs)
}

# --- CNVkit (.cns or .bed) ---
read_cnvkit <- function(path_file) {
  df <- suppressMessages(suppressWarnings(read_tsv(path_file, comment = "#", show_col_types = FALSE)))
  if (is.null(df) || nrow(df) == 0 || !any(nchar(names(df)))) {
    df <- read.delim(path_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }

  names(df) <- tolower(names(df))
  if ("chrom" %in% names(df)) names(df)[names(df) == "chrom"] <- "chromosome"
  if ("chr" %in% names(df))   names(df)[names(df) == "chr"]   <- "chromosome"
  if ("chromosome_start" %in% names(df)) names(df)[names(df) == "chromosome_start"] <- "start"
  if ("chromosome_end"   %in% names(df)) names(df)[names(df) == "chromosome_end"]   <- "end"

  if (!"cn" %in% names(df)) {
    if ("copy_number" %in% names(df)) names(df)[names(df) == "copy_number"] <- "cn"
    if ("copynumber"  %in% names(df)) names(df)[names(df) == "copynumber"]  <- "cn"
  }

  if (!"cn" %in% names(df) && "log2" %in% names(df)) {
    df$cn <- round(2 * 2^df$log2)
  }

  needed_any <- c("chromosome", "start", "end")
  if (!all(needed_any %in% names(df))) {
    stop("CNVkit missing mandatory columns (expected: chromosome/start/end). Found: ", paste(names(df), collapse = ", "))
  }

  df %>%
    mutate(chromosome = as.character(chromosome),
           start = as.integer(start),
           end = as.integer(end))
}

# --- Placeholder readers for other callers ---
read_facets   <- function(path_file) read.delim(path_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
read_ascat    <- function(path_file) read.delim(path_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
read_purple   <- function(path_file) read.delim(path_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
read_sequenza <- function(path_file) read.delim(path_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
read_gatk_gcnv <- function(path_file) read.delim(path_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ───────────────────────────────────────────────
# 4️⃣ Input dispatcher
# ───────────────────────────────────────────────
raw <- NULL
if (source_type == "ogm") {
  raw <- read_ogm(inputdata)
} else if (source_type == "wes") {
  if (caller == "cnvkit") {
    raw <- read_cnvkit(inputdata)
  } else if (caller == "purple") {
    raw <- read_purple(inputdata)
  } else if (caller == "facets") {
    raw <- read_facets(inputdata)
  } else if (caller == "ascat") {
    raw <- read_ascat(inputdata)
  } else if (caller == "sequenza") {
    raw <- read_sequenza(inputdata)
  } else if (caller == "gatk_gcnv") {
    raw <- read_gatk_gcnv(inputdata)
  } else {
    stop("Unsupported WES CNV caller: ", caller)
  }
} else {
  stop("Source must be 'ogm' or 'wes'")
}

if (is.null(raw) || nrow(raw) == 0) {
  dfcnv <- data.frame(chromosome = character(), start = numeric(),
                      end = numeric(), segVal = numeric(), sample = character())
  write.table(dfcnv, file = output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  quit(save = "no")
}

# ───────────────────────────────────────────────
# 5️⃣ CNV normalization per source/caller
# ───────────────────────────────────────────────
norm <- NULL

if (source_type == "ogm") {
  cols <- tolower(names(raw))
  chr <- if ("chromosome" %in% cols) raw[[which(cols == "chromosome")]] else raw$Chromosome
  st  <- if ("start" %in% cols) raw[[which(cols == "start")]] else raw$Start
  en  <- if ("end"   %in% cols) raw[[which(cols == "end")]]   else raw$End
  CN_abs  <- if ("copynumber" %in% cols) raw[[which(cols == "copynumber")]] else raw$CopyNumber
  CN_frac <- if ("fractionalcopynumber" %in% cols) raw[[which(cols == "fractionalcopynumber")]] else raw$fractionalCopyNumber

  norm <- data.frame(
    chromosome = norm_chr(chr),
    start = as.integer(round(as_num(st))),
    end   = as.integer(round(as_num(en))),
    cn_abs = as_num(CN_abs),
    cn_frac = as_num(CN_frac),
    stringsAsFactors = FALSE
  )

} else if (source_type == "wes" && caller == "cnvkit") {
  chr <- norm_chr(raw$chromosome)
  st  <- as.integer(raw$start)
  en  <- as.integer(raw$end)
  log2v <- if ("log2" %in% names(raw)) as_num(raw$log2) else NA
  cn_abs <- if ("cn" %in% names(raw)) as_num(raw$cn) else if (!all(is.na(log2v))) round(2 * 2^log2v) else NA

  ploidy <- suppressWarnings(as.numeric(ploidy_in))
  cn_frac <- if (!is.na(ploidy)) {
    cn_abs / ploidy
  } else if (!all(is.na(log2v))) {
    2^(log2v)
  } else if (!all(is.na(cn_abs))) {
    cn_abs / 2
  } else NA

  norm <- data.frame(
    chromosome = chr,
    start = st,
    end = en,
    cn_abs = cn_abs,
    cn_frac = cn_frac,
    stringsAsFactors = FALSE
  )

} else {
  stop("Caller normalization not yet implemented: ", caller)
}

# ───────────────────────────────────────────────
# 6️⃣ Compute segVal depending on method
# ───────────────────────────────────────────────
if (method == "tao") {
  segVal <- if (!all(is.na(norm$cn_abs))) {
    norm$cn_abs
  } else if (!all(is.na(norm$cn_frac))) {
    round(norm$cn_frac * 2)
  } else NA
} else { # drews
  segVal <- if (!all(is.na(norm$cn_frac))) {
    norm$cn_frac
  } else if (!all(is.na(norm$cn_abs))) {
    norm$cn_abs / ifelse(is.na(as.numeric(ploidy_in)), 2, as.numeric(ploidy_in))
  } else NA
}

# ───────────────────────────────────────────────
# 7️⃣ Build final output table
# ───────────────────────────────────────────────
out_df <- data.frame(
  sample = rep(patient, nrow(norm)),
  chromosome = norm$chromosome,
  start = norm$start,
  end = norm$end,
  segVal = as_num(segVal),
  stringsAsFactors = FALSE
)

# Basic cleanup and sorting
out_df <- out_df %>%
  filter(!is.na(segVal), !is.na(start), !is.na(end), start < end) %>%
  arrange(chromosome, start)

write.table(out_df, file = output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("✅ CNV formatting completed for", patient, "→", output, "\n")

