#!/usr/bin/env Rscript
# ==================================================
# Script: ogm_to_sigprofiler.R
# Description: Converts OGM SV export files (.smap) into
#              SigProfiler-compatible SV format (bedpe).
# Author: Marta Portasany (v2025-10-31)
# Pipeline: PULPO
# ==================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript ogm_to_sigprofiler.R <bionanodirectory> <patient_id> <out_sigprofiler_bedpe>")
}
bionanodir <- args[1]          
patient_id <- args[2]          
outfile    <- args[3]         



suppressPackageStartupMessages({
  library(dplyr)
})

# local helpers
mk <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# 1) localizar .smap como haces tú (prioriza Rare_Variant_Analysis)
ogm_dir <- file.path(bionanodir, patient_id, "OGMdata")
files <- if (dir.exists(ogm_dir)) list.files(ogm_dir, pattern="\\.smap$", full.names=FALSE) else character(0)
if (length(files) == 0) {
  message("No .smap para ", patient_id, "; generando salida vacía.")
  mk(dirname(outfile))
  write.table(data.frame(chrom1=character(), start1=integer(), end1=integer(),
                         chrom2=character(), start2=integer(), end2=integer(),
                         sample=character(), svclass=character()),
              file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
  quit(save="no")
}
preferred <- files[grepl("Rare_Variant_Analysis", files)]
smap_file <- if (length(preferred)>0) preferred[1] else files[1]
full_path <- file.path(ogm_dir, smap_file)

# 2) leer .smap "a tu manera" (header = FALSE)
df <- tryCatch(read.table(full_path, header=FALSE),
               error=function(e) NULL)

mk(dirname(outfile))
if (is.null(df) || nrow(df)==0) {
  write.table(data.frame(chrom1=character(), start1=integer(), end1=integer(),
                         chrom2=character(), start2=integer(), end2=integer(),
                         sample=character(), svclass=character()),
              file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
  quit(save="no")
}

# Construcción según tu createbedpe (end = start tal cual)
chrom1 <- df$V3
start1 <- round(as.numeric(df$V7))
end1   <- round(as.numeric(df$V7))
chrom2 <- df$V4
start2 <- round(as.numeric(df$V8))
end2   <- round(as.numeric(df$V8))
sample <- patient_id

# OJO: aquí asumo que la clase viene en V10 tal y como mostraste en tu ejemplo.
# Si en tus .smap reales la clase está en otra columna, cambia df$V10 por la que toque.
sv_raw <- tolower(gsub("_","-", as.character(df$V10)))

# Normalización → CLASES LARGAS
sv_long <- dplyr::case_when(
  sv_raw %in% c("translocation-interchr","translocation-intrachr","translocation") ~ "translocation",
  grepl("^dup", sv_raw) | sv_raw %in% c("duplication","duplication-inverted","duplication-split","tandem-duplication") ~ "tandem-duplication",
  grepl("^inv", sv_raw) | sv_raw %in% c("inversion","inversion-paired","inversion-partial") ~ "inversion",
  grepl("^del", sv_raw) | sv_raw %in% c("deletion","del") ~ "deletion",
  TRUE ~ sv_raw
)

# Filter unwanted SV classes
drop_vals <- c("insertion", "insertion-nbase", "insertion_nbase", "inversion-partial")
keep <- !(sv_long %in% drop_vals)

# Build initial filtered data frame
df_temp <- data.frame(
  chrom1 = as.character(chrom1[keep]),
  start1 = as.integer(start1[keep]),
  end1   = as.integer(end1[keep]),
  chrom2 = as.character(chrom2[keep]),
  start2 = as.integer(start2[keep]),
  end2   = as.integer(end2[keep]),
  sample = rep(sample, sum(keep)),
  svclass = as.character(sv_long[keep]),
  stringsAsFactors = FALSE
)

# --- Conditional filter for inversions with -1 coordinates ---
invalid_inversions <- subset(df_temp, svclass == "inversion" & (chrom2 == "-1" | start2 == -1 | end2 == -1))
n_invalid <- nrow(invalid_inversions)

if (n_invalid > 5) {
  message("⚠️  ", n_invalid, " inversions with -1 coordinates detected → they will be removed")
  df_temp <- subset(df_temp, !(svclass == "inversion" & (chrom2 == "-1" | start2 == -1 | end2 == -1)))
} else if (n_invalid > 0) {
  message("✅  Only ", n_invalid, " inversions with -1 coordinates → they will be kept")
} else {
  message("✅  No inversions with -1 coordinates detected")
}

# Final output
out <- df_temp

write.table(out, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
message("OK: ", patient_id, " → ", outfile)
