# ==================================================
# Script: 4_FORMATBEDPE.R
# Description: First step to formatting DATA needed to SigProfiler
# Author: Marta Portasany
# Created on: 2025-02-27
# Last modified: 2025-06-13
# Pipeline: PULPO
# Dependencies: stringr, readr
#########################################
args <- commandArgs(trailingOnly = TRUE)
bionanodirectory <- args[1]
outputdirbaseSVs <- args[2]
############ LOADS ########################
library(gtools)
library(doParallel)
library(R.utils)
library(magrittr)
library(tidyverse)
###########################################
# bionanodirectory <- "/home/user/MARTA/LABO/PULPO/PULPO_synthetic/results/DATA/Patients/"
# file_path <- bionanodirectory
# outputdirbaseSVs <- "/home/user/MARTA/LABO/PULPO/PULPO_synthetic/results/SVs/Patients/"

output_folder <- "/ProcessData/"
outputfolder <- output_folder
createbedpe <- function(file_path, outputfolder, outputdirbaseSVs) {
  base_directory <- file_path
  patients <- list.files(base_directory)
  for (i in patients) {
    intermediarydir <- paste0(base_directory, "/", i, "/OGMdata/")
    print(intermediarydir)
    
    # Listar todos los archivos .smap
    files <- list.files(intermediarydir, pattern = "\\.smap$", full.names = FALSE)
    
    # Priorizar los archivos que contengan "Rare_Variant_Analysis"
    preferred_files <- files[grepl("Rare_Variant_Analysis", files)]
    
    # Elegir archivo apropiado
    if (length(preferred_files) > 0) {
      smap_file <- preferred_files[1]
    } else if (length(files) > 0) {
      smap_file <- files[1]
      message("Warning: File with “Rare_Variant_Analysis” was not found for the patient ", i, ". Using file:", smap_file)
    } else {
      message("No .smap file found for patient ", i, ". Skipping.")
      next
    }
    
    full_dir <- paste0(intermediarydir, smap_file)
    output_dir <- paste0(outputdirbaseSVs, "/", i, outputfolder)
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    output_file <- paste0(output_dir, i, ".bedpe")
    
    bedpe <- tryCatch(
      {
        read.table(full_dir, header = TRUE) #this was false before synthethic data
      },
      error = function(e) {
        message("The file cannot be read or is empty. Creating empty .bedpe file for the patient", i)
        bedpe_final_order <- data.frame()
        write.table(bedpe_final_order, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        return(NULL)
      }
    )
    
    if (is.null(bedpe) || nrow(bedpe) == 0) {
      message("Skipping patient ", i, " per empty file.")
      next
    }
    
    chrom1 <- bedpe$V3
    start1 <- round(as.numeric(bedpe$V7))
    end1 <- round(as.numeric(bedpe$V7)) + 10000
    chrom2 <- bedpe$V4
    start2 <- round(as.numeric(bedpe$V8))
    end2 <- round(as.numeric(bedpe$V8)) + 10000
    name <- i
    score <- bedpe$V10
    strand <- bedpe$V24
    
    bedpefinal <- data.frame(
      chrom1 = chrom1, start1 = start1, end1 = end1,
      chrom2 = chrom2, start2 = start2, end2 = end2,
      name = name, score = score, strand = strand
    )
    
    bedpefinal <- bedpefinal %>%
      separate(strand, into = c("strand1", "strand2"), sep = "/") %>%
      mutate(
        strand1 = ifelse(is.na(strand1), ".", strand1),
        strand2 = ifelse(is.na(strand2), ".", strand2)
      )
    
    bedpe_final_order <- bedpefinal[order(
      bedpefinal$chrom1, bedpefinal$start1,
      bedpefinal$chrom2, bedpefinal$start2
    ), ]
    
    bedpe_final_order <- unique(bedpe_final_order)
    
    write.table(bedpe_final_order, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}

createbedpe(file_path = bionanodirectory, outputfolder = output_folder, outputdirbaseSVs = outputdirbaseSVs)
