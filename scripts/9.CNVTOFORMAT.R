# ==================================================
# Script: 9.CNVTOFORMAT.R
# Description: Descompression of OGM raw data as a first step in the PULPO pipeline
# Author: Marta Portasany
# Created on: 2025-02-27
# Last modified: 2025-02-27
# Pipeline: PULPO
# Dependencies: dplyr
# ==================================================

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
inputdata <- args[1]
output <- args[2]
patient <- args[3]

pattern <- "CNV"
file_path <- list.files(path = inputdata, pattern = pattern, full.names = TRUE)

# Check if the file has only headers (i.e. no data)
if (length(file_path) == 0) {
  # Create an empty file with only the headers
  dfcnv <- data.frame(sample = character(), chromosome = character(), chromosome_start = numeric(), chromosome_end = numeric(), copy_number = numeric(), mutation_type = character())
  if (!dir.exists(dirname(output))) {
    dir.create(dirname(output), recursive = TRUE)
  }
  write.table(dfcnv, file = output, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  ###FORMATITINGHEADERS###
  lines <- readLines(file_path)
  header_lines <- grep("^#", lines, value = TRUE)
  last_header_line <- tail(header_lines, 1)
  column_names <- strsplit(last_header_line, split = "\t")[[1]]
  column_names <- gsub("^#", "", column_names) # REMOVE '#'
  variant_lines <- grep("^#", lines, value = TRUE, invert = TRUE)
  
  if (length(variant_lines) == 0){
    dfcnv <- data.frame(sample = character(), chromosome = character(), chromosome_start = numeric(), chromosome_end = numeric(), copy_number = numeric(), mutation_type = character())
    write.table(dfcnv, file = output, sep = "\t", row.names = FALSE, col.names = TRUE)
  } else {
    
  txt <- list.files(path = inputdata, pattern = "txt", full.names = TRUE)
  csv <- list.files(path = inputdata, pattern = ".csv", full.names = TRUE)
  if (length(txt) == 1) {
    filecnv <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    colnames(filecnv) <- column_names 
    chromosome <- filecnv$Chromosome
    chromosome_start <- round(filecnv$Start)
    chromosome_end <- round(filecnv$End)
    copy_number <- filecnv$CopyNumber
    mutation_type <- gsub("(loss|gain)_masked", "\\1", filecnv$Type)
    sample <- patient
    dfcnv <- data.frame(sample, chromosome, chromosome_start, chromosome_end, copy_number, mutation_type)
    dfcnv <- dfcnv %>%
      filter(!(copy_number== 2 & mutation_type == "loss")) %>% 
      mutate(copy_number = ifelse(copy_number == 0,1, copy_number))
   # print(dfcnv)  
    if (!dir.exists(dirname(output))) {
      dir.create(dirname(output), recursive = TRUE)
    }
    write.table(dfcnv, file = output, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  else if (length(csv) == 1){
    filecnv <- read.csv(file_path,header = TRUE)
    chromosome <- filecnv$Chromosome
    chromosome_start <- round(filecnv$Start)
    chromosome_end <- round(filecnv$End)
    copy_number <- filecnv$CopyNumber
    mutation_type <- gsub("(loss|gain)_masked", "\\1", filecnv$Type)
    sample <- patient
   
    dfcnv <- data.frame(sample,chromosome,chromosome_start, chromosome_end, copy_number, mutation_type)
   
    dfcnv <- dfcnv %>%
      filter(!(copy_number== 2 & mutation_type == "loss")) %>% 
      mutate(copy_number = ifelse(copy_number == 0,1, copy_number))

    if (!dir.exists(dirname(output))) {
      dir.create(dirname(output), recursive = TRUE)
    }
    write.table(dfcnv, file = output, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  else {
    print("The extension of your CNV file is not correct")
  }
 
} }

