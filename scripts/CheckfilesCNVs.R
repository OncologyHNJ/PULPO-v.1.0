############################################################################################################
args <- commandArgs(trailingOnly = TRUE)
patient <- args[1]
inputcnv <- args[2]
error_log <- file.path(inputcnv, "check_cnvs_error.txt")  # Error file
############################################################################################################
# Clear previous error file if exists
if (file.exists(error_log)) {
  file.remove(error_log)
}

# Search for CNV files with .csv or .txt extensions
csv_files <- list.files(inputcnv, pattern = "CNV.*\\.csv$", full.names = TRUE)
txt_files <- list.files(inputcnv, pattern = "CNV.*\\.txt$", full.names = TRUE)

# If there are no valid files, log error and exit
if (length(csv_files) == 0 && length(txt_files) == 0) {
  message <- paste0("ERROR: No valid CNV files (.csv or .txt) found for ", patient, ". Exiting.")
  write(message, file = error_log, append = TRUE)
  stop(message)
}

# Function to validate CNV files
check_cnv_file <- function(file_path, extension) {
  all_lines <- readLines(file_path)
  variant_lines <- grep("^[^#]", all_lines, value = TRUE)
  
  if (length(variant_lines) == 0) {
    message <- paste0("ERROR: The CNV ", extension, " file for ", patient, " is empty.")
    write(message, file = error_log, append = TRUE)
    stop(message)
  }
}

# Review CNV files (CSV and TXT)
if (length(csv_files) > 0) {
  check_cnv_file(csv_files, "CSV")
}
if (length(txt_files) > 0) {
  check_cnv_file(txt_files, "TXT")
}

