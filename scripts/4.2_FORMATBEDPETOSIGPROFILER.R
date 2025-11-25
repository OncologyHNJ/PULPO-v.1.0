#######CREATED BY:MARTA PORTASANY########
#########################################
args <- commandArgs(trailingOnly = TRUE)
bedpedirectory <- args[1]
patientid <- args[2]
inputdata <- args[3]
outputdirectoryfile <- args[4]
#########################################
# bedpedirectory <-"/home/user/MARTA/LABO/PULPO/PULPO_synthetic/results/SVs/Patients/Syntheticpatient1/SigProfiler/data/SigProfilerSVdf.bedpe"
# patientid <- "Syntheticpatient1"
# inputdata <- "/home/user/MARTA/LABO/PULPO/PULPO_synthetic/results/SVs/Patients/Syntheticpatient1/SigProfiler/data/"
# # outputdirectory <- "/home/user/MARTA/PULPO_RevisadoBionano/results/Patients/Patient-109/SigProfiler/results/"
# outputdirectory <- "/home/user/MARTA/LABO/PULPO/PULPO_synthetic/results/SVs/Patients/Syntheticpatient1/SigProfiler/results/"
# #outputdirectoryfile <- "/home/user/MARTA/PULPO_RevisadoBionano/results/Patients/Patient-109/SigProfiler/data/SigProfilerSVdf.bedpe"
# #########################################
library(reticulate)
library(devtools)
library(dplyr)
library(magrittr)
#########################################

# Read the .bedpe file into a dataframe
#df <- read.table(bedpedirectory, header = FALSE)


df <- tryCatch(
  {
    read.table(bedpedirectory, header = FALSE)
  },
  error = function(e) {
    message("The file is empty or cannot be read. Creating an empty output file.")
    df_total <- data.frame()  # Crear un data frame vacío
    write.table(df_total, file = outputdirectoryfile, sep = "\t", row.names = FALSE, quote = FALSE)
    return(NULL)  # Devolver NULL para evitar más procesamiento
  }
)

# Si df es NULL, termina el script aquí
if (is.null(df) || nrow(df) == 0) {
  message("Skipping processing due to empty file.")
  quit(save = "no")  # Finaliza la ejecución del script
}


end1 <- df[3] - 10000
end2 <- df[6] - 10000

df1 <- cbind(df[c(1,2)], end1, df[c(4,5)], end2) 
df3 <- patientid
df4 <- df[8]

df_total <- cbind(df1,df3,df4)
colnames(df_total) <- c("chrom1","start1","end1","chrom2","start2","end2","sample","svclass")
df_total <- df_total %>% filter(svclass != "inversion_partial")


df_total <- df_total %>%
  mutate(svclass = case_when(
    svclass == "translocation_interchr" ~ "translocation",
    svclass == "translocation_intrachr" ~ "translocation",
    svclass == "duplication" ~ "tandem-duplication",
    svclass == "duplication_inverted" ~ "tandem-duplication",
    svclass == "duplication_split" ~ "tandem-duplication",
    svclass == "inversion_paired" ~ "inversion",
    svclass == "inversion_partial" ~ "inversion",
    

    TRUE ~ svclass  # Hold values that do not match any conditions
  ))

df_total <- subset(df_total, svclass != "insertion")
df_total <- subset(df_total, svclass != "insertion_nbase")


# Create the directory if it does not exist
if (!dir.exists(inputdata)) {
  dir.create(inputdata, recursive = TRUE)
}

write.table(df_total, file = outputdirectoryfile, sep = "\t", row.names = FALSE, quote = FALSE)
