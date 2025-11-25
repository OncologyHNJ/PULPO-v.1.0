#############################################
## Script R: Descarga CNVs/SVs 

## Usando GenomicDataCommons
## Marta Portasany
## Fecha: 19/09/2025
#############################################

# Instalar paquetes si no están presentes
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("GenomicDataCommons", quietly = TRUE)) BiocManager::install("GenomicDataCommons")
if(!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if(!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if(!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
if(!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

library(GenomicDataCommons)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# 1. Comprobar conexión a GDC
stopifnot(GenomicDataCommons::status()$status == "OK")

###Preludio
projects() %>% GenomicDataCommons::filter(~ project_id == "TARGET-ALL-P2") %>% results_all()

files() %>%
  GenomicDataCommons::filter(~ cases.project.project_id == "TARGET-ALL-P2") %>%
  facet(c("type", "data_category")) %>%
  aggregations()

# Campos que queremos
desired_fields <- c(
  "cases.project.project_id",
  default_fields("files"),
  grep_fields("files", "associated_entities"),
  "analysis.analysis_type",
  "analysis.workflow_type",
  "analysis.workflow_version"
)

# Filtrar por tipos CNV/SV existentes
qfiles <- files(fields = desired_fields) %>%
  GenomicDataCommons::filter(
    ~ type %in% c(
      "structural_variation",
      "copy_number_segment",
      "copy_number_estimate",
      "copy_number_auxiliary_file"
    ) &
      cases.project.project_id == "TARGET-ALL-P2"
  )

# Obtener todos los resultados
res.cnvsv <- results_all(qfiles)

# Convertir a data frame para revisar
res_df <- as.data.frame(res.cnvsv)
dim(res_df)
head(res_df)

length(res.cnvsv$file_id)
# o
nrow(res_df)

# Ahora sí podemos usar count() de dplyr
library(dplyr)
res_df %>% count(type) 
# Convertir a data.frame
res_df <- as.data.frame(res.cnvsv)

# Ahora sí podemos usar count() de dplyr
library(dplyr)
res_df %>% count(type) 
qfiles %>% count()

# 3. Obtener resultados completos y crear ID map
res.cnvsv <- results_all(qfiles)

# Eliminar archivos con múltiples IDs asociados
idx <- sapply(res.cnvsv$associated_entities, nrow) %>% grep(2, .)
ID.map <- res.cnvsv[!grepl("list|data.frame", sapply(res.cnvsv, class))] %>%
  as.data.frame() %>%
  slice(-idx) %>%
  mutate(project.project_id = unlist(res.cnvsv$cases[-idx])) %>%
  bind_cols(., bind_rows(res.cnvsv$associated_entities[-idx]))

head(ID.map)
dim(ID.map)

# 4. Crear manifest de descarga
manifest_df <- manifest(qfiles)
write.table(manifest_df, "TARGET-ALL-P2_CNV_SV_manifest.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

# Configurar cache/download folder
dir.create("CNV_SV_Data", showWarnings = FALSE)
gdc_set_cache("CNV_SV_Data/")

# Descargar los archivos
fnames <- gdcdata(manifest_df$id, progress = TRUE)


# 5. Configurar directorio de descarga
dir.create("CNV_SV_Data", showWarnings = FALSE)
gdc_set_cache(directory = "CNV_SV_Data/")

# 6. Descargar archivos
fnames <- gdcdata(manifest_df$id, progress = TRUE, access_method = "api", use_cached = FALSE)
head(fnames)

# 7. Reorganizar archivos CNV/SV
files_to_cat <- list.files("CNV_SV_Data", recursive = TRUE, full.names = TRUE)
# Dependiendo del tipo de archivo, aquí solo movemos todos a la carpeta base
for(f in files_to_cat) {
  file.rename(f, file.path("CNV_SV_Data", basename(f)))
}

# 8. Descargar datos clínicos del proyecto TARGET-ALL-P2
case_ids <- cases() %>%
  filter(~ project.project_id == "TARGET-ALL-P2") %>%
  ids()

clin_res <- gdc_clinical(case_ids)
sapply(clin_res, dim)

# Combinar todos los datos clínicos
full_clin <- clin_res$main %>%
  left_join(clin_res$demographic, by = "case_id") %>%
  left_join(clin_res$exposures, by = "case_id") %>%
  left_join(clin_res$diagnoses, by = "case_id")

head(full_clin)
dim(full_clin)

# 9. Guardar resultados
write.csv(full_clin, "TARGET_ALL_P2_Clinical_Data.csv", row.names = FALSE)
cat("Descarga y organización de CNVs/SVs y datos clínicos completadas.\n")
