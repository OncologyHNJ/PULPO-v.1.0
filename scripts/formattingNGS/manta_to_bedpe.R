#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("
manta_to_sigprofiler_bedpe.R
--------------------------------
Convierte pbta-sv-manta.tsv(.gz) (Manta+AnnotSV) a BEDPE mínimo por paciente
para SigProfiler (svclass: deletion/duplication/tandem-duplication/inversion/translocation).

USO:
  Rscript manta_to_sigprofiler_bedpe.R \\
    --in pbta-sv-manta.tsv.gz \\
    --map /ruta/samplesall.tsv \\
    --outdir /.../results/SVs/Patients \\
    [--pass-only] [--dup-as-tandem] [--consolidated /ruta/all.bedpe]

FLAGS:
  --in             Fichero TSV/TSV.GZ de OpenPBTA (Manta anota con AnnotSV)
  --map            samplesall.tsv con columnas: sample↔anonymised/anonymized
  --outdir         Carpeta base: .../Patients/<anon>/SigProfiler/data/SigProfilerSVdf.bedpe
  --pass-only      (opcional) filtra SVs con FILTER=PASS
  --dup-as-tandem  (opcional) mapea cualquier DUP como tandem-duplication
  --consolidated   (opcional) escribe también un BEDPE con toda la cohorte

Ejemplo:
  Rscript manta_to_sigprofiler_bedpe.R \\
    --in pbta-sv-manta.tsv.gz \\
    --map config/samplesall.tsv \\
    --outdir /home/user/MARTA/LABO/PULPO/PULPO_OpenPBTA/results/SVs/Patients \\
    --pass-only --consolidated /tmp/openpbta_all_sv.bedpe
\n")
  quit(status = 1)
}

# --- parse args sencillito ---
get_flag <- function(key) any(grepl(paste0("^", key, "$"), args))
get_opt  <- function(key, default = NA_character_) {
  hit <- which(args == key)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[hit + 1]
}

if (get_flag("--help") || length(args) == 0) usage()

infile  <- get_opt("--in")
mapfile <- get_opt("--map")
outdir  <- get_opt("--outdir")
consol  <- get_opt("--consolidated", NA_character_)
pass_only     <- get_flag("--pass-only")
dup_as_tandem <- get_flag("--dup-as-tandem")

if (is.na(infile) || is.na(mapfile) || is.na(outdir)) usage()

message("[INFO] infile:      ", infile)
message("[INFO] mapfile:     ", mapfile)
message("[INFO] outdir base: ", outdir)
if (!is.na(consol)) message("[INFO] consolidated: ", consol)
if (pass_only) message("[INFO] pass-only activado")
if (dup_as_tandem) message("[INFO] dup-as-tandem activado")

# --- helpers ---
to_num_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x %in% c("X")]     <- "23"
  x[x %in% c("Y")]     <- "24"
  x[x %in% c("M","MT")]<- "25"
  x
}

# Soporta los 4 patrones de BND; extrae primero chr:pos entre corchetes
parse_alt_mate <- function(alt) {
  # ejemplo válidos: N[chr10:126891954[, ]chr22:41050406]N, T]chr3:176966137], [chr7:133415006[A
  m <- str_match(alt, "(?:\\[|\\])?\\s*chr?([0-9XYMT]+)\\s*:(\\d+)\\s*(?:\\[|\\])?")
  if (is.na(m[1,1])) return(list(chr = NA_character_, pos = NA_integer_))
  list(chr = to_num_chr(m[1,2]), pos = as.integer(m[1,3]))
}

get_info_field <- function(info, key) {
  p <- paste0(";", key, "=([^;]+)")
  m <- str_match(paste0(";", info), p)
  fifelse(is.na(m[,2]), NA_character_, m[,2])
}

safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# --- read input ---
dt <- fread(infile)
setnames(dt, make.names(names(dt)))  # nombres seguros tipo SV.type, AnnotSV.type, etc.

if (pass_only && "FILTER" %in% names(dt)) {
  dt <- dt[FILTER == "PASS"]
}

# trabajamos solo AnnotSV.type == "full" si existe la columna
if ("AnnotSV.type" %in% names(dt)) {
  dt <- dt[is.na(AnnotSV.type) | AnnotSV.type == "full"]
}

# columna de muestra en el TSV de OpenPBTA
sample_col <- "Kids.First.Biospecimen.ID.Tumor"
stopifnot(sample_col %in% names(dt))

# normaliza cromosoma principal y posiciones
dt[, chrom1 := to_num_chr(SV.chrom)]
dt[, pos1   := as.integer(SV.start)]
dt[, pos2raw:= as.integer(SV.end)]

# mapea clases
raw_type <- dt$SV.type
svmap <- c(
  "DEL"         = "deletion",
  "INV"         = "inversion",
  "DUP"         = "tandem-duplication",
  "DUP:TANDEM"  = "tandem-duplication",
  "BND"         = "translocation",
  "TRA"         = "translocation"
)
dt[, svclass := svmap[raw_type]]
dt <- dt[!is.na(svclass)]

# --- intrachrom: DEL/DUP/INV ---
intra <- dt[SV.type %in% c("DEL","DUP","DUP:TANDEM","INV")]
intra_out <- intra[, .(
  chrom1 = chrom1,
  start1 = pos1, end1 = pos1,
  chrom2 = chrom1,
  start2 = pos2raw, end2 = pos2raw,
  sample = get(sample_col),
  svclass = svclass
)]

# --- BND / TRA: parse ALT ---
bnd <- dt[SV.type %in% c("BND","TRA")]
if (nrow(bnd)) {
  bnd[, event  := get_info_field(INFO, "EVENT")]
  tmp <- bnd[, {
    m <- parse_alt_mate(ALT)
    .(sample = get(sample_col),
      event2 = fifelse(!is.na(event), event, ID),
      chrom1 = chrom1,
      pos1   = pos1,
      chrom2 = m$chr,
      pos2   = m$pos)
  }, by = ID]
  
  # avisos por mates no parseables
  n_na <- sum(is.na(tmp$chrom2) | is.na(tmp$pos2))
  if (n_na > 0) warning(sprintf("[WARN] %d BND sin mate parseable por ALT; se omiten.", n_na))
  tmp <- tmp[!is.na(chrom2) & !is.na(pos2)]
  
  # dedupe por evento (si dos mitades, quédate con 1) + clave no dirigida final
  setorder(tmp, sample, event2, chrom1, pos1, chrom2, pos2, na.last = TRUE)
  tmpu <- tmp[!duplicated(paste(sample, event2))]
  
  bnd_out <- tmpu[, .(
    chrom1 = chrom1, start1 = pos1, end1 = pos1,
    chrom2 = chrom2, start2 = pos2, end2 = pos2,
    sample = sample,
    svclass = "translocation"
  )]
} else {
  bnd_out <- data.table()
}

# --- unir ---
out <- rbind(intra_out, bnd_out, fill = TRUE)
out <- out[complete.cases(out[, .(chrom1,start1,chrom2,start2)])]

# --- segunda capa de dedupe no dirigida para TRA (por si falta EVENT) ---
is_tra <- out$svclass == "translocation"
if (any(is_tra)) {
  tra <- out[is_tra]
  swap <- (tra$chrom1 > tra$chrom2) | (tra$chrom1 == tra$chrom2 & tra$start1 > tra$start2)
  a_chr <- ifelse(swap, tra$chrom2, tra$chrom1)
  a_pos <- ifelse(swap, tra$start2, tra$start1)
  b_chr <- ifelse(swap, tra$chrom1, tra$chrom2)
  b_pos <- ifelse(swap, tra$start1, tra$start2)
  tra[, key := paste(sample, a_chr, a_pos, b_chr, b_pos, sep="|")]
  tra <- tra[!duplicated(key)][, key := NULL]
  out <- rbind(out[!is_tra], tra, fill = TRUE)
}

# --- mapping a anonymized ---
map <- fread(mapfile)
cn <- tolower(names(map))

cand_sample_cols <- c("sample","biospecimen_id","kids.first.biospecimen.id.tumor",
                      "kids.first.biospecimen.id.normal","tumor_sample","tumor_id")
sample_col_map <- names(map)[match(cand_sample_cols, cn, nomatch = 0)]
if (length(sample_col_map) == 0) {
  stop("No encuentro columna de sample en samplesall.tsv. Esperaba una de: ",
       paste(cand_sample_cols, collapse=", "))
}
sample_col_map <- sample_col_map[1]

cand_anon_cols <- c("anonymised","anonymized","anonymised_id","anonymized_id",
                    "patient_anon","patient_anonymised","patient_anonymized")
anon_col_map <- names(map)[match(cand_anon_cols, cn, nomatch = 0)]
if (length(anon_col_map) == 0) {
  stop("No encuentro columna anonymized en samplesall.tsv. Esperaba una de: ",
       paste(cand_anon_cols, collapse=", "))
}
anon_col_map <- anon_col_map[1]

map_vec <- setNames(as.character(map[[anon_col_map]]), as.character(map[[sample_col_map]]))
out[, patient := map_vec[sample]]
n_na_pat <- sum(is.na(out$patient))
if (n_na_pat > 0) {
  warning(sprintf("[WARN] %d muestras sin mapping a anonymized; usaré el ID de sample.", n_na_pat))
  out[is.na(patient), patient := sample]
}
out[, sample := patient]

# --- escribir por paciente ---
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cols <- c("chrom1","start1","end1","chrom2","start2","end2","sample","svclass")

for (p in sort(unique(out$patient))) {
  sub <- out[patient == p, ..cols]
  if (nrow(sub) == 0) next
  pdir <- file.path(outdir, safe_name(p), "SigProfiler", "data")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  of <- file.path(pdir, "SigProfilerSVdf.bedpe")
  fwrite(sub, of, sep = "\t", quote = FALSE)
  message(sprintf("[OK] %s: %d SVs → %s", p, nrow(sub), of))
}

# --- opcional consolidado ---
if (!is.na(consol)) {
  setcolorder(out, cols)
  fwrite(out[, ..cols], consol, sep = "\t", quote = FALSE)
  message(sprintf("[OK] Consolidated BEDPE: %s (%d filas)", consol, nrow(out)))
}
