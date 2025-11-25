#!/usr/bin/env Rscript
# ==================================================
# Script: spiralize_signatures.R
# Description: Generates a spiral visualization of CNV or SV
#              signature contributions per sample.
#              Output paths are provided by Snakemake.
# Author: Marta Portasany-Rodríguez
# Pipeline: PULPO 🐙
# ==================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: spiralize_signatures.R <activities_file> <datatype: SV|CNV> <out_svg> [out_pdf]")
}

sig_file  <- args[1]
data_type <- toupper(args[2])
out_svg   <- args[3]
out_pdf   <- ifelse(length(args) >= 4,
                    args[4],
                    sub("\\.svg$", ".pdf", out_svg))

if (!data_type %in% c("SV", "CNV")) {
  stop("Data type must be 'SV' or 'CNV'.")
}

if (!file.exists(sig_file)) {
  stop("Input activities file not found: ", sig_file)
}

suppressPackageStartupMessages({
  library(spiralize)
  library(data.table)
  library(grid)
  library(grDevices)
})

message("📂 Reading ", data_type, " signature activities file: ", sig_file)
message("📝 Output SVG: ", out_svg)
message("📝 Output PDF: ", out_pdf)

# ── Load and normalize input data ─────────────────────────
signatures <- fread(sig_file)

if (!"Samples" %in% colnames(signatures)) {
  stop("Input file must contain a 'Samples' column.")
}

sample_names <- signatures$Samples
mat <- as.matrix(signatures[, -which(colnames(signatures) == "Samples"), with = FALSE])

if (!is.numeric(mat)) {
  mode(mat) <- "numeric"
}

# Handle rows with zero sum to avoid NaN
rs <- rowSums(mat, na.rm = TRUE)
rs[rs == 0] <- NA
mat_norm <- mat / rs
rownames(mat_norm) <- sample_names

# Remove rows that are all NA (originally all zeros)
keep_rows <- rowSums(!is.na(mat_norm)) > 0
mat_norm <- mat_norm[keep_rows, , drop = FALSE]

if (nrow(mat_norm) == 0) {
  stop("All rows have zero activity; nothing to plot.")
}

n_signatures <- ncol(mat_norm)
sig_names <- colnames(mat_norm)

if (n_signatures == 0) {
  stop("No signature columns found in activities file.")
}

# ── Color palette ─────────────────────────────────────────
cols <- setNames(hcl.colors(n_signatures, "Dark 3"), sig_names)

# ── Legend helper ─────────────────────────────────────────
draw_vertical_legend <- function(cols, labels, x = 0.98, y = 0.95) {
  n <- length(cols)

  # Altura de cada fila de leyenda
  row_h    <- unit(1.1, "lines")  # espacio por entrada
  box_size <- unit(0.9, "lines")  # tamaño del cuadrado de color

  # Viewport anclado arriba a la derecha
  vp <- viewport(
    x = x,
    y = y,
    width  = unit(0.22, "npc"),
    height = row_h * n,
    just   = c("right", "top")
  )
  pushViewport(vp)

  for (i in seq_along(cols)) {
    # Centro vertical de la fila i (1 = arriba)
    y_center <- unit(1, "npc") - (i - 0.5) * row_h

    # Cuadrado de color
    grid.rect(
      x = unit(0, "npc"),
      y = y_center,
      width  = box_size,
      height = box_size,
      just   = c("left", "center"),
      gp = gpar(fill = cols[i], col = NA)
    )

    # Texto alineado al centro vertical del cuadrado
    grid.text(
      labels[i],
      x = box_size + unit(0.4, "lines"),
      y = y_center,
      just = c("left", "center"),
      gp = gpar(cex = 0.9)
    )
  }

  popViewport()
}


# ── Spiral plot helper ────────────────────────────────────
plot_spiral <- function(title) {
  n_samples <- nrow(mat_norm)

  spiral_initialize(xlim = c(0.5, n_samples + 0.5))
  #spiral_track(height = 0.9)
  spiral_track(
  height    = 0.9
)

  # spiral_bars: matrix of values, one column per signature
  spiral_bars(
    value = mat_norm,
    pos = 1:n_samples,
    gp = gpar(fill = cols[colnames(mat_norm)], col = NA),
  #  bar_width = 0.8
    bar_width = 0.9
  )

  grid.text(title,
            x = 0.5,
            y = 0.97,
            gp = gpar(fontsize = 14, fontface = "bold"))

  draw_vertical_legend(cols, names(cols))
}

plot_title <- paste("Distribution of", data_type, "signatures per sample")

# ── SVG ───────────────────────────────────────────────────
svg(out_svg, width = 10, height = 9)
plot_spiral(plot_title)
dev.off()

# ── PDF ───────────────────────────────────────────────────
pdf(out_pdf, width = 10, height = 9)
plot_spiral(plot_title)
dev.off()

message("✅ Spiral plots generated successfully.") 

