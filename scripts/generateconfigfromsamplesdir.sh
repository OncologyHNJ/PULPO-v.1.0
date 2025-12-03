#!/bin/bash
OUT_DIR="/home/user/MARTA/LABO/PULPO/PULPO_OpenBTA/processed/CNVs/cnvkit_per_sample"
SAMPLES_TSV="/home/user/MARTA/LABO/PULPO/PULPO_OpenBTA/config/samples.tsv"
GENOME="hg38"

mkdir -p "$(dirname "$SAMPLES_TSV")"
echo -e "sample\tsource\tcnv_caller\tcnv_path\tgenome_build\tpurity\tploidy" > "$SAMPLES_TSV"
for f in "$OUT_DIR"/*.cns; do
  s=$(basename "$f" .cns)
  echo -e "${s}\twes\tcnvkit\t${f}\t${GENOME}\tNA\t2" >> "$SAMPLES_TSV"
done
echo ">> Escrito: $SAMPLES_TSV"
