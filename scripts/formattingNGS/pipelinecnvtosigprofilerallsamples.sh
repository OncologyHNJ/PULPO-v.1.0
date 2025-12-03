OUTROOT=/home/user/MARTA/LABO/PULPO/PULPO_OpenPBTA
SAMPLES=/home/user/MARTA/LABO/PULPO/PULPO_OpenPBTA/config/samplesall.tsv
LOHDIR=/home/user/refs/loh_beds   # opcional

# sample<TAB>anonymised<TAB>cnv_path
tail -n +2 "$SAMPLES" | awk -F'\t' -v OFS='\t' '{print $1,$2,$5}' | \
while IFS=$'\t' read -r SAMPLE ANON CNS; do
  [ -s "$CNS" ] || { echo "[WARN] no CNS para $SAMPLE -> $CNS"; continue; }
  OUT="$OUTROOT/results/CNVs/Patients/$ANON/SigProfiler/data/SigProfilerCNVdf.tsv"
  mkdir -p "$(dirname "$OUT")"
  LOH_ARG=""
  [ -n "$LOHDIR" ] && [ -s "$LOHDIR/$ANON.loh.bed" ] && LOH_ARG="--loh-bed $LOHDIR/$ANON.loh.bed"

  python /home/user/MARTA/LABO/PULPO/PULPO_OpenPBTA/processed/cnv_to_sigprofiler.py \
    --input "$CNS" \
    --input-type cnvkit \
    --sample "$ANON" \
    --genome-build hg38 \
    $LOH_ARG \
    --output "$OUT"
  echo "[OK] $SAMPLE -> $OUT"
done

