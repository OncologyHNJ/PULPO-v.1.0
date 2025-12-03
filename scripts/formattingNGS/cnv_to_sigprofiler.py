#!/usr/bin/env python3
import argparse, os, gzip, sys
import pandas as pd
from collections import defaultdict

def open_maybe_gz(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')

def load_loh_bed(loh_bed, sample):
    """
    Carga y combina todos los intervalos LOH para una muestra.
    - loh_bed puede ser:
      * un fichero concreto (*.bed) con intervalos para UNA muestra
      * una carpeta que contenga SAMPLE.loh.bed por muestra
    Devuelve: dict chr -> lista de (start,end) no solapados (merge simple).
    """
    if loh_bed is None:
        return {}

    if os.path.isdir(loh_bed):
        cand = os.path.join(loh_bed, f"{sample}.loh.bed")
        if not os.path.exists(cand):
            return {}
        bed_path = cand
    else:
        bed_path = loh_bed

    try:
        df = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0,1,2],
                         names=["chrom","start","end"], dtype={"chrom":str})
    except Exception as e:
        sys.stderr.write(f"[WARN] No pude leer LOH BED {bed_path}: {e}\n")
        return {}

    # normaliza cromosomas sin 'chr'
    df["chrom"] = df["chrom"].astype(str).str.replace("^chr","", regex=True)
    merged = defaultdict(list)
    for chrom, sub in df.groupby("chrom", sort=False):
        s = sub.sort_values(["start","end"]).reset_index(drop=True)
        out = []
        for _, row in s.iterrows():
            if not out:
                out.append([int(row.start), int(row.end)])
            else:
                if int(row.start) <= out[-1][1]:  # solapa
                    out[-1][1] = max(out[-1][1], int(row.end))
                else:
                    out.append([int(row.start), int(row.end)])
        merged[chrom] = [(a,b) for a,b in out]
    return merged

def overlaps_any(chrom, s, e, loh_intervals):
    if chrom not in loh_intervals:
        return False
    # búsqueda lineal (suficiente para cientos/ pocos miles de tramos)
    for a,b in loh_intervals[chrom]:
        if not (e <= a or s >= b):  # intersección estricta
            return True
    return False

def cn_from_cnvkit(row, ploidy=2):
    # Si trae CN directo (copy.num / cn), úsalo; si no, convierte desde log2:
    if "copy.num" in row and pd.notna(row["copy.num"]):
        try:
            return int(round(float(row["copy.num"])))
        except:
            pass
    if "cn" in row and pd.notna(row["cn"]):
        try:
            return int(round(float(row["cn"])))
        except:
            pass
    # CNVkit típico: seg.mean (log2) -> CN ≈ round(ploidy * 2**log2)
    for cand in ["seg.mean","log2","segmean"]:
        if cand in row and pd.notna(row[cand]):
            try:
                return int(round(ploidy * (2 ** float(row[cand]))))
            except:
                continue
    # Último recurso: 2
    return 2

def cn_from_consensus(row):
    # OpenPBTA consensus trae copy.num
    if "copy.num" in row and pd.notna(row["copy.num"]):
        try:
            return int(round(float(row["copy.num"])))
        except:
            pass
    # Si no, intenta seg.mean
    for cand in ["seg.mean","segmean"]:
        if cand in row and pd.notna(row[cand]):
            try:
                return int(round(2 * (2 ** float(row[cand]))))  # suposición ploidy=2
            except:
                continue
    return 2

def main():
    ap = argparse.ArgumentParser(
        description="Convierte segmentos CNV (CNVkit .cns o OpenPBTA consensus SEG) "
                    "a formato PCAWG para SigProfiler (opcionalmente anotando LOH en CN=2)."
    )
    ap.add_argument("--input", required=True,
                    help="Ruta a .cns (CNVkit) o pbta-cnv-consensus.seg(.gz)")
    ap.add_argument("--input-type", required=True, choices=["cnvkit","consensus"],
                    help="Tipo de entrada: 'cnvkit' o 'consensus'")
    ap.add_argument("--sample", required=True,
                    help="ID de muestra a escribir en 'sample'")
    ap.add_argument("--ploidy", type=float, default=2.0,
                    help="Ploidía a usar si hay que convertir log2 -> CN (default: 2)")
    ap.add_argument("--loh-bed", default=None,
                    help="(Opcional) BED de LOH (o carpeta con SAMPLE.loh.bed) para marcar CN=2 como 'copy neutral LOH'")
    ap.add_argument("--genome-build", default="hg38",
                    help="Solo para documentar en cabecera (no altera coordenadas)")
    ap.add_argument("--output", required=True,
                    help="TSV de salida: sample,chromosome,chromosome_start,chromosome_end,copy_number,mutation_type")
    args = ap.parse_args()

    # lee CNV
    with open_maybe_gz(args.input) as fh:
        df = pd.read_csv(fh, sep="\t", comment="#", dtype=str).fillna("")

    # normaliza nombres de columnas
    cols = {c.lower(): c for c in df.columns}

    def pick(*names):
        for n in names:
            if n in cols:
                return df[cols[n]]
        raise SystemExit(f"Falta alguna columna requerida ({'/'.join(names)}) en {args.input}; columnas={list(df.columns)}")

    if args.input_type == "cnvkit":
        chrom = pick("chrom","chr","chromosome").astype(str)
        start = pick("start","start.bp","loc.start","locstart")
        end   = pick("end","end.bp","loc.end","locend")
        df["_chrom"]=chrom
        df["_start"]=pd.to_numeric(start, errors="coerce").astype("Int64")
        df["_end"]  =pd.to_numeric(end,   errors="coerce").astype("Int64")
        # calcula CN por fila
        df["_cn"] = df.apply(lambda r: cn_from_cnvkit({k.lower(): r[k] for k in df.columns}, args.ploidy), axis=1)

    else:  # consensus
        chrom = pick("chrom","chr","chromosome").astype(str)
        start = pick("loc.start","start","chromosome_start")
        end   = pick("loc.end","end","chromosome_end")
        df["_chrom"]=chrom
        df["_start"]=pd.to_numeric(start, errors="coerce").astype("Int64")
        df["_end"]  =pd.to_numeric(end,   errors="coerce").astype("Int64")
        df["_cn"] = df.apply(lambda r: cn_from_consensus({k.lower(): r[k] for k in df.columns}), axis=1)

    # limpia chr
    df["_chrom"] = df["_chrom"].str.replace("^chr","", regex=True)

    # carga LOH (si hay)
    loh_map = load_loh_bed(args.loh_bed, args.sample) if args.loh_bed else {}

    # asigna mutation_type
    muts = []
    for c, s, e, cn in df[["_chrom","_start","_end","_cn"]].itertuples(index=False, name=None):
        if pd.isna(s) or pd.isna(e) or pd.isna(cn):
            muts.append("copy neutral")
            continue
        cn = int(cn)
        if cn <= 1:
            muts.append("loss")
        elif cn >= 3:
            muts.append("gain")
        else:  # cn == 2
            if loh_map and overlaps_any(str(c), int(s), int(e), loh_map):
                muts.append("copy neutral LOH")
            else:
                muts.append("copy neutral")
    df["_mutation_type"] = muts

    out = pd.DataFrame({
        "sample": args.sample,
        "chromosome": df["_chrom"].astype(str),
        "chromosome_start": df["_start"].astype("int64", errors="ignore"),
        "chromosome_end": df["_end"].astype("int64", errors="ignore"),
        "copy_number": df["_cn"].astype("int64", errors="ignore"),
        "mutation_type": df["_mutation_type"].astype(str)
    })
    out = out.dropna(subset=["chromosome_start","chromosome_end","copy_number"])

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    sys.stderr.write(f"[OK] Escrito {len(out)} filas en {args.output}\n")

if __name__ == "__main__":
    main()

