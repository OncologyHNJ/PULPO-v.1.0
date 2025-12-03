#!/usr/bin/env python3
import argparse, os, sys, gzip
import pandas as pd

def open_maybe_gz(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')

def merge_intervals(df):
    # df: columns chr,start,end (enteros)
    df = df.sort_values(["chr","start","end"]).reset_index(drop=True)
    out = []
    cur_chr, cur_s, cur_e = None, None, None
    for _,r in df.iterrows():
        c, s, e = r["chr"], int(r["start"]), int(r["end"])
        if cur_chr is None:
            cur_chr, cur_s, cur_e = c, s, e
            continue
        if c == cur_chr and s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            out.append((cur_chr, cur_s, cur_e))
            cur_chr, cur_s, cur_e = c, s, e
    if cur_chr is not None:
        out.append((cur_chr, cur_s, cur_e))
    return pd.DataFrame(out, columns=["chr","start","end"])

def main():
    ap = argparse.ArgumentParser(
        description="Construye BEDs de LOH por muestra a partir de pbta-all-gene-loh.tsv.gz y un BED de genes."
    )
    ap.add_argument("--pbta-gene-loh", required=True,
                    help="Ruta a pbta-all-gene-loh.tsv.gz (columns: sample_id, gene_symbol, loh)")
    ap.add_argument("--genes-bed", required=True,
                    help="BED de genes (chr, start, end, gene_symbol ...)")
    ap.add_argument("--outdir", required=True, help="Carpeta de salida para SAMPLE.loh.bed")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # lee gene LOH
    with open_maybe_gz(args.pbta_gene_loh) as fh:
        loh = pd.read_csv(fh, sep="\t", dtype=str)

    # normaliza
    cols = {c.lower(): c for c in loh.columns}
    for need in ["sample_id","gene_symbol","loh"]:
        if need not in cols:
            raise SystemExit(f"Falta columna {need} en {args.pbta_gene_loh}")
    loh = loh.rename(columns={
        cols["sample_id"]: "sample_id",
        cols["gene_symbol"]: "gene_symbol",
        cols["loh"]: "loh"
    })
    loh = loh[loh["loh"].str.upper().isin(["TRUE","T","1","YES"])].copy()

    # lee genes BED
    genes = pd.read_csv(args.genes_bed, sep="\t", header=None,
                        usecols=[0,1,2,3],
                        names=["chr","start","end","gene_symbol"],
                        dtype={"chr":str})
    genes["chr"] = genes["chr"].astype(str).str.replace("^chr","", regex=True)

    # join por gene_symbol
    m = loh.merge(genes, how="inner", on="gene_symbol")
    if m.empty:
        sys.stderr.write("[WARN] No hay solape entre pbta LOH y lista de genes.\n")

    # para cada muestra, merge de intervalos y emitir BED
    for sid, sub in m.groupby("sample_id"):
        sub = sub[["chr","start","end"]].copy()
        sub["start"] = pd.to_numeric(sub["start"], errors="coerce").astype("Int64")
        sub["end"]   = pd.to_numeric(sub["end"], errors="coerce").astype("Int64")
        sub = sub.dropna()
        sub = sub.sort_values(["chr","start","end"])
        # por cromosoma
        merged_all = []
        for c, g in sub.groupby("chr"):
            merged = merge_intervals(g.assign(chr=c))
            merged_all.append(merged)
        if not merged_all:
            continue
        bed = pd.concat(merged_all, ignore_index=True)
        outpath = os.path.join(args.outdir, f"{sid}.loh.bed")
        bed.to_csv(outpath, sep="\t", header=False, index=False)
        sys.stderr.write(f"[OK] Escrito {outpath} ({len(bed)} tramos)\n")

if __name__ == "__main__":
    main()

