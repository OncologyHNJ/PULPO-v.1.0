### IMPORTS
import pandas as pd
import os
import yaml
from pathlib import Path

#configfile: "./config/configpulpoOGM.yaml"

# Load YAML
with open("./config/configpulpoOGM.yaml", "r") as fh:
    config = yaml.safe_load(fh)

# Helpers
def read_tsv(p):
    return pd.read_table(p, sep="\t", dtype=str).fillna("")

def ensure_col(df, col, value=None, default=""):
    if col not in df.columns:
        if value is not None and value in df.columns:
            df[col] = df[value]
        else:
            df[col] = default
    return df

# Read sample table
sample_path = config["input"]["samples"]
sample_table = read_tsv(sample_path)

# Optional anonymisation map (raw_id -> anonymised)
anon_map_path = config.get("input", {}).get("anonymisation_map", "")
anon_df = None
if anon_map_path and Path(anon_map_path).exists():
    anon_df = read_tsv(anon_map_path)
    needed = {"raw_id", "anonymised"}
    if not needed.issubset(set(anon_df.columns)):
        raise KeyError(f"anonymisation_map must contain columns: {needed}")

# Detect whether we truly need strict OGM columns this run
ft = config.get("file_types", {})
analysis_type = config.get("analysis", {}).get("analysis_type")

uses_ogm = (
    ((analysis_type in ["SVs", "SVs_and_CNVs"]) and (ft.get("SVs") in ["smap"])) or
    ((analysis_type in ["CNVs", "SVs_and_CNVs"]) and (ft.get("CNVs") in ["csv", "txt"]))
)


# Always require at least 'sample'
if "sample" not in sample_table.columns:
    raise KeyError("samples.tsv must contain at least a 'sample' column.")


# Build/validate anonymised (always present)
if "anonymised" not in sample_table.columns or (sample_table["anonymised"] == "").all():
    if anon_df is not None:
        m = anon_df.set_index("raw_id")["anonymised"]
        sample_table["anonymised"] = sample_table["sample"].map(m)
        # fallback to identity if any missing
        sample_table["anonymised"] = sample_table["anonymised"].fillna(sample_table["sample"])
    else:
        # No map provided -> identity mapping
        print("WARN: anonymisation_map not set; using sample as anonymised.")
        sample_table["anonymised"] = sample_table["sample"]

# Ensure other legacy columns exist (to keep OGM rules happy even si hoy solo WES)
sample_table = ensure_col(sample_table, "family", default="NA")
sample_table = ensure_col(sample_table, "bionano", value="sample")  # identity by default

# Strict check only if truly using OGM on this run
if uses_ogm:
    expected_columns = {'sample', 'family', 'bionano', 'anonymised'}
    if not expected_columns.issubset(sample_table.columns):
        raise KeyError(f"The sample file does not contain the necessary columns: {expected_columns}")

# Dictionaries
bionano_to_anonymised = sample_table.set_index('bionano')['anonymised'].to_dict()
anonymised_to_bionano = {v: k for k, v in bionano_to_anonymised.items()}

# Wildcards
wildcard_constraints:
    bionano = "|".join(sample_table['bionano'].unique()),
    anonymised = "|".join(sample_table['anonymised'].unique()),
    method = "cosmic|drews|tao"

# Working dirs
workdirectory = config['directories']['workdirectory']
logsdir = f"{workdirectory}/logs"
resultsdir = f"{workdirectory}/results"
os.makedirs(workdirectory, exist_ok=True)
os.makedirs(logsdir, exist_ok=True)
os.makedirs(resultsdir, exist_ok=True)




def validate_config(config):
    print("PULPO v1.0 is running. Please take care of the analysis!")
    valid_types = {"SVs", "CNVs", "SVs_and_CNVs"}

    analysis_config = config.get('analysis', {})
    if 'analysis_type' not in analysis_config:
        raise KeyError("The configuration file is missing the 'analysis_type' under 'analysis'.")
    if analysis_config['analysis_type'] not in valid_types:
        raise ValueError(f"Invalid analysis_type: {analysis_config['analysis_type']}. Must be one of {valid_types}.")

    if analysis_config['analysis_type'] == "SVs":
        print("You have selected 'SVs' as the analysis type. Proceeding with SVs analysis...")
    elif analysis_config['analysis_type'] == "CNVs":
        print("You have selected 'CNVs' as the analysis type. Proceeding with CNVs analysis...")
    else:
        print("You have selected 'SVs_and_CNVs' as the analysis type. Proceeding with both SVs and CNVs analysis...")

    run_cohort_analysis = analysis_config.get("run_cohort_analysis", False)
    zip_flag = bool(config.get("analysis",{}).get("zip",True))

    if run_cohort_analysis:
        print("Cohort analysis is enabled.")
    else:
        print("Running individual analysis only.")
    print("Configuration is valid! Proceeding with the pipeline...")

validate_config(config)


################################################################################################################################################################################
###########################################################################DEFINITIONSOFPULPOS'OPTIONS##########################################################################
################################################################################################################################################################################
pythonenvdir = config['directories']['pythonenvdir']
snakefiles_to_include = []
all_targets = []

in_CNVs = config.get("inputs", {}).get("CNVs", {})
in_SVs  = config.get("inputs", {}).get("SVs", {})
pp_CNVs = config.get("pipelines", {}).get("CNVs", {})
pp_SVs  = config.get("pipelines", {}).get("SVs", {})

cnv_source = in_CNVs.get("source", "ngs")      # "ogm" | "ngs"
cnv_format = in_CNVs.get("format", "cns")      # ogm: csv|txt ; ngs: cns|bed
sv_source  = in_SVs.get("source", "ngs")       # "ogm" | "ngs"
sv_format  = in_SVs.get("format", "bedpe")     # ogm: smap ; ngs: bedpe

cnv_methods = [m.lower() for m in pp_CNVs.get("methods", [])]   # p.ej., ["cosmic","drews","tao"]
sv_methods  = [m.lower() for m in pp_SVs.get("methods", [])]    # p.ej., ["cosmic"]

analysis_type = config.get("analysis", {}).get("analysis_type")
run_cohort = bool(config.get("analysis", {}).get("run_cohort_analysis",  True))
zip_flag = bool(config.get("analysis", {}).get("zip", True))   # 👈 NUEVA LÍNEA


ogm_for_svs  = (analysis_type in ["SVs","SVs_and_CNVs"]) and (sv_source == "ogm")  and (sv_format == "smap")
ogm_for_cnvs = (analysis_type in ["CNVs","SVs_and_CNVs"]) and (cnv_source == "ogm") and (cnv_format in ["csv","txt"])
uses_ogm_run = ogm_for_svs or ogm_for_cnvs

if uses_ogm_run:
    include: "rules/0_prepare_ogm_data.smk"
    ruleorder: prepare_ogm_data > decompressOGM

uses_ogm = uses_ogm_run

# ---------------- SVs ----------------
if analysis_type in ["SVs", "SVs_and_CNVs"] and ("cosmic" in sv_methods):
    if (sv_source == "ogm") and (sv_format == "smap"):
        snakefiles_to_include.append("rules/1.1_Preprocessing_SVs.smk")
        all_targets.extend([
            f"{resultsdir}/SVs/Patients/{{anonymised}}/Check/check_svs_done.txt",
          #  f"{resultsdir}/SVs/Patients/{{anonymised}}/ProcessData/{{anonymised}}.bedpe",
            f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/data/SigProfilerSVdf.bedpe",
        ])
        snakefiles_to_include.append("rules/2.1_Individualanalysis_SVs.smk")
        all_targets.extend([
            f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.SV32.matrix.tsv",
            f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/Extractor/",
            f"{logsdir}/SVs/Patients/{{anonymised}}/SigProfiler/sigprofilermatrix.log",
        ])
        if run_cohort:
            snakefiles_to_include.append("rules/3.1_Cohortanalysis_SVs.smk")
            all_targets.extend([
                f"{resultsdir}/SVs/Cohort/Cohort.SV32.matrix.tsv",
                f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor/",
                f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor/SV32/Suggested_Solution/COSMIC_SV32_Decomposed_Solution/Activities/COSMIC_SV32_Activities.txt",
                f"{resultsdir}/SVs/Cohort/SigProfiler/COSMIC_SV32_Activities_SV_spiral.svg"
            ])

    elif (sv_source == "ngs") and (sv_format == "bedpe"):
        snakefiles_to_include.append("rules/2.1_Individualanalysis_SVs.smk")
        all_targets.extend([
            f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.SV32.matrix.tsv",
            f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/Extractor/",
            f"{logsdir}/SVs/Patients/{{anonymised}}/SigProfiler/sigprofilermatrix.log",
        ])
        if run_cohort:
            snakefiles_to_include.append("rules/3.1_Cohortanalysis_SVs.smk")
            all_targets.extend([
                f"{resultsdir}/SVs/Cohort/Cohort.SV32.matrix.tsv",
                f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor/SV32/All_solutions_stat.csv",
                f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor/SV32/Suggested_Solution/COSMIC_SV32_Decomposed_Solution/Activities/COSMIC_SV32_Activities.txt",
                f"{resultsdir}/SVs/Cohort/SigProfiler/COSMIC_SV32_Activities_SV_spiral.svg"
            ])
    else:
        raise ValueError(f"SVs config not supported: source={sv_source}, format={sv_format}")


if analysis_type in ["CNVs", "SVs_and_CNVs"] and ("cosmic" in cnv_methods):
    if (cnv_source == "ogm") and (cnv_format in ["csv","txt"]):
        snakefiles_to_include.append("rules/1.2_Preprocessing_CNVs.smk")
        all_targets.extend([
            f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/data/SigProfilerCNVdf.tsv",
            f"{logsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/sigprofilercnvmatrixgenerator.log",
            f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.CNV48.matrix.tsv",
        ])
        snakefiles_to_include.append("rules/2.2_Individualanalysis_CNVs.smk")
        all_targets.extend([
            f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.CNV48.matrix.tsv",
            f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/Extractor/CNV48/All_solutions_stat.csv"
        ])
        if run_cohort:
            snakefiles_to_include.append("rules/3.2_Cohortanalysis_CNVs.smk")
            all_targets.extend([
                f"{resultsdir}/CNVs/Cohort/SigProfiler/MatrixGenerator/Cohort.CNV48.matrix.tsv",
                f"{resultsdir}/CNVs/Cohort/SigProfiler/Extractor/CNV48/All_solutions_stat.csv",
                f"{resultsdir}/CNVs/Cohort/SigProfiler/Extractor/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/Activities/COSMIC_CNV48_Activities.txt",
                f"{resultsdir}/CNVs/Cohort/SigProfiler/COSMIC_CNV48_Activities_CNV_spiral.svg"
            ])

    elif (cnv_source == "ngs") and (cnv_format in ["cns","bed"]):
        # NGS CNV → sin preproc OGM, pero sí SigProfiler CNV48
        snakefiles_to_include.append("rules/2.2_Individualanalysis_CNVs.smk")
        all_targets.extend([
            f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.CNV48.matrix.tsv",
        ])
        if run_cohort:
            snakefiles_to_include.append("rules/3.2_Cohortanalysis_CNVs.smk")
            all_targets.extend([
                f"{resultsdir}/CNVs/Cohort/SigProfiler/MatrixGenerator/Cohort.CNV48.matrix.tsv",
                f"{resultsdir}/CNVs/Cohort/SigProfiler/results/Extractor/CNV48/All_solutions_stat.csv"
            ])
    else:
        raise ValueError(f"CNVs COSMIC config not supported: source={cnv_source}, format={cnv_format}")

# ---------------- CNVs: Drews / Tao (independientes de COSMIC) ----------------
if analysis_type in ["CNVs", "SVs_and_CNVs"]:
    cnv_sig_methods = [m for m in cnv_methods if m in ["drews","tao"]]
    if cnv_sig_methods:
        include: "rules/1.3_Format_CNVs.smk"
        all_targets.extend(
            expand(
                f"{resultsdir}/CNVs/Patients/{{anonymised}}/{{method}}/{{method}}formatted.tsv",
                anonymised=sample_table['anonymised'],
                method=cnv_sig_methods
            )
        )
        if "drews" in cnv_sig_methods:
            include: "rules/Drews.smk"
            all_targets.extend([
                f"{resultsdir}/CNVs/Cohort/Drews/drews_signatures_results.rds",
                f"{resultsdir}/CNVs/Cohort/Drews/activities.png",
                f"{resultsdir}/CNVs/Cohort/Drews/predictions.csv"
            ])
        if "tao" in cnv_sig_methods:
            include: "rules/Tao.smk"
            all_targets.extend([
                f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_signatures_raw.csv",
                f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_exposures_raw.csv",
                f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_signatures.png",
                f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_exposures.png" 
                ])


for snakefile in snakefiles_to_include:
    include: snakefile

# ---- rule all ----
rule all:
    input:
        expand(all_targets, anonymised=sample_table['anonymised'], allow_missing=True)

# helpers (si no los tienes ya arriba)
import os, glob
bionano_dir = config["input"]["bionanodata"]

def sample_folder(anon):
    # carpeta del bruto: p.ej., /.../Revisado_Bionano/Cells_802
    return os.path.join(bionano_dir, anonymised_to_bionano[anon])

rule decompressOGM:
    input:
        ogmdirectory = config['input']['bionanodata'],
        workdirectory = f"{workdirectory}",
        samplesfile   = config['input']['samples']
    output:
        ogm_dones = expand(f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata/.done",
                           anonymised=sample_table['anonymised'])
    params:
        script  = f"{config['directories']['scriptsdir']}/1_DECOMPRESSOGM.R",
        use_zip = zip_flag
    log:
        logfile = f"{logsdir}/DATA/decompressOGM/descompress.log"
    shell:
        r"""
        set -uo pipefail

        if [ "{params.use_zip}" = "True" ]; then
            echo "[decompressOGM] Running decompression with OGM ZIPs..." >> "{log.logfile}"
            Rscript "{params.script}" "{input.ogmdirectory}" "{input.workdirectory}" "{input.samplesfile}"
        else
            echo "[decompressOGM] zip flag is FALSE, mocking OGMdata from DATA/<sample>/*." >> "{log.logfile}"
        fi

        # Para cada muestra: crear OGMdata, copiar/enlazar CNVs y SVs y marcar .done
        awk -v FS='\t' '
          NR==1 {{
            for (i=1;i<=NF;i++) {{
              if ($i=="anonymised") ca=i
              if ($i=="sample")     cs=i
            }}
            next
          }}
          {{print $ca "\t" $cs}}
        ' "{input.samplesfile}" \
        | while read -r anon samp; do
            src_dir="{input.ogmdirectory}/$samp"
            dst_dir="{resultsdir}/DATA/Patients/$anon/OGMdata"
            mkdir -p "$dst_dir"

            if [ "{params.use_zip}" != "True" ]; then
              if [ -d "$src_dir" ]; then
                # copiamos/enlazamos los CNVs (.txt/.csv) y SVs (.smap) desde DATA/SyntheticX
                for f in "$src_dir"/*.txt "$src_dir"/*.csv "$src_dir"/*.smap; do
                  [ -e "$f" ] || continue
                  ln -sf "$f" "$dst_dir/$(basename "$f")"
                done
              else
                echo "[decompressOGM] WARNING: source dir '$src_dir' not found for anonymised=$anon" >> "{log.logfile}"
              fi
            fi

            [[ -f "$dst_dir/.done" ]] || date +"%F %T  OK $anon" > "$dst_dir/.done"
          done
        """




################################################################################################################################################################################
#import subprocess, os, signal, time

#pulpo_dir = os.path.join(workdirectory, ".pulpo"); os.makedirs(pulpo_dir, exist_ok=True)
#pid_anim  = os.path.join(pulpo_dir, "anim.pid")
#pid_tick  = os.path.join(pulpo_dir, "tick.pid")

#def _kill(pidfile):
#    try:
#        if os.path.exists(pidfile):
#            with open(pidfile) as f: pid = int(f.read().strip())
#            os.kill(pid, signal.SIGTERM); time.sleep(0.2); os.remove(pidfile)
#    except Exception:
#        pass

#onstart:
#    if not os.environ.get("PULPO_NO_ANIM"):
#        p1 = subprocess.Popen(["python","animate_pulpo.py","start"])
#        with open(pid_anim,"w") as f: f.write(str(p1.pid))
#        p2 = subprocess.Popen(["python","pulpo_progress.py"])
#        with open(pid_tick,"w") as f: f.write(str(p2.pid))

#onsuccess:
#    _kill(pid_tick); _kill(pid_anim)
#    subprocess.run(["python","animate_pulpo.py","success"])

#onerror:
#    _kill(pid_tick); _kill(pid_anim)
#    subprocess.run(["python","animate_pulpo.py","error"])


