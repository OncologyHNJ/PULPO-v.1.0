# ==================================================
# File: rules/1.2_Preprocessing_CNVs.smk
# Description:
#   Preprocesses Optical Genome Mapping (OGM) Copy Number Variant (CNV)
#   export files and converts them into SigProfiler-compatible TSV format.
#
#   For each patient:
#     1️⃣ Detects CNV export files (.txt or .csv) within the OGMdata directory
#     2️⃣ Selects the most recently modified CNV file
#     3️⃣ Executes ogm_cnv_to_sigprofiler.R to generate:
#         results/CNVs/Patients/<anonymised>/SigProfiler/data/SigProfilerCNVdf.tsv
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: ogm_cnv_to_sigprofiler.R
# ==================================================

import os
import glob

rule ogm_cnv_to_sigprofiler_data:
    input:
        ready = f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata/.done"
    output:
        out = f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/data/SigProfilerCNVdf.tsv"
    params:
        script = f"{config['directories']['scriptsdir']}/ogm_cnv_to_sigprofiler.R"
    log:
        log = f"{logsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/ogm_cnv_to_sigprofiler.log"
    run:
        ogm_dir = os.path.dirname(input.ready)
        preferred = glob.glob(os.path.join(ogm_dir, "**", "*CNV*.txt"), recursive=True) \
                  + glob.glob(os.path.join(ogm_dir, "**", "*CNV*.csv"), recursive=True)
        candidates = preferred if preferred else \
            (glob.glob(os.path.join(ogm_dir, "**", "*.txt"), recursive=True)
           + glob.glob(os.path.join(ogm_dir, "**", "*.csv"), recursive=True))

        if not candidates:
            raise FileNotFoundError(f"[CNVs] No CNV .txt/.csv under {ogm_dir}. Does the ZIP include CNV export?")

        candidates.sort(key=lambda p: os.path.getmtime(p), reverse=True)
        cnv_file = candidates[0]

        os.makedirs(os.path.dirname(output.out), exist_ok=True)
        os.makedirs(os.path.dirname(log.log), exist_ok=True)

        shell(f"""
        Rscript {params.script} {cnv_file} {output.out} &> {log.log}
        """)




