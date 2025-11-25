# ==================================================
# File: rules/1.1_Preprocessing_SVs.smk
# Description:
#   Structural Variant (SV) preprocessing module for the PULPO pipeline.
#   This stage performs:
#     1️⃣ Validation of OGM .smap files to ensure integrity (rule: check_svs)
#     2️⃣ Conversion of validated .smap files to SigProfiler-compatible BEDPE
#         format for mutational signature extraction (rule: smap_to_sigprofiler)
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: CheckfilesSV.R, ogm_to_sigprofiler.R
# ==================================================

rule check_svs:
    input:
        # Wait for the sentinel per patient (produced by descompressOGM)
        ready = f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata/.done",
    output:
        output = f"{resultsdir}/SVs/Patients/{{anonymised}}/Check/check_svs_done.txt",
    params:
        ogmdir = f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata/",
        script = f"{config['directories']['scriptsdir']}/CheckfilesSV.R",
        general_log = f"{resultsdir}/SVs/Checkpoints/check_svs_all.txt"
    log:
        f"{logsdir}/SVs/Patients/{{anonymised}}/Check/check_svs.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {resultsdir}/SVs/Checkpoints "$(dirname {log})"
        echo "Processing sample {wildcards.anonymised}..." >> {params.general_log}
        if ! Rscript "{params.script}" "{wildcards.anonymised}" "{params.ogmdir}" >> "{params.general_log}" 2>&1; then
            echo "Sample {wildcards.anonymised} failed. See {params.ogmdir}/check_svs_error.txt" >> "{params.general_log}"
            exit 1
        fi
        echo "Finished processing sample {wildcards.anonymised}" >> "{params.general_log}"
        touch "{output.output}"
        """

rule smap_to_sigprofiler:
    input:
        ready = f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata/.done",
    output:
        out = f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/data/SigProfilerSVdf.bedpe"
    params:
        bionanodir = f"{resultsdir}/DATA/Patients/",
        script = f"{config['directories']['scriptsdir']}/ogm_to_sigprofiler.R"
    log:
        f"{logsdir}/SVs/Patients/{{anonymised}}/SigProfiler/smap_to_sigprofiler.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.out})" "$(dirname {log})"
        Rscript "{params.script}" "{params.bionanodir}" "{wildcards.anonymised}" "{output.out}" &> "{log}"
        """

