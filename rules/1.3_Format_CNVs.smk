# ==================================================
# File: rules/1.3_Format_CNVs.smk
# Description:
#   Formats Copy Number Variant (CNV) data for each patient
#   into the specific input structure required by CNV signature
#   extraction methods (e.g. Tao, Drews, Steele, Sigminer).
#
#   For each anonymised patient:
#     1️⃣ Loads decompressed OGM CNV data (.done marker ready)
#     2️⃣ Executes CNV_formatting.R to normalize and reformat CNV tables
#     3️⃣ Generates per-method formatted CNV files for downstream analysis
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: CNV_formatting.R
# ==================================================

rule format_cnvs:
    input:
        cnv_done = f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata/.done"
    output:
        formatted = f"{resultsdir}/CNVs/Patients/{{anonymised}}/{{method}}/{{method}}formatted.tsv"
    params:
        script = f"{config['directories']['scriptsdir']}/CNV_formatting.R",
        patient = "{anonymised}",
        base_dir = f"{resultsdir}/CNVs/Patients/{{anonymised}}",
        cnv_dir = f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata"  
    log:
        logfile = f"{logsdir}/CNVs/Patients/{{anonymised}}/{{method}}/format_cnvs.log"
    shell:
        """
        Rscript {params.script} {params.cnv_dir} {params.base_dir} {params.patient} {wildcards.method}
        """

