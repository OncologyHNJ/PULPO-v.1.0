# ============================================================
# Module: Tao CNV Signature Extraction
# File: rules/tao.smk
# ============================================================
# Description:
#   Executes the CNV signature extraction method described by
#   Tao et al. (2023) using the Sigminer framework.
#
#   This step is optional and is executed only when 
#   `cnvsignatures.method: tao` is enabled in the configuration
#   file (`configpulpoOGM.yaml`).
#
#   The rule aggregates all individual patient CNV-formatted 
#   files (`taoformatted.tsv`) and performs CNV signature
#   discovery and exposure quantification across the cohort.
#
# Notes:
#   - This module identifies de novo CNV signatures using 
#     non-negative matrix factorization (NMF) implemented 
#     in the Sigminer R package, following Tao et al. (2023).
#   - Outputs include extracted signature profiles, exposure
#     matrices, and diagnostic plots.
#
# Author: Marta Portasany-Rodríguez
# Pipeline: PULPO 🐙 v1.0
# Last modified: 2025-11-05
# ============================================================


rule tao:
    resources:
        mem_mb = config["resources"]["per_rule"]["taomethod"],
        threads = 4
    input:
        expand(f"{resultsdir}/CNVs/Patients/{{anonymised}}/tao/taoformatted.tsv",
               anonymised=sample_table['anonymised'])
    output:
        output_dir = directory(f"{resultsdir}/CNVs/Cohort/Sigminer/"),
        signatures = f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_signatures_raw.csv",
        exposures = f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_exposures_raw.csv",
        sig_plot = f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_signatures.png",
        exp_plot = f"{resultsdir}/CNVs/Cohort/Sigminer/Tao_exposures.png"
    params:
        base_dir = f"{resultsdir}/CNVs/Patients/",
        script = f"{config['directories']['scriptsdir']}/taomethod.R",
        nruns = f"{config['signatures']['replicates']}",
    log:
        logfile = f"{logsdir}/CNVs/Cohort/tao_method.log"

    benchmark:
        f"{logsdir}/benchmarks/tao_benchmark.txt"

    shell:
        """
        Rscript {params.script} {params.base_dir} {output.output_dir} {params.nruns} 
        """

