# ============================================================
# File: rules/Drews.smk
# Description:
#   Executes the Drews et al. (2022) methodology for extracting 
#   CNV mutational signatures at the cohort level. 
#
#   This step is optional and is only executed if 
#   `cnvsignatures.method: drews` is enabled in the 
#   configuration file (`configpulpoOGM.yaml`).
#
#   The rule aggregates all individual patient CNV-formatted 
#   files (`drewsformatted.tsv`) and runs the CNV signature 
#   extraction and prediction workflow implemented in 
#   `drewsmethod.R`.
#
# Author: Marta Portasany-Rodríguez
# Pipeline: PULPO 🐙 
# Last modified: 2025-11-05
# ============================================================

rule drews:
    resources:
        mem_mb = config["resources"]["per_rule"]["drewsmethod"],
        threads = 4
    input:
        expand(f"{resultsdir}/CNVs/Patients/{{anonymised}}/drews/drewsformatted.tsv",
               anonymised=sample_table['anonymised'])
    output:
        output_dir = directory(f"{resultsdir}/CNVs/Cohort/Drews/"),
        rds = f"{resultsdir}/CNVs/Cohort/Drews/drews_signatures_results.rds",
        activities = f"{resultsdir}/CNVs/Cohort/Drews/activities.png",
        predictions = f"{resultsdir}/CNVs/Cohort/Drews/predictions.csv"
    params:
        base_dir = f"{resultsdir}/CNVs/Patients/",
        script = f"{config['directories']['scriptsdir']}/drewsmethod.R"
    log:
        logfile = f"{logsdir}/CNVs/Cohort/drews_method.log"
    benchmark:
        f"{logsdir}/benchmarks/drews_benchmark.txt"
    shell:
        """
        Rscript {params.script} {params.base_dir} {output.output_dir}
        """

