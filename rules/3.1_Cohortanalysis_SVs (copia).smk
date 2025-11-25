# ==================================================
# File: rules/3.1_Cohortanalysis_SVs.smk
# Description:
#   Performs cohort-level Structural Variant (SV) mutational signature analysis
#   across all anonymised patient samples.
#
#   Workflow overview:
#     1️⃣ Merge per-sample SV32 matrices into a single cohort matrix
#     2️⃣ Extract cohort-level mutational signatures using SigProfilerExtractorR
#     3️⃣ Generate summary visualization with spiralize (SVG + PDF)
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: SigProfilerMatrixGeneratorR, SigProfilerExtractorR, spiralize
# ==================================================

rule mergepatients:
    input:
        log_file = expand(f"{logsdir}/SVs/Patients/{{anonymised}}/SigProfiler/sigprofilermatrix.log", anonymised=sample_table['anonymised'])
    output:
        outputcohort = f"{resultsdir}/SVs/Cohort/Cohort.SV32.matrix.tsv",
    params:
        directorypatients = f"{resultsdir}/SVs/Patients",
        script = f"{config['directories']['scriptsdir']}/7_Mergecohortsvs.R",
    log:
        log_file = f"{logsdir}/SVs/Cohort/mergepatients/mergepatients.log"
    shell:
        """
        Rscript {params.script} {params.directorypatients} {output.outputcohort} &> {log.log_file} || true
        """

rule sigprofilerextractorcohort:
    resources:
    	mem_mb = config["resources"]["per_rule"]["sigprofiler_cohort"],
    	threads = 6
    input:
        data = f"{resultsdir}/SVs/Cohort/Cohort.SV32.matrix.tsv",
    output:
        output = directory(f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor"),
        svfile = f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor/SV32/Suggested_Solution/COSMIC_SV32_Decomposed_Solution/Activities/COSMIC_SV32_Activities.txt"
    params:
        pythondirectory = config['directories']['pythonenvdir'],
        script = f"{config['directories']['scriptsdir']}/8_Sigprofilerextractorcohort.R",
        minimum_signatures = config['signatures']['min_signatures'],
        maximum_signatures = config['signatures']['max_signatures'],
        nmf_replicates = config['signatures']['replicates'],
        output = directory(f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor"),
    log:
        log_file = f"{logsdir}/SVs/Cohort/SigProfiler/sigprofilerextractorcohort.log"
    benchmark:
        f"{logsdir}/benchmarks/svscohortsigprofiler_benchmark.txt"
    shell:
        """
        Rscript --vanilla {params.script} \
            {input.data} \
            {params.pythondirectory} \
            {output.output} \
            {params.minimum_signatures} \
            {params.maximum_signatures} \
            {params.nmf_replicates} \
            &> {log.log_file}
        """

rule spiralize_svs:
    resources:
        mem_mb = config["resources"]["per_rule"]["spiralizeplot"],
        threads = 2
    input:
        svfile = f"{resultsdir}/SVs/Cohort/SigProfiler/Extractor/SV32/Suggested_Solution/COSMIC_SV32_Decomposed_Solution/Activities/COSMIC_SV32_Activities.txt"
    output:
        svg = f"{resultsdir}/SVs/Cohort/SigProfiler/COSMIC_SV32_Activities_SV_spiral.svg",
        pdf = f"{resultsdir}/SVs/Cohort/SigProfiler/COSMIC_SV32_Activities_SV_spiral.pdf"
    params:
        script = f"{config['directories']['scriptsdir']}/spiralize.R",
        datatype = "SV"
    log:
        f"{logsdir}/SVs/Cohort/spiralize_SVs.log"
    shell:
        "Rscript {params.script} {input.svfile} {params.datatype} {output.svg} {output.pdf} > {log} 2>Rscript {params.script} {input.svfile} {params.datatype} > {log} 2>&11"

