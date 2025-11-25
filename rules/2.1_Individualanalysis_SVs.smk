# ==================================================
# File: rules/2.1_Individualanalysis_SVs.smk
# Description:
#   Performs Structural Variant (SV) signature analysis at the individual
#   sample level using SigProfilerMatrixGenerator and SigProfilerExtractor.
#
#   For each anonymised patient:
#     1️⃣ Generates the SV32 mutational matrix using SigProfilerMatrixGeneratorR
#     2️⃣ Extracts de novo SV mutational signatures using SigProfilerExtractorR
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: SigProfilerMatrixGeneratorR, SigProfilerExtractorR
# ==================================================

rule sigprofilermatrixgenerator:
    threads: 1
    input:
        data = f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/data/SigProfilerSVdf.bedpe"
    output:
        output = directory(f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/"),
        data= f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.SV32.matrix.tsv",
    params:
        pythondirectory = f"{config['directories']['pythonenvdir']}",
        script = f"{config['directories']['scriptsdir']}/5_Sigprofilermatrixgenerator.R",
    log:
        log_file = f"{logsdir}/SVs/Patients/{{anonymised}}/SigProfiler/sigprofilermatrix.log"
    shell:
        """
        Rscript {params.script} {input.data} {params.pythondirectory} {output.output} {wildcards.anonymised} &> {log.log_file} || true
        """

rule sigprofilerextractor:
    resources: 
        mem_mb = config["resources"]["per_rule"]["sigprofiler_individual"],
        threads = 4
    input:
        data = f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.SV32.matrix.tsv",
    output:
        output = directory(f"{resultsdir}/SVs/Patients/{{anonymised}}/SigProfiler/results/Extractor/"),

    params:
        pythondirectory = f"{config['directories']['pythonenvdir']}",
        script = f"{config['directories']['scriptsdir']}/6_Sigprofilerextractor.R",
        minimum_signatures = config["signatures"]["min_signatures"],
        maximum_signatures = config["signatures"]["max_signatures"],
        nmf_replicates = config["signatures"]["replicates"]
    log:
        log_file = f"{logsdir}/SVs/Patients/{{anonymised}}/SigProfiler/sigprofilerextractor.log"
    benchmark:
        f"{logsdir}/benchmarks/Patients/{{anonymised}}/svsindividualsigprofiler_benchmark.txt"
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




