# ==================================================
# File: rules/2.2_Individualanalysis_CNVs.smk
# Description:
#   Performs Copy Number Variant (CNV) mutational signature analysis
#   at the individual sample level using SigProfilerMatrixGeneratorR
#   and SigProfilerExtractorR.
#
#   For each anonymised patient:
#     1️⃣ Generates the CNV48 mutational matrix using SigProfilerMatrixGeneratorR
#     2️⃣ Extracts de novo CNV mutational signatures using SigProfilerExtractorR
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Dependencies: SigProfilerMatrixGeneratorR, SigProfilerExtractorR
# ==================================================
rule sigprofilermatrixgeneratorcnv:
    threads: 1
    input:
        sigprofilercnv=f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/data/SigProfilerCNVdf.tsv"

    output:
        output=directory(f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/"),
        cnvmatrix=f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.CNV48.matrix.tsv",
        log_file=f"{logsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/sigprofilercnvmatrixgenerator.log",
    params:
        pythondirectory=f"{config['directories']['pythonenvdir']}",
        script=f"{config['directories']['scriptsdir']}/10_Sigprofilermatrixgeneratorcnvs.R",
        log_file = f"{logsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/sigprofilercnvmatrixgenerator.log",
    shell:
        """
        (touch {output.log_file}
        Rscript {params.script} {input.sigprofilercnv} {params.pythondirectory} {output.output} {wildcards.anonymised} {output.cnvmatrix})  > {output.log_file} 2>&1 || echo "..." >> {output.log_file}
        """

rule sigprofilerextractorcnv:
    resources: 
        mem_mb = config["resources"]["per_rule"]["sigprofiler_individual"],
        threads = 4
    input:
        cnvmatrix=f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.CNV48.matrix.tsv",
    output:
        output= directory(f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/Extractor/"),
        stats =  f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/Extractor/CNV48/All_solutions_stat.csv"
    params:
        min_signatures= f"{config['signatures']['min_signatures']}",
        max_signatures= f"{config['signatures']['max_signatures']}",
        nmf_replicates= f"{config['signatures']['replicates']}",
        pythondirectory=f"{config['directories']['pythonenvdir']}",
        script=f"{config['directories']['scriptsdir']}/11_Sigprofilerextractorcnvs.R",
        error_log=f"{logsdir}/CNVs/SigProfilerExtractor/Error_log.tsv", # Centralised error file

    benchmark:
        f"{logsdir}/benchmarks/Patients/{{anonymised}}/cnvsindividualsigprofiler_benchmark.txt"

    shell:
        """
        Rscript --vanilla {params.script} {input.cnvmatrix} {params.pythondirectory} \
        {output.output} {params.min_signatures} {params.max_signatures} \
        {params.nmf_replicates} {params.error_log} || true
        """


