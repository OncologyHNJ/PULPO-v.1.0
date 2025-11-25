# ==================================================
# File: 3.2_Cohortanalysis_CNVs.smk
# Description:
#   Performs cohort-level CNV signature analysis using SigProfilerExtractor.
#   This module merges all single-sample CNV48 matrices, executes SigProfilerExtractor
#   in cohort mode, and generates a spiral plot summarizing CNV signature activities
#   across all samples.
#
# Workflow summary:
#   CNV48 matrices (per sample)
#        └──▶ merged into cohort matrix
#                └──▶ SigProfilerExtractor (cohort mode)
#                        └──▶ Spiral visualization (SVG + PDF)
#
# Author: Marta Portasany-Rodríguez
# Created on: 2025-11-05
# Pipeline: PULPO 🐙 v1.0
# Dependencies:
#   - R ≥ 4.2
#   - SigProfilerExtractorR
#   - spiralize, data.table, ComplexHeatmap, grid, grDevices
#   - Snakemake ≥ 7.0
# ==================================================


rule mergepatientcnvs:
    input:
        cnvmatrix = expand(f"{resultsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/results/MatrixGenerator/{{anonymised}}.CNV48.matrix.tsv",anonymised=sample_table['anonymised']),
        log_file = expand(f"{logsdir}/CNVs/Patients/{{anonymised}}/SigProfiler/sigprofilercnvmatrixgenerator.log", anonymised=sample_table['anonymised']),
    output:
        outputcohort =  f"{resultsdir}/CNVs/Cohort/SigProfiler/MatrixGenerator/Cohort.CNV48.matrix.tsv",
    params:
        directorypatients = f"{resultsdir}/CNVs/Patients",
        script = f"{config['directories']['scriptsdir']}/12_Mergecohortcnvs.R",
    log:
        log_file = f"{logsdir}/CNVs/Cohort/mergepatients/mergepatientscnvs.log"
    shell:
        """
        touch {output.outputcohort}
        Rscript {params.script} {params.directorypatients} {output.outputcohort} &> {log.log_file} || true
        """

rule sigprofilerextractorcnvscohort:
    resources:
    	mem_mb = config["resources"]["per_rule"]["sigprofiler_cohort"],
    	threads = 6
    input:
        cnvmatrix =  f"{resultsdir}/CNVs/Cohort/SigProfiler/MatrixGenerator/Cohort.CNV48.matrix.tsv",
    output:
        output= directory(f"{resultsdir}/CNVs/Cohort/SigProfiler/Extractor/"),
        stats = f"{resultsdir}/CNVs/Cohort/SigProfiler/Extractor/CNV48/All_solutions_stat.csv",
        cnvfile = f"{resultsdir}/CNVs/Cohort/SigProfiler/Extractor/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/Activities/COSMIC_CNV48_Activities.txt"
        
    params:
        pythondirectory = f"{config['directories']['pythonenvdir']}",
        script =  f"{config['directories']['scriptsdir']}/13_Sigprofilerextractorcnvscohort.R",
        min_signatures = f"{config['signatures']['min_signatures']}",
        max_signatures = f"{config['signatures']['max_signatures']}",
        nmf_replicates = f"{config['signatures']['replicates']}",

    log:
        log_file= f"{logsdir}/CNVs/Cohort/SigProfiler/sigprofilercnvextractorcnvscohort.log"

    benchmark:
        f"{logsdir}/benchmarks/cnvscohortsigprofiler_benchmark.txt"

    shell:
        """
        Rscript --vanilla {params.script} \
            {input.cnvmatrix} \
            {params.pythondirectory} \
            {output.output} \
            {params.min_signatures} \
            {params.max_signatures} \
            {params.nmf_replicates} \
            &> {log}
        """
        


rule spiralize_cnvs:
    resources:
        mem_mb = config["resources"]["per_rule"]["spiralizeplot"],
        threads = 2
    input:
        cnvfile = f"{resultsdir}/CNVs/Cohort/SigProfiler/Extractor/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/Activities/COSMIC_CNV48_Activities.txt"
    output:
        svg = f"{resultsdir}/CNVs/Cohort/SigProfiler/COSMIC_CNV48_Activities_CNV_spiral.svg",
        pdf = f"{resultsdir}/CNVs/Cohort/SigProfiler/COSMIC_CNV48_Activities_CNV_spiral.pdf"
    params:
        script = f"{config['directories']['scriptsdir']}/spiralize.R",
        datatype = "CNV"
    log:
        f"{logsdir}/CNVs/Cohort/spiralize_CNVs.log"
    shell:
        "Rscript {params.script} {input.cnvfile} {params.datatype} {output.svg} {output.pdf} > {log} 2>&1"


