<!-- Logo aligned to the right -->
<p align="right">
  <img src="figs/PULPOlogo.png" width="220">
</p>

<!-- Title and description -->
<h1 align="left"> PULPO  v2.0.0 </h1>
<h3 align="left"> Pipeline of Understanding Large-scale Patterns of Oncogenomic signatures</h3>

<!-- Badges aligned to the right -->
<p align="right">
  <a href="LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-green">
  </a>
  <img alt="Snakemake v7.32.4" src="https://img.shields.io/badge/Snakemake-v7.32.4-blue">
  <img alt="Conda environment" src="https://img.shields.io/badge/Conda-environment-green">
  <a href="https://doi.org/10.1101/2025.07.02.661487">
    <img alt="bioRxiv preprint DOI: 10.1101/2025.07.02.661487"
         src="https://img.shields.io/badge/preprint-bioRxiv-blue">
  </a>
  <a href="(https://doi.org/10.5281/zenodo.17700531)">
    <img alt="Zenodo DOI:"
         src="https://zenodo.org/badge/DOI/10.5281/zenodo.17700531.svg">
  </a>
</p>


###### Implemented by: 
***Marta Portasany-Rodríguez***

*Gonzalo Soria-Alcaide*

*Jorge García-Martínez*

## Introduction

PULPO is a **Snakemake-based workflow** for analysing **structural variants (SVs)** and **copy-number variants (CNVs)** and extracting **mutational signatures** from:

- **Optical Genome Mapping (OGM)** data (Bionano SMAP + CNV export), and
- Optionally, **NGS-based** SV/CNV calls formatted as BEDPE / CNV tables.

The pipeline builds **SV32** and **CNV48** matrices, runs *SigProfiler* for de novo and refitting analyses, and supports **multiple CNV signature catalogues** (COSMIC, Drews et al., Tao et al.), generating cohort-level visualisations and summary tables.

---

## 1. Key features

- **End-to-end SV/CNV workflow**
  - From raw OGM exports (SMAP, CNV txt/csv) or NGS calls to signature activities per sample and per cohort.
- **Multiple signature catalogues**
  - SV: COSMIC SV32.
  - CNV: COSMIC CNV48, Drews CIN signatures, Tao/Sigminer CNV signatures.
- **OGM-aware preprocessing**
  - Automatic organisation of OGM exports into a standard `results/DATA/Patients/<anonymised>/OGMdata` layout.
  - QC of OGM files and robust conversion to SigProfiler-compatible formats.
- **Cohort-level analyses**
  - Merges per-sample matrices into cohort matrices.
  - Runs SigProfilerExtractor for cohort signatures.
  - Generates spiral plots and heatmaps for signature activities across the genome.
- **Config-driven and modular**
  - All behaviour controlled via a YAML config file.
  - Snakemake rules split into logical modules (0_, 1.x, 2.x, 3.x, Drews, Tao) for easy extension.
- **Reproducible and script-based**
  - All heavy lifting done in R and Python scripts under `scripts/`.
  - Same pipeline can be run on OGM and non-OGM cohorts by changing the configuration and input formats.

---

## 2. Repository structure

```text
PULPO/
├── Snakefile                  # Main Snakemake entry point
├── config/
│   ├── configpulpoOGM.yaml    # Main configuration file
│   └── samples_*.tsv          # Sample sheets (OGM, synthetic, etc.)
├── rules/
│   ├── 0_prepare_ogm_data.smk
│   ├── 1.1_Preprocessing_SVs.smk
│   ├── 1.2_Preprocessing_CNVs.smk
│   ├── 1.3_Format_CNVs.smk
│   ├── 2.1_Individualanalysis_SVs.smk
│   ├── 2.2_Individualanalysis_CNVs.smk
│   ├── 3.1_Cohortanalysis_SVs.smk
│   ├── 3.2_Cohortanalysis_CNVs.smk
│   ├── Drews.smk
│   └── Tao.smk
├── scripts/                   # R scripts for conversion, SigProfiler, Drews, Tao, visualisation
├── DATA/                      # Example input data (e.g., synthetic OGM cohort)
├── results/                   # Output directory created by the pipeline
└── logs/                      # Log files per rule and per sample
