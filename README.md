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

PULPO v2.0 is a major update of our Snakemake-based pipeline for analysing **structural variants (SVs)** and **copy number variants (CNVs)** extracting **mutational signatures** from both **Optical Genome Mapping (OGM)** and **NGS-based** data. 

Unlike v1.0 (now deprecated and kept only for reproducibility), PULPO v2.0:

- Supports **multiple CNV signature catalogues** (COSMIC CNV48, Drews CIN signatures, Tao/Sigminer CNV signatures).
- Handles **SV and CNV data from OGM or NGS callers** via flexible input formats (SMAP/BEDPE for SVs, CNV exports/cns/bed for CNVs).
- Includes **cohort-level analyses and visualisations** (spiral plots, exposure heatmaps, summary tables).
- Uses a **modular Snakemake design** and a single YAML configuration file to control the full workflow.

From raw SV/CNV calls to interpretable mutational signatures, PULPO v2.0 aims to make advanced structural and copy-number signature analysis accessible and reproducible, while remaining compatible with paediatric and non-paediatric cancer cohorts.

## 📚 Table of Contents

1. [🚀 Quick Start](#quick-start)  
2. [🛠️ Installation](#installation)  
3. [⚙️ Configuration](#configuration)  
   - [📁 Input Files and Directories](#1-input-files-and-directories)  
   - [📂 Directory Paths](#2-directory-paths)  
   - [🧪 Analysis Configuration](#3-analysis-configuration)  
   - [📄 Input Sources and Formats](#4-input-sources-and-formats)  
   - [🧬 Signature Catalogues](#5-signature-catalogues)  
   - [💻 Resources](#6-resources)  
4. [📈 Usage](#usage)  
   - [🎯 Running Specific Steps](#running-specific-steps)  
   - [🧬 Main Rules](#main-rules)  
   - [🐞 Error Handling and Debugging](#error-handling-and-debugging)  
5. [🗂 Repository Structure](#repository-structure)  
6. [🔗 External Tools and References](#external-tools-and-references)  
7. [📬 Contact](#contact)  
8. [💬 Community Support](#community-support-)  
9. [📄 License](#license)  
10. [🧾 Citation](#citation)

---
## Quick Start 🚀

Minimal example to get PULPO v2.0 running:

```bash
# 1. Clone the repository
git clone https://github.com/OncologyHNJ/PULPO-v.2.0.0.git
cd PULPO-v.2.0.0

# 2. Create a Conda environment with Snakemake
conda create -n PULPO python=3.10 snakemake=7.32.4 -c conda-forge -c bioconda
conda activate PULPO

# 3. Edit the main config file to point to your data
nano config/configpulpoOGM.yaml

# 4. Run the pipeline
snakemake --cores 4
```

This will:

- Read `config/configpulpoOGM.yaml` and your sample sheet.
- Prepare OGM/NGS data.
- Perform SV and/or CNV preprocessing.
- Run per-sample SigProfiler analyses.
- Optionally run cohort-level COSMIC, Drews and Tao CNV signatures (if enabled).

## Installation 🛠️

To run PULPO, you need:

- **Python 3.8+**  
- **Snakemake** (≥ 7.x)  
- **R** (≥ 4.x) with the following key packages:
  - `SigProfilerMatrixGeneratorR`
  - `SigProfilerExtractor`
  - `CINSignatureQuantification`
  - `limSolve`
  - `sigminer`
  - `spiralize`
  - plus standard tidyverse / `data.table` dependencies

### 1. Clone the repository

```bash
git clone https://github.com/OncologyHNJ/PULPO-v.2.0.0.git
cd PULPO-v.2.0.0
```
### 2. Install Conda and Snakemake

If you don’t have Conda installed, you can download Miniconda from:
https://docs.conda.io/en/latest/miniconda.html

Then install Snakemake in a new environment:

```bash
conda create -n PULPO python=3.10 snakemake=7.32.4 -c conda-forge -c bioconda
conda activate PULPO
```
### 3. Install R dependencies

From within the activated environment:

```R

Inside R:

install.packages(c("data.table", "dplyr", "ggplot2"))

# Install BiocManager and then SigProfilerMatrixGeneratorR, SigProfilerExtractor, etc.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("SigProfilerMatrixGeneratorR"))

# Install additional packages such as:
# CINSignatureQuantification, sigminer, spiralize, ...
# (see manuscript / documentation for the full list)
```

You may choose to manage R packages in a separate R/Conda environment if preferred.

## Configuration ⚙️

The pipeline is controlled by a single YAML configuration file located in config/:

- configpulpoOGM.yaml – main configuration file for PULPO v2.0.

You should copy and edit this file to adapt it to your own cohort.

### 1. Input Files and Directories
```yaml
input:
  samples: "/path/to/config/samples_mycohort.tsv"
  bionanodata: "/path/to/OGM_or_CNV_data_root"
- input.samples: path to a TSV file with your sample metadata.
    - Must contain at least a sample column.
    - If anonymised is missing or empty, PULPO can fall back to using sample as anonymised ID.

- input.bionanodata: root directory containing per-sample folders with OGM exports or CNV/SV files (for OGM use cases).
```

For NGS-only analyses you can still use the same sample sheet and point PULPO to directories with BEDPE/CNV tables.

### 2. Directory Paths

```yaml
directories:
  workdirectory: "/path/to/PULPO-v.2.0.0"
  pythonenvdir: "/path/to/python3.10"
  scriptsdir: "/path/to/PULPO-v.2.0.0/scripts"
```
- directories.workdirectory: base working directory where results/ and logs/ will be created.

- directories.pythonenvdir: path to the Python interpreter used by the R reticulate bridge (needed by SigProfilerMatrixGeneratorR).

- directories.scriptsdir: path to the scripts/ folder inside this repository.

⚠️ After cloning, update the default /home/user/... paths to match your local system.


### 3. Analysis Configuration
```yaml

analysis:
  zip: false
  analysis_type: "SVs_and_CNVs"   # "SVs", "CNVs" or "SVs_and_CNVs"
  run_cohort_analysis: true
```

- analysis.zip:

    - true → OGM exports are zipped; PULPO expects a prior decompression step.

    - false → OGM exports are already decompressed (one folder per sample).

- analysis.analysis_type:

    - "SVs" → only structural variants.

    - "CNVs" → only copy-number variants.

    - "SVs_and_CNVs" → both branches are executed.

- analysis.run_cohort_analysis:

    - true → run individual and cohort-level modules.
    - false → only per-sample analysis is performed.
    - 
### 4. Input Sources and Formats

PULPO v2.0 can work with both OGM and NGS calls:

```yaml
inputs:
  SVs:
    source: "ogm"        # "ogm" or "ngs"
    format: "smap"       # ogm: "smap" ; ngs: "bedpe"
  CNVs:
    source: "ogm"        # "ogm" or "ngs"
    format: "csv"        # ogm: "csv"/"txt" ; ngs: "cns"/"bed"
```
- For OGM SVs:

    - source: "ogm", format: "smap" → .smap files will be converted to BEDPE.

- For NGS SVs:

    - source: "ngs", format: "bedpe" → PULPO uses BEDPE directly.

- For OGM CNVs:

    - source: "ogm", format: "csv" (or "txt") → CNV export tables are reformatted.

- For NGS CNVs:

    - source: "ngs", format: "cns" or "bed" depending on your CNV caller.
 
### 5. Signature Catalogue

You can enable different SV/CNV signature frameworks:
```yaml
pipelines:
  SVs:
    methods: ["COSMIC"]                # COSMIC SV32
  CNVs:
    methods: ["COSMIC", "drews", "tao"]  # any subset of COSMIC, drews, tao
```
- SVs:

    - Currently supports COSMIC SV32 signatures via SigProfiler.

- CNVs:

    - COSMIC CNV48 via SigProfiler.

    - Drews CIN signatures via CINSignatureQuantification.

    - Tao CNV signatures via sigminer.

PULPO can run one or multiple CNV catalogues in parallel on the same cohort.

### 6. Resources

```yaml
resources:
  total_mem_mb: 100000
  total_cores: 4
  per_rule:
    sigprofiler_individual: 16000
    sigprofiler_cohort: 64000
    spiralizeplot: 16000
    taomethod: 24000
    drewsmethod: 24000
```
- total_mem_mb and total_cores describe your machine/cluster limits.
- per_rule entries can be used in Snakemake resources: clauses, especially for memory-hungry modules such as SigProfilerExtractor, Drews and Tao methods.
