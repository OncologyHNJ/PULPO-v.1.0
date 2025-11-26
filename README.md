Summary
PULPO v2.0 is a major update of the pipeline, corresponding to the revised version of the PULPO: Pipeline of understanding large-scale patterns of oncogenomic signatures manuscript. This release introduces new CN signature catalogues and several usability and documentation improvements. It supersedes v1.0, which is now deprecated and kept only for reproducibility.

New features
Added support for multiple CN signature catalogues (COSMIC, Drews et al., Tao et al.).
Included new QC reports and visualisations (spiralize, signature exposure plots, summary tables).
Pipeline and configuration changes
Simplified and modularised the Snakemake workflow for easier extension and maintenance.
Updated configuration files with clearer defaults and comments for OGM and non-OGM cohorts.
Improved logging, error messages and intermediate output structure to facilitate debugging and reuse.
Bug fixes and refactoring
Fixed minor bugs affecting edge cases in CNV/SV formatting and feature generation.
Refactored R and Python helper scripts for clarity, consistency and better performance on larger cohorts.
Cleaned up deprecated options and harmonised naming across modules and outputs.
Reproducibility
This release is associated with the v2.0 code snapshot used in the revised manuscript.
For reproducibility of previous analyses, the original implementation remains available as v1.0 (deprecated).
A versioned, archived snapshot of this release is available via Zenodo (10.5281/zenodo.17306077).
