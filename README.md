# Pseudomonas CRISPRi Library Analysis

Genome-wide CRISPRi screen in *P. aeruginosa* PA14 during murine pneumonia infection. Found 178 host-essential genes, validating ispD and pgsA as key targets.

## Pipeline
- Raw data processing (>5M reads, >85% mapping)
- edgeR v4.2.0 statistical analysis
- Population bottleneck analysis (Abel et al. 2015)
- STRING-DB v12.0.2 enrichment

## Requirements
edgeR 4.2.0, data.table 1.15.4, ggplot2 3.5.1, limma 3.60.0

## Data
- Sequencing: NCBI BioProject PRJNA1178310 
- Annotations: [STRING-DB PA14](https://version-12-0.string-db.org/organism/STRG0A01FJP)
