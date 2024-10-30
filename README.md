# Pseudomonas CRISPRi Library Analysis

## Overview
This repository contains the complete analysis pipeline for processing and analyzing CRISPRi library screening data in *Pseudomonas aeruginosa* PA14. The study investigates essential gene perturbations that create vulnerabilities in the mammalian host environment through a comprehensive CRISPRi knockdown library screen in a murine pneumonia model. We identified 178 genes showing significant fitness defects in mice relative to axenic culture, with particular focus on validating the roles of ispD (isoprenoid precursor biosynthesis) and pgsA (phospholipid synthesis).

## Data Analysis Pipeline

### 1. Sequencing Data Processing
The pipeline begins with processing raw amplicon sequencing data from the CRISPRi library screens:
- Quality control and filtering of sequencing data (threshold: >5 million reads, >85% mapping rate)
- Processing of dual-indexed samples from multiple experimental conditions
- Guide RNA sequence abundance quantification and normalization
- Integration of technical and biological replicates

### 2. Statistical Analysis
Comprehensive statistical framework for identifying significant fitness effects:
- Implementation of edgeR (v4.2.0) with quasi-likelihood negative binomial models
- Multi-condition comparison across:
 - Initial inoculum baseline
 - In vitro growth conditions
 - In vivo host environment (high-inoculum samples)
- Integration of multiple guides per gene using Stouffer's method
- Robust statistical testing with FDR correction for multiple comparisons
- Generation of normalized count matrices and differential abundance statistics

### 3. Population Analysis
Detailed assessment of population dynamics and bottleneck effects:
- Implementation of frequency-based bottleneck calculations (Abel et al. 2015)
- Tracking guide frequency changes between experimental conditions
- Calculation of effective population sizes (Nb) with error estimation
- Separate analysis tracks for control and knockdown guide populations
- Visualization of population complexity across experimental stages

### 4. Functional Analysis
Integration with biological pathway and functional annotations:
- Gene-set enrichment analysis using STRING-DB (v12.0.2) annotations
- Competitive gene set testing implementation using camera (limma v3.60.0)
- Integration with P. aeruginosa UCBPP-PA14 proteome functional data
- Enrichment analysis of Gene Ontology terms, protein domains, and pathways
- Visualization of functional networks and enriched pathways

## Dependencies
### R packages:
- edgeR (v4.2.0): Differential abundance analysis
- data.table (v1.15.4): Efficient data manipulation
- ggplot2 (v3.5.1): Data visualization
- limma (v3.60.0): Gene set testing and functional analysis

## Data Availability
- Raw sequencing data: NCBI BioProject PRJNA1178310
- Proteome annotations: [STRING-DB PA14](https://version-12-0.string-db.org/organism/STRG0A01FJP)
