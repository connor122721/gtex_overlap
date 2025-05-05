# GTEx eQTL Overlap Analysis Pipeline

This pipeline analyzes the overlap between disease eQTLs (TOPCHEF) and GTEx eQTLs across different tissue types.

## Pipeline Overview

- Reads the TOPCHEF eQTL data and processes GTEx eQTLs
- Calculates the overlap between TOPCHEF and GTEx eQTLs
- Colocalization with GTEx eQTL and HF GWAS
- Generates overlap statistics and visualizations

## Usage and prerequisites

To run the pipeline, you will need both Apptainer v1.3.4 and Nextflow v23.04.1:

```bash
nextflow run main.nf -profile slurm
```
