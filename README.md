# GTEx eQTL Overlap Analysis Pipeline

This pipeline analyzes the overlap between disease eQTLs (TOPCHEF) and GTEx eQTLs across different tissue types.

## Pipeline Overview

1. Reads the TOPCHEF eQTL data
2. Processes GTEx eQTL data from Parquet files
3. Calculates the overlap between TOPCHEF and GTEx eQTLs
4. Generates overlap statistics and visualizations

## Usage

Run the pipeline with default parameters:

```bash
nextflow run main.nf
```

Run the pipeline with custom parameters:

```bash
nextflow run main.nf --gtex_dir /path/to/gtex_files --chef_file /path/to/chef_file.rds --output_dir /path/to/output
```

## Parameters

- `gtex_dir`: Directory containing GTEx eQTL Parquet files
- `chef_file`: Path to TOPCHEF eQTL RDS file
- `output_dir`: Directory to save output files

## Output

The pipeline produces:

1. `overlapGTEX_eqtl_results.csv`: CSV file with overlap statistics for each GTEx tissue
2. `overlap_plot.png`: Horizontal bar plot visualizing the overlap percentages

## Author

Connor Murray
