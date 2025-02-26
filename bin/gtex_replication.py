#!/usr/bin/env python3
# Connor Murray - Modified by GitHub Copilot
# Analyzing overlap of disease eQTLs and GTEx across tissue types

import os
import argparse
import pandas as pd
import numpy as np
import pyarrow.parquet as pq
import pyreadr as pr
import time
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze overlap of disease eQTLs and GTEx across tissue types')
    parser.add_argument('--gtex_dir', required=True, help='Directory containing GTEx eQTL Parquet files')
    parser.add_argument('--chef_file', required=True, help='Path to TOPCHEF eQTL RDS file')
    parser.add_argument('--output_dir', required=True, help='Directory to save output files')
    return parser.parse_args()

def process_file(file, chef_df, output_file):
    """
    Process a single GTEx eQTL file and compute overlap with TOPCHEF eQTLs
    """
    start_time = time.time()
    
    # Load only the `variant_id` column
    gtex_df = pq.read_table(file, columns=['variant_id']).to_pandas()
    gtex_df = gtex_df.assign(
        chrom=gtex_df['variant_id'].str.split("_").str[0],
        position=gtex_df['variant_id'].str.split("_").str[1],
        snp=lambda x: x['chrom'] + ":" + x['position'] + "[b37]")
    
    # Calculate proportion of variant IDs in GTEx SNPs
    prop_table = chef_df['variant_id'].isin(gtex_df['snp'])
    proportion = np.round(prop_table.value_counts(normalize=True)*100, 2)

    # Ensure both True and False proportions are present
    true_prop = proportion.get(True, 0.0)
    false_prop = proportion.get(False, 0.0)
    
    # Message
    print(file, true_prop)
    elapsed_time = time.time() - start_time
    print(f"Processed {os.path.basename(file)} in {elapsed_time:.2f} seconds.")
    
    # Write result to the output file
    with open(output_file, "a") as f:
        f.write(f"{os.path.basename(file)},{true_prop},{false_prop}\n")

    return {"file": os.path.basename(file), "snp_overlap_true": true_prop, "snp_overlap_false": false_prop}

def main():
    args = parse_args()
    
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Output file path
    output_file = os.path.join(args.output_dir, "overlapGTEX_eqtl_results.csv")
    
    # GTEx eQTLs V10 list of Parquet files
    gtex_files = [os.path.join(args.gtex_dir, f) 
                  for f in os.listdir(args.gtex_dir) 
                  if f.endswith(".parquet")]
    
    print(f"Found {len(gtex_files)} GTEx eQTL files")
    
    # Find indices of files containing: heart
    heart_indices = [i for i, file in enumerate(gtex_files) if "heart" in file.lower()]
    heart_files = [gtex_files[i] for i in heart_indices]
    print(f"Heart-related files: {len(heart_files)}")
    
    # TOPCHEF eQTLs - mapped using tensorqtl
    print(f"Reading TOPCHEF eQTL file: {args.chef_file}")
    chef_df = pr.read_r(args.chef_file)
    chef_df = chef_df[None]  # Pull out relevant data frame from RDS
    
    # Write header to the output file
    with open(output_file, "w") as f:
        f.write("file,true,false\n")
    
    # Process files sequentially
    results = []
    for file in gtex_files:
        results.append(process_file(file, chef_df, output_file))
    
    # Convert results to a DataFrame
    proportions_df = pd.DataFrame(results)
    
    # Simplify name of tissue
    proportions_df['simp'] = proportions_df['file'].str.replace('.v10.eQTLs.signif_pairs.parquet', '')
    
    # Sort
    proportions_df = proportions_df.sort_values('snp_overlap_true', ascending=True)
    
    # Create a stacked bar plot
    plt.figure(figsize=(8, 16))
    
    # Add bars for the "true" proportion
    plt.barh(proportions_df["simp"], proportions_df["snp_overlap_true"], color="lightblue")
    
    # Formatting the plot
    plt.title("TOPChef Overlap with GTEx v10 eQTLs")
    plt.xlabel("Percentage", fontsize=16)
    plt.ylabel("GTEx Tissue", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(args.output_dir, "overlap_plot.png"))
    print(f"Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()

