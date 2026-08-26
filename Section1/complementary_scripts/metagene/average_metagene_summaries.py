#!/usr/bin/env python3
"""
average_metagene_summaries.py

Averages multiple metagene summary CSV files (the outputs of summarize_metagene_coverage.py)
into a single consolidated metagene coverage file.

It averages the 100 bin columns (bin_1 ... bin_100) and the mean_n_transcripts column
grouped by each length_bin. It also sums or sets the n_samples to represent the total
number of samples.

Usage:
  python3 average_metagene_summaries.py --inputs <file1.csv> <file2.csv> ... --out <average_summary.csv>
Or using a FOFN of summary CSVs:
  python3 average_metagene_summaries.py --fofn <summaries.fofn> --out <average_summary.csv>
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Average multiple length-binned metagene summaries.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--inputs", nargs="+", help="Space-separated paths to individual metagene summary CSV files")
    group.add_argument("--fofn", help="Path to a FOFN containing paths to individual metagene summary CSV files (one per line)")
    parser.add_argument("--out", required=True, help="Path to save the consolidated average summary CSV")
    args = parser.parse_args()
    
    # 1. Collect Input Paths
    input_paths = []
    if args.inputs:
        input_paths = args.inputs
    elif args.fofn:
        if not os.path.exists(args.fofn):
            print(f"Error: FOFN file not found: {args.fofn}", file=sys.stderr)
            sys.exit(1)
        with open(args.fofn, "r") as f:
            for line in f:
                line_clean = line.strip().split("\t")[0].strip() # Handles tab-separated metadata if present
                if line_clean and not line_clean.startswith("#"):
                    input_paths.append(line_clean)
                    
    # Validate input files exist
    valid_paths = []
    for path in input_paths:
        if os.path.exists(path):
            valid_paths.append(path)
        else:
            print(f"Warning: Input file not found: {path}. Skipping...", file=sys.stderr)
            
    if not valid_paths:
        print("Error: No valid input summary CSV files found.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Averaging {len(valid_paths)} metagene summary files...")
    
    # 2. Load and Combine Data
    all_dfs = []
    for path in valid_paths:
        try:
            df = pd.read_csv(path)
            # Basic validation of columns
            required_cols = ["length_bin"] + [f"bin_{i}" for i in range(1, 101)]
            for col in required_cols:
                if col not in df.columns:
                    print(f"Error: File {path} is missing required column: {col}", file=sys.stderr)
                    sys.exit(1)
            all_dfs.append(df)
        except Exception as e:
            print(f"Error reading file {path}: {e}", file=sys.stderr)
            sys.exit(1)
            
    # Combine all dataframes
    combined = pd.concat(all_dfs, ignore_index=True)
    
    # Get bin columns
    bin_cols = [f"bin_{i}" for i in range(1, 101)]
    
    # 3. Group and Average
    # If the files have 'mean_n_transcripts', we average them. If not, default to 0.
    agg_dict = {
        "n_samples": "sum"
    }
    if "mean_n_transcripts" in combined.columns:
        agg_dict["mean_n_transcripts"] = "mean"
    elif "n_transcripts" in combined.columns:
        # Compatibility with different headers
        combined["mean_n_transcripts"] = combined["n_transcripts"]
        agg_dict["mean_n_transcripts"] = "mean"
        
    # Add bin columns to aggregation (using mean)
    for col in bin_cols:
        agg_dict[col] = "mean"
        
    # Group by length_bin and apply aggregations
    summary_df = combined.groupby("length_bin", as_index=False).agg(agg_dict)
    
    # 4. Standardize Column Order
    col_order = ["length_bin", "n_samples", "mean_n_transcripts"] + bin_cols
    # If mean_n_transcripts was not present, we create it
    if "mean_n_transcripts" not in summary_df.columns:
        summary_df["mean_n_transcripts"] = 0
        
    summary_df = summary_df[col_order]
    
    # Write Output
    print(f"Saving consolidated average summary to: {args.out}")
    # Ensure parent dir exists
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    summary_df.to_csv(args.out, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
