#!/usr/bin/env python3
"""
summarize_metagene_coverage.py

Aggregates individual transcript coverage profiles into metagene summaries,
optionally stratifying by transcript length intervals.

Usage:
  python3 summarize_metagene_coverage.py --input <profiles.csv> --out <summary.csv> \
      [--intervals "0-1000,1000-2000,2000-3000,3000-6000,>=6000"]
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np

def parse_intervals(intervals_str):
    """
    Parses a comma-separated list of interval strings and returns a list of
    tuples: (label, filter_function).
    
    Supported formats:
      - 'min-max' (e.g. '0-1000' -> L <= 1000; '1000-2000' -> 1000 < L <= 2000)
      - '>=val' (e.g. '>=6000' -> L >= 6000)
      - '>val' (e.g. '>6000' -> L > 6000)
      - '<=val' (e.g. '<=1000' -> L <= 1000)
      - '<val' (e.g. '<1000' -> L < 1000)
    """
    intervals = []
    if not intervals_str:
        return None
        
    chunks = [c.strip() for c in intervals_str.split(",") if c.strip()]
    for chunk in chunks:
        if chunk.startswith(">="):
            val = float(chunk[2:])
            intervals.append((chunk, lambda L, v=val: L >= v))
        elif chunk.startswith(">"):
            val = float(chunk[1:])
            intervals.append((chunk, lambda L, v=val: L > v))
        elif chunk.startswith("<="):
            val = float(chunk[2:])
            intervals.append((chunk, lambda L, v=val: L <= v))
        elif chunk.startswith("<"):
            val = float(chunk[1:])
            intervals.append((chunk, lambda L, v=val: L < v))
        elif "-" in chunk:
            parts = chunk.split("-")
            if len(parts) == 2:
                try:
                    min_val = float(parts[0])
                    max_val = float(parts[1])
                    if min_val == 0:
                        intervals.append((chunk, lambda L, mn=min_val, mx=max_val: L <= mx))
                    else:
                        intervals.append((chunk, lambda L, mn=min_val, mx=max_val: mn < L <= mx))
                except ValueError:
                    print(f"Warning: Could not parse interval chunk: {chunk}", file=sys.stderr)
            else:
                print(f"Warning: Could not parse interval chunk: {chunk}", file=sys.stderr)
        else:
            print(f"Warning: Unknown interval format: {chunk}", file=sys.stderr)
            
    return intervals

def main():
    parser = argparse.ArgumentParser(description="Summarize transcript coverage profiles into metagene curves.")
    parser.add_argument("--input", required=True, help="Path to input transcript profiles CSV (from calculate_transcript_coverage.py)")
    parser.add_argument("--out", required=True, help="Path to output summarized metagene CSV")
    parser.add_argument("--intervals", help="Optional comma-separated length intervals (e.g., '0-1000,1000-2000,2000-3000,3000-6000,>=6000')")
    parser.add_argument("--method", choices=["mean", "sum"], default="mean",
                        help="Aggregation method for profiles: 'mean' or 'sum' (default: 'mean')")
    parser.add_argument("--normalize", action="store_true",
                        help="Normalize the final aggregated profile relative to its maximum bin value (max = 1.0)")
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Reading transcript profiles from {args.input}...")
    df = pd.read_csv(args.input)
    
    # Check if necessary columns exist
    if "length" not in df.columns or "transcript_id" not in df.columns:
        print("Error: Input CSV must contain 'transcript_id' and 'length' columns.", file=sys.stderr)
        sys.exit(1)
        
    # Get bin columns (bin_1 to bin_100)
    bin_cols = [col for col in df.columns if col.startswith("bin_")]
    if len(bin_cols) != 100:
        print(f"Error: Expected 100 bin columns (bin_1 ... bin_100), found {len(bin_cols)}.", file=sys.stderr)
        sys.exit(1)
        
    # Parse intervals
    intervals = parse_intervals(args.intervals)
    
    summary_rows = []
    
    if intervals is None:
        # Calculate single global metagene summary for all transcripts
        print("No intervals provided. Calculating a single global metagene summary...")
        n_transcripts = len(df)
        if n_transcripts > 0:
            if args.method == "sum":
                profile = df[bin_cols].sum().values
            else:
                profile = df[bin_cols].mean().values
                
            if args.normalize:
                max_val = max(profile)
                if max_val > 0:
                    profile = [x / max_val for x in profile]
                    
            row = {
                "length_bin": "all",
                "n_samples": 1,
                "mean_n_transcripts": n_transcripts
            }
            for i, val in enumerate(profile):
                row[f"bin_{i+1}"] = val
            summary_rows.append(row)
    else:
        # Group by specified intervals
        print(f"Summarizing by intervals: {args.intervals}")
        for label, filt in intervals:
            # Apply filter
            sub_df = df[df["length"].apply(filt)]
            n_transcripts = len(sub_df)
            
            if n_transcripts > 0:
                if args.method == "sum":
                    profile = sub_df[bin_cols].sum().values
                else:
                    profile = sub_df[bin_cols].mean().values
                    
                if args.normalize:
                    max_val = max(profile)
                    if max_val > 0:
                        profile = [x / max_val for x in profile]
            else:
                profile = np.zeros(100)
                
            row = {
                "length_bin": label,
                "n_samples": 1,
                "mean_n_transcripts": n_transcripts
            }
            for i, val in enumerate(profile):
                row[f"bin_{i+1}"] = val
            summary_rows.append(row)
            
            print(f"  Interval '{label}': {n_transcripts} transcripts summarized ({args.method} aggregation).")
            
    # Write output
    out_df = pd.DataFrame(summary_rows)
    # Ensure correct column ordering
    col_order = ["length_bin", "n_samples", "mean_n_transcripts"] + [f"bin_{i}" for i in range(1, 101)]
    out_df = out_df[col_order]
    
    print(f"Writing summarized metagenes to {args.out}...")
    out_df.to_csv(args.out, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
