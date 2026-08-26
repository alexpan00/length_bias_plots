#!/usr/bin/env python3
"""
calculate_transcript_coverage.py

Calculates average coverage per bin (100 bins) for each transcript in a BAM file.
Outputs a 102-column CSV file: transcript_id, length, bin_1, ..., bin_100.

Usage:
  python3 calculate_transcript_coverage.py --bam <input.bam> --out <output.csv> [--normalize]
"""

import os
import sys
import math
import argparse
import subprocess
from collections import defaultdict

def get_transcript_lengths(bam_path):
    """Parses BAM header to get transcript IDs and their exact lengths."""
    lengths = {}
    cmd = ["samtools", "view", "-H", bam_path]
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        print("Error: samtools is not installed or not in PATH.", file=sys.stderr)
        sys.exit(1)
        
    for line in proc.stdout:
        if line.startswith("@SQ"):
            parts = line.strip().split("\t")
            sn = None
            ln = None
            for p in parts:
                if p.startswith("SN:"):
                    sn = p[3:]
                elif p.startswith("LN:"):
                    ln = int(p[3:])
            if sn and ln:
                lengths[sn] = ln
                
    proc.wait()
    if proc.returncode != 0:
        err_msg = proc.stderr.read()
        print(f"Error reading BAM header: {err_msg}", file=sys.stderr)
        sys.exit(1)
        
    return lengths

def compute_bin_sizes(length):
    """Computes exact coordinate sizes of the 100 bins for a transcript of a given length."""
    bin_sizes = [0] * 100
    for pos in range(1, length + 1):
        bin_idx = min(99, math.floor((pos - 1) * 100 / length))
        bin_sizes[bin_idx] += 1
    return bin_sizes

def parse_coverage(bam_path, lengths, min_len=100):
    """
    Parses BAM coverage using samtools depth.
    Streams non-zero coverage positions, tracks max depth per transcript,
    and groups them into 100 bins per transcript.
    """
    cmd = ["samtools", "depth", bam_path]
    print(f"Executing: {' '.join(cmd)}")
    
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    coverage_sums = defaultdict(lambda: [0.0] * 100)
    max_depths = defaultdict(float)
    count = 0
    
    for line in proc.stdout:
        count += 1
        if count % 5000000 == 0:
            print(f"  Processed {count} depth lines...")
            
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
            
        tx_id, pos_str, depth_str = parts
        pos = int(pos_str)
        depth = float(depth_str)
        
        if tx_id not in lengths:
            continue
            
        L = lengths[tx_id]
        if L < min_len:
            continue
            
        # Track max depth
        if depth > max_depths[tx_id]:
            max_depths[tx_id] = depth
            
        # Determine the bin index (0 to 99)
        bin_idx = min(99, math.floor((pos - 1) * 100 / L))
        coverage_sums[tx_id][bin_idx] += depth
        
    proc.wait()
    if proc.returncode != 0:
        err_msg = proc.stderr.read()
        print(f"Error running samtools depth: {err_msg}", file=sys.stderr)
        sys.exit(1)
        
    print(f"  Finished streaming depth. Extracted coverage for {len(coverage_sums)} expressed transcripts.")
    return coverage_sums, max_depths

def calculate_profiles(coverage_sums, lengths, max_depths, normalize=False, min_coverage=0.0, min_len=100):
    """
    Calculates average coverage per bin (with optional normalization to max bin value = 1.0).
    Discards transcripts unless at least one position has a coverage of at least min_coverage.
    """
    profiles = {}
    for tx_id, sums in coverage_sums.items():
        L = lengths[tx_id]
        if L < min_len:
            continue
            
        # Discard transcripts if no position has a coverage of at least min_coverage
        if max_depths[tx_id] < min_coverage:
            continue
            
        bin_sizes = compute_bin_sizes(L)
        
        # Calculate average depth per bin
        bin_avgs = [0.0] * 100
        for i in range(100):
            if bin_sizes[i] > 0:
                bin_avgs[i] = sums[i] / bin_sizes[i]
                
        # Optional relative normalization
        if normalize:
            max_avg = max(bin_avgs)
            if max_avg > 0:
                bin_avgs = [x / max_avg for x in bin_avgs]
            else:
                continue # Skip transcripts with zero absolute average coverage
                
        profiles[tx_id] = bin_avgs
        
    return profiles

def main():
    parser = argparse.ArgumentParser(description="Calculate average coverage per bin for each transcript in a BAM file.")
    parser.add_argument("--bam", required=True, help="Path to input transcriptome-aligned BAM file")
    parser.add_argument("--out", required=True, help="Path to output CSV file")
    parser.add_argument("--normalize", action="store_true", 
                        help="Normalize coverage profiles relative to each transcript's maximum bin value")
    parser.add_argument("--min-coverage", type=float, default=0.0,
                        help="Minimum coverage at least at one position to retain a transcript")
    parser.add_argument("--min-len", type=int, default=100,
                        help="Minimum transcript length to calculate coverage (default: 100)")
    args = parser.parse_args()
    
    if not os.path.exists(args.bam):
        print(f"Error: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)
        
    # 1. Load lengths
    lengths = get_transcript_lengths(args.bam)
    print(f"Loaded {len(lengths):,} transcript lengths from BAM header.")
    
    # 2. Parse coverage
    coverage_sums, max_depths = parse_coverage(args.bam, lengths, min_len=args.min_len)
    
    # 3. Calculate profiles
    print("Calculating bin coverages...")
    profiles = calculate_profiles(coverage_sums, lengths, max_depths, normalize=args.normalize, min_coverage=args.min_coverage, min_len=args.min_len)
    
    # 4. Write CSV
    print(f"Writing profiles to {args.out}...")
    with open(args.out, "w") as f:
        bin_headers = ",".join(f"bin_{i}" for i in range(1, 101))
        f.write(f"transcript_id,length,{bin_headers}\n")
        for tx_id, profile in profiles.items():
            profile_str = ",".join(f"{x:.10f}" for x in profile)
            f.write(f"{tx_id},{lengths[tx_id]},{profile_str}\n")
            
    print("Done!")

if __name__ == "__main__":
    main()
