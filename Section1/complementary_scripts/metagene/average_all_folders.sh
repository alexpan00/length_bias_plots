#!/bin/bash
# average_all_folders.sh
#
# Recursively traverses subdirectories under a base directory, finds all metagene
# summary CSV files in each leaf folder, and runs average_metagene_summaries.py
# to generate consolidated average files.
#
# Usage:
#   ./average_all_folders.sh [base_dir] [output_dir]
#
# Defaults:
#   base_dir:   <repo>/data/metagene_coverage/metagenes_transcripts
#   output_dir: <repo>/data/metagene_coverage

set -e

# Locate the helper scripts relative to this script's own directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default to the repository's expected data layout (two levels up from metagene/)
if [ "$#" -eq 0 ]; then
    DEFAULT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    BASE_DIR="${DEFAULT_ROOT}/data/metagene_coverage/metagenes_transcripts"
    OUT_DIR="${DEFAULT_ROOT}/data/metagene_coverage"
else
    BASE_DIR=$1
    OUT_DIR=$2
fi

# Locate Python averaging script (co-located with this script)
AVERAGE_SCRIPT="${SCRIPT_DIR}/average_metagene_summaries.py"

if [ ! -f "$AVERAGE_SCRIPT" ]; then
    echo "Error: Required averaging script not found: average_metagene_summaries.py" >&2
    exit 1
fi

if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Base directory does not exist: $BASE_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "=========================================================="
echo "Batch Averaging Metagene Folders"
echo "Base Directory  : $BASE_DIR"
echo "Output Directory: $OUT_DIR"
echo "=========================================================="

# Find all subdirectories under BASE_DIR
find "$BASE_DIR" -type d | while read -r dir; do
    # Check if there are any CSV files in this directory
    # Using nullglob to avoid literal *.csv expansion when no matches exist
    shopt -s nullglob
    all_files=( "$dir"/*.csv )
    shopt -u nullglob
    
    # Process only if directory contains at least one CSV file
    if [ ${#all_files[@]} -gt 0 ]; then
        # Filter out intermediate _profiles.csv files, leaving only the summary CSVs
        csv_files=()
        for f in "${all_files[@]}"; do
            if [[ "$f" != *_profiles.csv ]]; then
                csv_files+=( "$f" )
            fi
        done
        
        # Calculate relative path from BASE_DIR to name the output file
        # E.g. metagenes/nih/isoseq -> nih_isoseq.csv
        rel_path="${dir#$BASE_DIR/}"
        
        # Avoid processing BASE_DIR itself if it contains CSVs
        if [ "$rel_path" != "$dir" ] && [ -n "$rel_path" ] && [ ${#csv_files[@]} -gt 0 ]; then
            # Replace slashes with underscores for a clean output filename
            out_name=$(echo "$rel_path" | tr '/' '_')
            out_file="${OUT_DIR}/${out_name}.csv"
            
            echo "Processing folder: $rel_path (${#csv_files[@]} files)"
            python3 "$AVERAGE_SCRIPT" --inputs "${csv_files[@]}" --out "$out_file"
            echo "--------------------------------------------------------"
        fi
    fi
done

echo "Batch processing completed successfully!"
