#!/bin/bash
#SBATCH --job-name=metagene_array
#SBATCH --output=logs/metagene_%A_%a.out
#SBATCH --error=logs/metagene_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# submit_metagene_array.sh
#
# Process a single BAM file from a FOFN in parallel using Slurm arrays.
# Sorts the BAM, calculates normalized coverages (saved to output_dir),
# and summarizes them into length-binned metagene averages (saved to output_dir).
#
# Usage:
#   sbatch --array=1-N complementary_scripts/metagene/submit_metagene_array.sh <input.fofn> <output_dir>
#
# Notes:
#   - The helper Python scripts (calculate_transcript_coverage.py and
#     summarize_metagene_coverage.py) live in the same directory as this script.
#   - When run under Slurm, the helper directory is derived from SLURM_SUBMIT_DIR
#     (the directory where sbatch was invoked) plus the repo-relative path above.
#     Override with the environment variable METAGENE_SCRIPT_DIR if needed.

set -e
set -o pipefail

FOFN=$1
OUT_DIR=$2

if [ -z "$FOFN" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: sbatch --array=1-N $0 <input.fofn> <output_dir>" >&2
    exit 1
fi

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a Slurm job array (using sbatch --array=1-N)." >&2
    exit 1
fi

# 1. Extract BAM Path and Info
RAW_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FOFN")
BAM_PATH=$(echo "$RAW_LINE" | cut -f1 | xargs)

if [ -z "$BAM_PATH" ]; then
    echo "Error: No BAM path found on line $SLURM_ARRAY_TASK_ID of $FOFN" >&2
    exit 1
fi

if [ ! -f "$BAM_PATH" ]; then
    echo "Error: BAM file not found: $BAM_PATH" >&2
    exit 1
fi

# Get sample ID from the basename of the BAM file
SAMPLE_ID=$(basename "$BAM_PATH" .bam)
PROFILES_CSV="${OUT_DIR}/${SAMPLE_ID}_profiles.csv"
OUT_PATH="${OUT_DIR}/${SAMPLE_ID}.csv"

echo "=========================================================="
echo "Array Task ID   : $SLURM_ARRAY_TASK_ID"
echo "Sample ID       : $SAMPLE_ID"
echo "Input BAM       : $BAM_PATH"
echo "Profiles Output : $PROFILES_CSV"
echo "Summary Output  : $OUT_PATH"
echo "=========================================================="


CALC_SCRIPT="calculate_transcript_coverage.py"
SUM_SCRIPT="summarize_metagene_coverage.py"

# Validate scripts exist
if [ ! -f "$CALC_SCRIPT" ] || [ ! -f "$SUM_SCRIPT" ]; then
    echo "Error: Helper scripts not found" >&2
    exit 1
fi

mkdir -p "$OUT_DIR" logs

# Create a local temporary directory for the intermediate sorted BAM
TEMP_DIR=$(mktemp -d -t metagene_task_tmp.XXXXXXXXXX)
cleanup() {
    rm -rf "$TEMP_DIR"
}
trap cleanup EXIT

SORTED_BAM="${TEMP_DIR}/sorted.bam"

# 3. Sort BAM File
echo "Sorting BAM to $SORTED_BAM..."
samtools sort -@ "$SLURM_CPUS_PER_TASK" -T "${TEMP_DIR}/sort_tmp" -o "$SORTED_BAM" "$BAM_PATH"

# 4. Calculate Coverage Profiles
METHOD=${METAGENE_METHOD:-"mean"}
if [ "$METHOD" = "sum" ]; then
    echo "Computing raw transcript coverage (sum method)..."
    python3 "$CALC_SCRIPT" --bam "$SORTED_BAM" --out "$PROFILES_CSV" --min-coverage 10
else
    echo "Computing normalized transcript coverage (mean/default method)..."
    python3 "$CALC_SCRIPT" --bam "$SORTED_BAM" --out "$PROFILES_CSV" --normalize --min-coverage 10
fi

# 5. Summarize Coverage per Length Bin
echo "Summarizing metagene coverage by length bins using method: $METHOD..."
INTERVALS="0-1000,1000-2000,2000-3000,3000-6000,>=6000"

if [ "$METHOD" = "sum" ]; then
    python3 "$SUM_SCRIPT" \
        --input "$PROFILES_CSV" \
        --out "$OUT_PATH" \
        --intervals "$INTERVALS" \
        --method sum \
        --normalize
else
    python3 "$SUM_SCRIPT" \
        --input "$PROFILES_CSV" \
        --out "$OUT_PATH" \
        --intervals "$INTERVALS" \
        --method mean
fi

echo "Task completed successfully for $SAMPLE_ID"
