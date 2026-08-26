#!/bin/bash
#SBATCH --job-name=sample_read_lengths
#SBATCH --output=logs/sample_read_lengths_%A_%a.out
#SBATCH --error=logs/sample_read_lengths_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00

# ------------------------------------------------------------------------------
# SLURM Array Script: sample 1M reads per (dataset, technology) combination.
#
# Usage:
#   1. Create a 3-column (TAB separated) fofn: [Path] [Dataset] [Technology]
#        /path/to/sample.bam   NIH      ONT_cDNA
#        /path/to/seq.bam      LongBench Kinnex
#      - The first column points to a BAM/CRAM file.
#      - Provide EXACTLY ONE sample/BAM per (dataset, technology) combination.
#      - Technology uses underscores for multi-word names (e.g. ONT_cDNA,
#        ONT_dRNA, IsoSeq, Kinnex) to keep filenames clean.
#   2. Point FOFN below at your input, set N_FILES on the #SBATCH --array line
#      to the number of lines (rows) in the fofn.
#   3. Submit: sbatch helper_scripts/sample_read_lengths_array.sh [fofn] [output_dir]
#
# Output (one file per combination):
#     <output_dir>/<Dataset>_<Technology>_read_lengths.txt.gz
#   A gzip-compressed file with ~1,000,000 lines, one per sampled read,
#   each line being the read length in bp.
# ------------------------------------------------------------------------------

set -e
set -o pipefail

ml samtools

FOFN="${1:?Usage: sbatch --array=1-N sample_read_lengths_array.sh <input.fofn> <output_dir>}"
OUT_DIR="${2:?Usage: sbatch --array=1-N sample_read_lengths_array.sh <input.fofn> <output_dir>}"

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a Slurm job array (sbatch --array=1-N)." >&2
    exit 1
fi

mkdir -p "$OUT_DIR" logs

# Read this array task's row from the fofn
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FOFN")
BAM_PATH=$(echo "$LINE" | cut -f1 | xargs)
DATASET=$(echo "$LINE"   | cut -f2 | xargs)
TECH=$(echo "$LINE"      | cut -f3 | xargs)

if [ -z "$BAM_PATH" ] || [ -z "$DATASET" ] || [ -z "$TECH" ]; then
    echo "Error: Malformed row on line $SLURM_ARRAY_TASK_ID of $FOFN" >&2
    exit 1
fi

if [ ! -f "$BAM_PATH" ]; then
    echo "Error: BAM/CRAM file not found: $BAM_PATH" >&2
    exit 1
fi

echo "=========================================================="
echo "Array Task ID : $SLURM_ARRAY_TASK_ID"
echo "Input BAM     : $BAM_PATH"
echo "Dataset       : $DATASET"
echo "Technology    : $TECH"
echo "=========================================================="

OUT_FILE="${OUT_DIR}/${DATASET}_${TECH}_read_lengths.txt.gz"

# Extract read lengths and sample 1M of them.
#   - Flags filtered out: 0x100 (secondary), 0x800 (supplementary), 0x4 (unmapped)
#     -> FLAG = 0x904, keeping primary mapped reads only.
#   - Read length is taken as the length of the SEQ field (column 10).
#   - `shuf -n 1000000` performs streaming reservoir sampling, so memory stays
#     bounded even for very large BAMs.
samtools view -F 0x904 "$BAM_PATH" | \
    awk 'BEGIN{FS="\t"} {n=length($10); if (n>0) print n}' | \
    shuf -n 1000000 | \
    gzip -c > "$OUT_FILE"

N_READS=$(zcat "$OUT_FILE" | wc -l)
echo "Done: $OUT_FILE ($N_READS reads)"
