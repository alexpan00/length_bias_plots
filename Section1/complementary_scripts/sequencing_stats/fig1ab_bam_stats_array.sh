#!/bin/bash
#SBATCH --job-name=bam_stats
#SBATCH --output=logs/stats_%A_%a.out
#SBATCH --error=logs/stats_%A_%a.err
#SBATCH --array=1-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# ------------------------------------------------------------------------------
# SLURM Array Script for BAM Statistics (Multi-column support)
# Usage:
# 1. Create a 3-column fofn (TAB separated): [Path] [Dataset] [Technology]
#    Example: /path/to/B31.bam  NIH  Kinnex
# 2. Update 'N_FILES' in the #SBATCH --array line.
# 3. Submit: sbatch complementary_scripts/sequencing_stats/fig1ab_bam_stats_array.sh
#
# The parsed-stats helper (fig1_bam_mapped_stats.py) lives in the same directory
# as this script. Under Slurm its location is derived from SLURM_SUBMIT_DIR;
# override with the environment variable STATS_SCRIPT_DIR if needed.
# ------------------------------------------------------------------------------


ml samtools
FOFN="bams.fofn"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FOFN)

BAM_PATH=$(echo "$LINE" | awk '{print $1}')
DATASET=$(echo "$LINE" | awk '{print $2}')
TECH=$(echo "$LINE" | awk '{print $3}')

SAMPLE_NAME=$(basename "$BAM_PATH" .bam)
OUT_DIR="read_stats/${DATASET}/${TECH}"

mkdir -p "$OUT_DIR" logs

echo "Processing: $DATASET | $TECH | $SAMPLE_NAME"

MAPPED_SCRIPT="fig1_bam_mapped_stats.py"

if [ ! -f "$MAPPED_SCRIPT" ]; then
    echo "Error: Helper script not found: fig1_bam_mapped_stats.py" >&2
    exit 1
fi

# Run samtools stats (Full read sequence stats)
samtools stats -@ 4 "$BAM_PATH" > "${OUT_DIR}/${SAMPLE_NAME}.stats"

# Run custom mapped stats (Mapped read length excluding soft clipping)
python3 "$MAPPED_SCRIPT" "$BAM_PATH" > "${OUT_DIR}/${SAMPLE_NAME}.mapped_stats"

echo "Finished $SAMPLE_NAME"
