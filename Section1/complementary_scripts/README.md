# Complementary Scripts

Supplementary scripts used to generate additional data that are not part of the
main Snakemake pipeline.

## Repository layout

```
complementary_scripts/
├── README.md
├── metagene/
│   ├── calculate_transcript_coverage.py  # Per-transcript coverage from a BAM (bin_1..bin_100)
│   ├── summarize_metagene_coverage.py    # Aggregate profiles into length-binned metagene curves
│   ├── average_metagene_summaries.py     # Average multiple metagene summaries into one file
│   ├── average_all_folders.sh            # Batch-average all metagene folders under a base dir
│   └── submit_metagene_array.sh          # Slurm array: sort BAM -> coverage -> summarize
└── sequencing_stats/
    ├── fig1ab_bam_stats_array.sh         # Slurm array for per-sample BAM statistics
    ├── fig1_bam_mapped_stats.py          # Custom mapped-length stats (clipping excluded)
    └── fig1ab_parse_bam_stats.py         # Combine .stats files into a summary CSV
```

The Python helper scripts live in the same directory as the shell script that
calls them. Every shell script resolves these helpers relative to its own
location, so the repository can be cloned anywhere and run as-is.

## Dependencies

- Python 3 with `pandas` and `numpy`
- `samtools`
- A Slurm cluster (only for the two `*_array.sh` scripts)

## Metagene pipeline

The metagene workflow produces a per-transcript coverage profile and then
length-binned summaries, which feed the metagene coverage panels.

### 1. Per-transcript coverage

```bash
python3 calculate_transcript_coverage.py --bam <input.bam> --out <output.csv> [--normalize]
```

### 2. Length-binned metagene summaries

```bash
python3 summarize_metagene_coverage.py --input <profiles.csv> --out <summary.csv> \
    --intervals "0-1000,1000-2000,2000-3000,3000-6000,>=6000"
```

### 3. Average summaries across folders (batch)

```bash
./average_all_folders.sh [base_dir] [output_dir]
```

Defaults to the repository layout `<repo>/data/metagene_coverage/metagenes_transcripts`
and writes averages to `<repo>/data/metagene_coverage`.

### 4. Slurm array (per-sample pipeline)

```bash
sbatch --array=1-N complementary_scripts/metagene/submit_metagene_array.sh <input.fofn> <output_dir>
```

Each line of `input.fofn` contains a path to a BAM file (optionally with
tab-separated metadata, which is ignored). The helper scripts are located
relative to `SLURM_SUBMIT_DIR`; override with the environment variable
`METAGENE_SCRIPT_DIR` if they are not in the standard repo location.

## Sequencing statistics

### 1. Compute per-sample stats (Slurm array)

Create a 3-column tab-separated `bams.fofn`:

```
/path/to/B31.bam  NIH  Kinnex
/path/to/B32.bam  NIH  Kinnex
```

Update `N_FILES` in the `#SBATCH --array=1-N_FILES%10` line of
`fig1ab_bam_stats_array.sh`, then submit it (e.g. with `sbatch
complementary_scripts/sequencing_stats/fig1ab_bam_stats_array.sh`).
Output is written to `read_stats/<Dataset>/<Technology>/`.

The helper `fig1_bam_mapped_stats.py` is located relative to `SLURM_SUBMIT_DIR`;
override with the environment variable `STATS_SCRIPT_DIR` if needed.

### 2. Parse into a summary CSV

```bash
python3 fig1ab_parse_bam_stats.py
```

Reads all `.stats` files under `read_stats/` and writes
`analysis_results/sequencing_performance_summary.csv`.
