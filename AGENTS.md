# AGENTS.md

Reproduction code for the paper "Systematic evaluation of transcript normalization methods and quantification pipelines across long-read sequencing platforms". This is an R/Python research repo (plots + tables), not a software project — there is no build, lint, typecheck, or test suite to run.

## Layout & entrypoints

- Four `SectionN/` folders, each holding `.Rmd` scripts that read from `SectionN/data/` and write figures/tables into `SectionN/plots/`.
- `generate_manuscript_tables.py` (repo root) computes all manuscript statistics and writes to `manuscript_tables/`.
- `export_supplementary_table_s4.Rmd` (repo root) fits LMM ANOVAs and exports Supplementary Table S4.

## Running code

- **`.Rmd` path resolution is per-folder**: each script reads data via paths relative to its own `SectionN/` directory (e.g. `fig1.Rmd` uses `data/sequencing_performance_summary.csv`). Knit with the Rmd's folder as working directory — do not `setwd` to repo root.
- **Python script is repo-root-relative**: `python generate_manuscript_tables.py` must run from the repo root (it uses `__file__` paths, so it works from anywhere, but outputs to `manuscript_tables/`).

## Data gotchas

- **Large inputs are NOT in git.** Served from Zenodo (doi 10.5281/zenodo.22109031). `Section1/supplementary_metagenes.Rmd` and `Section1/upset.Rmd` will fail until `Section1/data/metagene_coverage/metagenes_transcripts/` and `Section1/data/UJC/*` are extracted in place (see README "Data"). Don't commit these.
- Only **PNG** figures are versioned. `*.tiff`, `*.svg`, `*.pdf` are generated on render and git-ignored (`.gitignore`).
- `.Rhistory`, `.RData`, `.Rproj.user` are git-ignored artifacts; leave them untracked.

## Environment

- R ≥ 4.1: `tidyverse`, `dplyr`, `lme4`, `lmerTest`, `emmeans`, `patchwork`, `svglite`, `ragg`, `readr`, plus `deeptime` (used by `fig1.Rmd`) and `knitr`.
- Python ≥ 3.8: `pandas`, `numpy`.
- Rmd chunks are set globally to suppress warnings/messages; figure/report errors surface as failed chunk paths referencing the per-section data.

## Conventions

- Tool names normalized to `FLAIR`, `IsoQuant`, `Kallisto`, `Oarfish`, `bambu`.
- Normalization names normalized to `CPM`, `TPM`, `CQN`, `EDASeq`, `SMR` (SMR = ratio-based correction, internal names `ratio_correction`/`ratio_counts`). Keep these exact labels when adding tables/plots.
- Dataset display names: `NIH` → `Mouse`, `Neuron` → `iNeuron`, plus `LongBench`. Illumina read counts are halved (fragments, not reads).
