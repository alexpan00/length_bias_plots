# Systematic evaluation of transcript normalization methods and quantification pipelines across long-read sequencing platforms

Code for plots in **"Systematic evaluation of transcript normalization methods and quantification pipelines across long-read sequencing platforms"**.

This repository contains the scripts, manifests, PNG figures, and small analysis inputs needed to reproduce every figure and table in the paper. Large input data is **not** committed to this repository — it is hosted externally and must be fetched before running (see [Data](#data)).

## Repository layout

The analysis is organized in four sections, matching the results:

| Folder | Contents |
|---|---|
| `Section1/` | Figure 1 + supplementary: sequencing throughput, read lengths, metagene coverage, transcript discovery, SQANTI3 classifications, UJC overlap (upset) plots |
| `Section2/` | Figure 2 + supplementary: length-expression relationship under different normalizations, coefficient-of-variation (CV) analyses |
| `Section3/` | Figure 3 + supplementary: SIRV/ERCC spike-in absolute quantification (RMSE) and sensitivity |
| `Section4/` | Figure 4 + supplementary: orthogonal validation against Illumina short reads, linear mixed model (LMM) and coefficients heatmap |
| `manuscript_tables/` | Regenerated statistics tables referenced by the manuscript draft |
| `generate_manuscript_tables.py` | Reproduces all manuscript statistics (Sections 3.1–3.4) into `manuscript_tables/` |
| `export_supplementary_table_s4.Rmd` | Fits the LMM ANOVAs and exports Supplementary Table S4 |


Each `SectionN/` folder contains `.Rmd` (or `.py`/`.sh`) scripts that read data from `SectionN/data/` and write figures/tables.

## Data

### Small data (committed)
All `.fofn` manifests, small summary CSVs, and the light analysis inputs are included in the repository.

### Large data (externally hosted — fetch before running)
The following large analysis inputs are **excluded** from Git (see `.gitignore`) because they exceed GitHub's practical limits. They are hosted on Zenodo ([10.5281/zenodo.22109031](https://doi.org/10.5281/zenodo.22109031)) as two compressed archives and must be restored into the paths shown before the corresponding scripts will run:

| Archive (from Zenodo) | Zip size | Restore to | Unzipped | Used by |
|---|---|---|---|---|
| `metagenes_transcripts.zip` | ~1.2 GB | `Section1/data/metagene_coverage/` | ~4.5 GB | `Section1/supplementary_metagenes.Rmd` |
| `UJC.zip` | ~20 MB | `Section1/data/` | ~405 MB (UJC `.fofn` manifests + `*_UJC_summary.tsv`) | `Section1/upset.Rmd` |

**To reproduce:**

```bash
# 1. Download the two archives from Zenodo (doi.org/10.5281/zenodo.22109031)
# 2. Extract them preserving their internal folder structure:
cd Section1/data/metagene_coverage && unzip metagenes_transcripts.zip   # creates metagenes_transcripts/
cd Section1/data                      && unzip UJC.zip                 # creates the UJC/*.fofn and *_UJC_summary.tsv files
```

The archives preserve the original internal structure, so extracting them at the locations above populates the exact paths the scripts expect. Keep the rest of the repository untouched.

## Reproducing the results

Requirements:
- R ≥ 4.1 with `tidyverse`, `dplyr`, `lme4`, `lmerTest`, `emmeans`, `patchwork`, `svglite`, `ragg`, `readr`
- Python ≥ 3.8 with `pandas`, `numpy` (for `generate_manuscript_tables.py`)

Recommended run order:

```bash
# 1. Restore large data (see Data section) into Section1/data/

# 2. Regenerate manuscript statistics tables
python generate_manuscript_tables.py

# 3. Export Supplementary Table S4 (LMM ANOVAs)
#    Knit export_supplementary_table_s4.Rmd

# 4. Render each figure Rmd (per-section)
#    e.g. knit Section1/fig1.Rmd, Section3/fig3_sensitivity.Rmd, etc.
```

Each `.Rmd` knits to its `SectionN/plots/` folder and prints named figure files.

## Figures

Only **PNG** versions of the figures are versioned to keep the repository lightweight. `*.tiff`, `*.svg`, and `*.pdf` plot outputs are generated on rendering but ignored by Git (see `.gitignore`).


