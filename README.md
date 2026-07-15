# sammy_r21

## Overview

Nextflow pipeline for scRNA-seq analysis of the Sammy R21 cohort — 15 Cell Ranger
samples, 10 `normal` vs 5 `obese`. The pipeline converts each sample's Cell Ranger
count matrix into a QC-annotated AnnData `.h5ad`; samples are clustered individually
downstream.

The samplesheet (`assets/samplesheet.csv`) has columns `sample,path`. Sample ids are
`<condition>_id_<study_id>` (e.g. `normal_id_1`), and that is split into three `obs`
columns on the resulting object — `sample` (`normal_id_1`), `id` (`1`), and
`condition` (`normal`) — so the design travels with the data.

## Directory structure

```
sammy_r21/
├── main.nf           # Single entry point; dispatches on --step
├── nextflow.config   # Params and the local / oscer profiles
├── modules/          # One .nf process module per step
├── bin/              # Python scripts invoked by the modules (argparse CLIs)
├── assets/           # Samplesheets
├── resources/        # Reference files (PPTX template, marker gene lists)
├── scripts/          # Quarto report decks (not yet written)
└── data/             # Gitignored; raw inputs live on OSCER scratch
```

## Setup

Runs use the `babiddy755/python_spatial:1.2.0` container via Apptainer, so no local
install is needed for the `local` or `oscer` profiles. To run with no profile, or for
ad hoc analysis of the produced h5ads:

```bash
conda env create -f environment.yml
conda activate sammy_r21
```

## Usage

```bash
# Build per-sample h5ads on OSCER (the real target — inputs are on /ourdisk)
nextflow run main.nf --step create_adata -profile oscer --samplesheet assets/samplesheet.csv

# Verify wiring without executing anything
nextflow run main.nf --step create_adata -stub --samplesheet assets/samplesheet.csv
```

Outputs land in `<out_root>/<run_id>/results/<sample>/create_adata/<sample>.h5ad`,
alongside a `create_adata_samplesheet.csv` handoff sheet listing them all as
`sample,path`.
