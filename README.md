# sammy_r21

## Overview

<!-- TODO: one-paragraph description — the biological question, the samples, and
     the deliverable. -->

Single-cell RNA-seq analysis for the Sammy R21. A short pipeline builds and
clusters an AnnData object from local CellRanger output, annotates cell types,
and renders a presentation deck of the results.

## Directory structure

```
sammy_r21/
├── data/             # Gitignored; not tracked
│   ├── raw/          # Original, immutable input data
│   └── processed/    # Cleaned and transformed data
├── resources/        # Reference files (marker gene lists, PPTX template, etc.)
└── scripts/          # Numbered Python scripts and Quarto documents for sequential analysis steps
```

## Setup

```bash
conda env create -f environment.yml
conda activate sammy_r21
```

Copy or symlink the raw CellRanger output into `data/raw/` — `data/` is
gitignored, so nothing under it is tracked.

## Usage

Run scripts in order from the repo root:

```bash
python scripts/01_create_anndata.py  # build AnnData object and write to data/processed/
```
