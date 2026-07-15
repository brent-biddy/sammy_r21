# sammy_r21

<!-- TODO: describe the biological question, the samples, and the deliverable. -->

Single-cell RNA-seq analysis for the Sammy R21.

## Environment

Conda env `sammy_r21` (`environment.yml`). Activate before running anything:

```bash
conda activate sammy_r21
```

## Pipeline

Linear; **run everything from the repo root** — all paths are relative to it
(the scripts don't anchor themselves, so the working directory matters).
Processed objects live in `data/processed/`.

1. **`scripts/01_create_anndata.py`** — build the AnnData from the raw
   CellRanger output in `data/raw/`; attach per-cell sample metadata, set obs
   color palettes, compute QC metrics. Output: `data/processed/adata.h5ad`.

<!-- Subsequent steps follow the oir-analysis pattern — add them here as they're
     written (e.g. 02_cluster.py, 03_annotate.py, 04_*.qmd for the deck). Quarto
     decks render with:
     quarto render scripts/NN_name.qmd --execute-dir . --output-dir ..
     `--execute-dir .` runs the deck's cells from the repo root so its paths
     match the scripts; `--output-dir ..` writes the PPTX to the repo root. -->

## Data & resources

- `data/` is **gitignored** in full — nothing under it is tracked. The pipeline
  scripts create `data/processed/` as needed; `data/raw/` holds the local
  instrument output and is populated by hand.
- `resources/` (tracked): reference files — marker gene lists, the OUHSC PPTX
  template.
- `docs/` is **gitignored** — documentation is kept locally, not committed.

## Conventions

- Work on feature branches and merge via PR (squash). Never commit directly to
  `main`.
- Rendering a deck to PPTX is the deliverable — don't convert or screenshot the
  PPTX.
- `*.pptx` and `*.quarto_ipynb` are gitignored (the template PPTX is the
  exception).
