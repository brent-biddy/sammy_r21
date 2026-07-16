# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

**Keep this file small.** It is loaded every session; code comments are not, so this
file earns its place only by carrying what would change what you do *before* you'd have
opened any file. Everything else — how a step works, why a plot looks the way it does,
what a config knob is for — goes in a comment next to that code, where it is reviewed
with the change and rots visibly. **If a note here restates a comment, delete the note.**
This file was 395 lines and cost two edits per correction; that is the failure mode to
avoid.

## What this pipeline does

`sammy_r21` is a Nextflow pipeline for scRNA-seq analysis of the Sammy R21 cohort:
15 Cell Ranger samples, **10 `normal`** vs **5 `obese`** (`assets/samplesheet.csv`),
raw matrices on OSCER scratch. One entry point, `main.nf`, dispatching on `--step` —
see its header for the step list, and each `bin/<step>.py` docstring for what that step
does.

`create_adata` and `cluster` are per-sample (one h5ad in, one out). `qc_report` is the
only fan-in step. `create_adata` and `qc_report` filter nothing; `cluster` is the only
step that filters.

## Settled — do not re-open without new evidence

Pointers only. The reasoning lives in the code, next to the thing it explains.

- **Mito cutoff: 50%, fixed cohort-wide.** Not per-sample or adaptive (miQC, median +
  3×MAD) — the condition confound was checked across all 15 samples and did not hold.
  → `cluster.py` (`MITO_MAX`), `qc_report.qmd` (`MITO_BINS`).
- **Use `doublet_score`, not `predicted_doublet`** — scrublet's auto-threshold is
  unstable across samples here. → `cluster.py`, at the scrublet call. *Open: revisit
  pinning a fixed threshold once all 15 have run.*
- **No concat step**, and none planned until per-sample clustering says merging is
  warranted. The per-sample clustered h5ads are the current endpoint.
- **Sample ids encode the design**, so renaming a sample silently changes its
  `condition`. → `create_adata.py` docstring.
- **No Quarto params machinery** — staging is the input contract, and both notebooks
  glob `*.h5ad` from their own directory. Revisited when `cluster_report` was added
  (the trigger that was written down) and still not warranted. → `qc_report.nf`.

## Commands

`--samplesheet` is always required. Steps do not chain inside Nextflow: every
artifact-producing step publishes a `<step>_samplesheet.csv` handoff sheet into
`outdir`, and you point the next step at it.

```bash
# 1. build per-sample h5ads
nextflow run main.nf --step create_adata -profile oscer --samplesheet assets/samplesheet.csv

# 2 & 3. qc_report and cluster are independent consumers of create_adata's handoff
# sheet and can run in either order — but read the report first: it is where the mito
# cutoff cluster applies was chosen.
nextflow run main.nf --step qc_report -profile oscer \
    --samplesheet /scratch/$USER/sammy_r21_out/<run_id>/results/create_adata_samplesheet.csv

nextflow run main.nf --step cluster -profile oscer --resolutions '0.4 0.8' \
    --samplesheet /scratch/$USER/sammy_r21_out/<run_id>/results/create_adata_samplesheet.csv

# 4. cluster_report — reads cluster's handoff sheet, not create_adata's
nextflow run main.nf --step cluster_report -profile oscer \
    --samplesheet /scratch/$USER/sammy_r21_out/<run_id>/results/cluster_samplesheet.csv

# wiring check and config parse check. Note `nextflow config .` does NOT compile
# main.nf — only a run or -stub run catches a syntax error there.
nextflow run main.nf --step create_adata -stub --samplesheet assets/samplesheet.csv
nextflow config .
```

`-profile oscer` is the real target; `-profile local` runs the 2-sample
`assets/test_samplesheet.csv` against a local copy of the raw data in `~/R21/`. Both
write outside the repo and are documented in `nextflow.config`. With no profile you
need the `environment.yml` conda env activated.

## Adding a step

1. `bin/<name>.py` with an `argparse` CLI (`parse_args()`), committed executable
   (`chmod +x` — `bin/` scripts run from `PATH`, and a non-executable one fails with a
   bare exit 126).
2. A process in `modules/<name>.nf` and a `--step` branch in `main.nf`. Copy an
   existing module's shape — the `script:` preamble and the publish-dir helper both
   carry constraints, and their comments say which.
3. Verify with `nextflow config .` and a `-stub` run.

## Conventions

- Work on feature branches and merge via PR (squash). Never commit directly to `main`.
- 4-space indentation in `.nf` files. Process names `UPPER_SNAKE_CASE`; params,
  variables, and CSV headers `snake_case`.
- Always `script:` blocks, never `exec:` — processes must run through SLURM.
- File-level header comments, docstrings on helpers, and WHY comments for non-obvious
  decisions. Annotate channel shape after non-obvious transformations.
- `data/`, `results*/`, `work/`, and `docs/` are gitignored. `resources/` is tracked
  (the OUHSC PPTX template).
