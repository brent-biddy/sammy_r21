# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**What belongs here:** decisions that span files, knowledge about code that no longer
exists, and cohort facts. **What does not:** how a given file works — that goes in a
comment *in that file*, where it is reviewed with the change and rots visibly. If a note
here restates a comment, delete the note, not the comment. This file was trimmed from
395 lines to roughly this size once the duplication started costing two edits per
correction.

## What this pipeline does

`sammy_r21` is a Nextflow pipeline for scRNA-seq analysis of the Sammy R21 cohort:
15 Cell Ranger samples, **10 `normal`** vs **5 `obese`** (`assets/samplesheet.csv`),
with the raw matrices living on OSCER scratch. All steps run through a single entry
point, `main.nf`, selected with `--step`:

- `create_adata` — Cell Ranger `filtered_feature_bc_matrix` → a sample-level
  `<sample>.h5ad` carrying raw counts, `sample`/`id`/`condition` in `obs`, and QC
  metrics. **Filters nothing** — it is the raw artifact.
- `qc_report` — renders `notebooks/qc_report.qmd` to `qc_report.pptx` over the whole
  cohort. Read-only: filters nothing, writes no h5ad, plots off the **raw** matrix.
  This is the step to look at **before** choosing thresholds. The only fan-in step.
- `cluster` — the scanpy pass over one `create_adata` h5ad (mito filter, scrublet,
  normalize/log1p, HVG, PCA, neighbours, UMAP, Leiden sweep). One sample in, one out.
  The only step that **filters**, and mito is the only threshold it applies.

Each script's own docstring and comments carry the detail. The rest of this file is
what those comments cannot say.

## Settled decisions

Re-opening these costs time and has already been paid for once each.

- **Mito cutoff: 50%, fixed, all samples.** Not per-sample, not adaptive (miQC,
  median + 3×MAD) — the condition confound was checked across all 15 samples and did
  not hold, so an adaptive threshold would buy nothing and cost comparability.
  **Do not re-raise without new evidence.** Full rationale in `qc_report.qmd`'s
  `MITO_BINS` comment; `cluster.py`'s `MITO_MAX` applies it.
- **Use `doublet_score`, not `predicted_doublet`.** Scrublet's auto-threshold is
  unstable across samples here (0.68 vs 0.30 on near-identical distributions). The
  score is comparable; the thresholded call is not. Rationale at `cluster.py`'s
  scrublet call. **Open:** revisit pinning a fixed threshold once all 15 have run and
  the score distributions can be compared — the same evidence path the mito call took.
- **No concat step, and samples are clustered individually first.** One is not planned
  until per-sample clustering says merging is warranted. Both `create_adata` and
  `cluster` are strictly one-sample-in / one-h5ad-out; the per-sample clustered h5ads
  are the current endpoint.
- **Sample ids encode the design** (`<condition>_id_<study_id>`), so the samplesheet
  stays `sample,path` and **renaming a sample silently changes its `condition`**. The
  plain three-way `split("_")` in `create_adata.py` is the whole validation, and
  deliberately so — a helper or regex to catch narrower typos was considered and
  rejected. If ids ever stop encoding the design, move `condition` to an explicit
  samplesheet column rather than making the parse cleverer.
- **No Quarto params machinery** (`notebook_registry.json`, `quarto_params.nf`). That
  exists in `xenium_nb` to validate params across many notebooks and is not worth it
  for one report — staging *is* the input contract. **A second notebook is the trigger
  to revisit this**, and one is likely: a before/after plotting doc. Note it would need
  both the raw and clustered h5ads, which share the name `<sample>.h5ad` and so cannot
  be staged flat — see the `stageAs` note in `cluster.nf`.

## Knowledge about code that no longer exists

The only record of work that was done and thrown away. Nothing to grep for — that is
the point of writing it here.

- **Mito-retention curves were built and then removed** (cells discarded vs cutoff, one
  line per sample, plus a per-sample grid). Do not rebuild them speculatively: they got
  the one call they were used for wrong, because they are marginal in mito while the
  question is joint. The scatter sees the same 60–70 cells at 100–300 genes — debris
  regardless of where the mito density thins. **When a marginal and a joint view
  disagree, the joint view wins.** If the question ever returns: the curve is the
  complement of the mito ECDF, so **its slope at a cutoff is the density of cells
  sitting there** — a flat stretch is a gap between populations, a steep one is a dense
  population being sliced.
- **Hidden subplots still reserve their space under `tight_layout`.** `set_visible(False)`
  hides an axes; it does not reclaim its grid cell. The removed small-multiples grid had
  a fixed column count, so with only the 2 test samples staged it squeezed its panels
  into the left two fifths of the slide and stranded `supxlabel` to the right of them.
  Any future grid wants `ncols = min(cap, len(sample_order))`.
- **Fixes here have twice outlived their cause** (a tick-count cap, a colorbar alpha
  workaround). When a panel changes, check whether its workarounds still apply.

## Reading a local render

Only `normal_id_20` and `obese_id_23` have h5ads built from the local `~/R21` copy, so
a `local`-profile render errors on every other sample's section — 2 content slides
each. That is `execute: error: true` working: each missing section fails on its own
slide instead of taking the deck down. The cost is that **a genuinely broken cell also
yields an error slide and exit 0** — read the slides, not the exit status, and check
the error count matches the number of unstaged samples exactly.

The per-sample sections are hardcoded (one `# <sample>` block each) while the summary
table iterates `sample_order`; only the sections need adding. Regenerate them from
`assets/samplesheet.csv` rather than typing ids by hand — the rationale for hardcoding
is in the qmd.

## Commands

### Run a step
`--samplesheet` is always required; columns vary by step (`create_adata` takes
`sample,path`).

```bash
# 1. build per-sample h5ads
nextflow run main.nf --step create_adata -profile oscer --samplesheet assets/samplesheet.csv

# 2. QC report over the cohort — point it at create_adata's handoff sheet
nextflow run main.nf --step qc_report -profile oscer \
    --samplesheet /scratch/$USER/sammy_r21_out/<run_id>/results/create_adata_samplesheet.csv

# 3. cluster each sample — same create_adata handoff sheet as the report reads
nextflow run main.nf --step cluster -profile oscer \
    --samplesheet /scratch/$USER/sammy_r21_out/<run_id>/results/create_adata_samplesheet.csv

# ...optionally overriding the Leiden sweep (quote it — it is one param, many values)
nextflow run main.nf --step cluster -profile oscer --resolutions '0.4 0.8' \
    --samplesheet /scratch/$USER/sammy_r21_out/<run_id>/results/create_adata_samplesheet.csv
```

`qc_report` and `cluster` are independent consumers of `create_adata` and read the same
handoff sheet, so they can run in either order — but read the report first: it is where
the mito cutoff `cluster` applies was chosen.

Every artifact-producing step publishes a handoff samplesheet into `outdir`
(`<step>_samplesheet.csv`) listing its outputs as `sample,path`, so the next step's
`--samplesheet` can be pointed straight at it instead of hand-building a CSV.

### Stub run (no script execution — verifies workflow wiring)
```bash
nextflow run main.nf --step create_adata -stub --samplesheet assets/samplesheet.csv
```

### Config parse check
```bash
nextflow config .
```

### Profiles
Defined in `nextflow.config`:

| Profile | Executor | Container |
|---------|----------|-----------|
| (none)  | local, no container | requires an activated conda env (`environment.yml`) |
| `local` | local, Apptainer | `babiddy755/python_spatial:1.2.0`, 8 CPUs, 16 GB |
| `oscer` | SLURM on OSCER HPC, Apptainer | same image, 16 CPUs, memory retries 48→96→144 GB |

`oscer` is the real target — the Cell Ranger matrices in `assets/samplesheet.csv` are
on `/ourdisk/hpc/lilab/babiddy/dont_archive/sammy_r21/`. The container is reused from
`xenium_nb`; it carries the scanpy + session-info stack, plus `leidenalg`/`igraph` and
the `scikit-image` scrublet needs. It does **not** carry `tabulate` or `scikit-misc`.

A **local copy of all 15 raw samples** lives in `~/R21/`, laid out exactly like the
OSCER tree (same sample dirs, same `filtered_feature_bc_matrix_*` subdirs — only the
root differs), so a test sheet is the OSCER one with the prefix rewritten.
`assets/test_samplesheet.csv` is that: two samples, one per condition
(`normal_id_20`, `obese_id_23` — the smallest of each), and it is what the `local`
profile defaults to:

```bash
nextflow run main.nf --step create_adata -profile local   # uses assets/test_samplesheet.csv
```

**Run directories.** The `local` and `oscer` profiles set their own `workDir` and
`outdir` so nothing lands in the repo. Each run gets one self-contained directory,
`<out_root>/<run_id>/{work,results}`, so a whole run is a single unit to size
(`du -sh`) or prune (`rm -rf`). The shared Apptainer cache is a sibling of the run
dirs (`<out_root>/apptainer_cache`), never nested under a `run_id`, so it survives
across runs:

- `local` → `~/sammy_r21_out/<run_id>/{work,results}`
- `oscer` → `/scratch/$USER/sammy_r21_out/<run_id>/{work,results}`

Keeping `work` and `results` under the same root also keeps them on one filesystem,
which the modules' hardlink publishing (`mode: 'link'`) relies on to avoid a second
full copy of each matrix.

Both profiles set `cleanup = true`, so the work dir is deleted once the run **completes
successfully** — leaving only the (hardlinked) results. A **failed** run keeps its work
dir, so resume-after-failure still works. Because `run_id` defaults to a fresh
timestamp, `-resume` across separate launches only works if you pin it with
`--run_id <name>`, and must be run from the same launch directory.

## Architecture

### Single entry point, one workflow per step
`main.nf` dispatches on `--step` to a named workflow, which reads a samplesheet, builds
a channel of tuples, and pipes it into a single process. There is no chaining between
steps inside Nextflow — to run steps in sequence, point the next step's `--samplesheet`
at the prior step's published handoff samplesheet. Each process emits a
`samplesheet_row` output whose published path comes from a per-module helper that also
drives that module's `publishDir` — so the convention is single-sourced in the module
and `main.nf` just `.map { it.text }` + `collectFile`s the rows (the `.text` read makes
`collectFile`'s `sort` deterministic). The row fragment is kept out of the publish dir
via `publishDir`'s `saveAs`.

### Scripts (`bin/`)
Each step runs a plain Python script with an `argparse` CLI (`bin/<step>.py`), invoked
directly from its module's `script:` block — so they must be committed executable
(`100755`), or the process fails with exit 126. `bin/timer.py` is a shared helper
providing the `timer` context manager and `timing_summary`, which every script uses to
emit `<sample>_timing.tsv`.

### Process conventions
- Always use `script:` blocks, never `exec:` — processes must run through SLURM.
- Every process script sets `XDG_CACHE_HOME=$PWD/.cache` and `TMPDIR=$PWD/tmp` to avoid
  writing to a read-only compute-node `/tmp`.
- Keep named input variables; do not inline maps into process call arguments.
- Build command lines with optional arguments using a Groovy list + conditional append.

## Adding a step

1. Create `bin/<name>.py` with an `argparse` CLI (`parse_args()` function), `chmod +x`.
2. Wire a new process into `modules/<name>.nf` and add a matching `--step` branch in
   `main.nf`, passing args directly.
3. Verify with `nextflow config .` and a `-stub` run.

## Data & resources

- `data/` and `results*/`/`work/` are **gitignored** — raw inputs live on OSCER, and run
  outputs land outside the repo under the profile's `out_root`.
- `resources/` (tracked): reference files — the OUHSC PPTX template, marker gene lists.
- `scripts/` — reserved for the eventual Quarto report deck (the oir-analysis pattern:
  `quarto render scripts/NN_name.qmd --execute-dir . --output-dir ..`). Empty for now.
- `docs/` is **gitignored** — documentation is kept locally, not committed.

## Conventions

- Work on feature branches and merge via PR (squash). Never commit directly to `main`.
- 4-space indentation in `.nf` files.
- Process names in `UPPER_SNAKE_CASE`; params, variables, CSV headers in `snake_case`.
- Add file-level header comments, docstrings on helper functions, section markers, and
  WHY comments for non-obvious decisions.
- Annotate channel shape after non-obvious transformations so the tuple structure is
  visible without tracing back through the chain.
