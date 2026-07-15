# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

`sammy_r21` is a Nextflow pipeline for scRNA-seq analysis of the Sammy R21 cohort:
15 Cell Ranger samples, **10 `normal`** vs **5 `obese`** (`assets/samplesheet.csv`),
with the raw matrices living on OSCER scratch. All steps run through a single entry
point, `main.nf`, selected with `--step`:

- `create_adata` — converts a Cell Ranger `filtered_feature_bc_matrix` directory into
  a sample-level `<sample>.h5ad`. Attaches `sample`, `id`, and `condition` to every
  cell in `obs` (all `Categorical`), and annotates per-cell/per-gene QC metrics
  (`total_counts`, `n_genes_by_counts`, `pct_counts_mt`, `percent_top`) but **filters
  nothing** — it produces the raw artifact and leaves thresholds to downstream
  analysis. Mito genes are auto-detected from the gene symbols with a case-insensitive
  `^[Mm][Tt]-` match, so human and mouse both work with no species flag. Ported from
  the `xenium_nb` pipeline, where it was the odd one out (scRNA-seq, not Xenium) and
  nothing consumed its h5ad.
- `qc_report` — renders `notebooks/qc_report.qmd` to `qc_report.pptx`: a cohort summary
  table, then per sample a section-header slide, a counts-vs-genes scatter coloured by
  mito %, and a 3-panel histogram (genes / counts / mito %). Read-only — it filters
  nothing and writes no h5ad. This is the step to look at **before** choosing
  clustering QC thresholds. Everything it plots comes off the **raw** count matrix.

**Sample sections are hardcoded**, one `# <sample>` block per sample, and currently
cover only `normal_id_20` and `obese_id_23` — the two with h5ads built from the local
`~/R21` copy. Add the other 13 when the cohort is built on OSCER. Hardcoding is
deliberate: generating headings dynamically needs `output: asis`, and under it an
inline figure lands between the markdown blocks and leaves the *next* heading not
starting its own line, so pandoc emits a literal `## sample` and drops the image. One
cell per sample means `plt.show()` just works.

### The QC report is a deck, and that constrains it

Output is PPTX via `resources/ouhsc_ppt_template.pptx` (`reference-doc`), the same
deliverable pattern as `oir-analysis`. These all fail quietly rather than loudly:

- **`slide-level: 2` is load-bearing.** It makes `#` a Section Header slide and `##` a
  Title-and-Content slide. Pandoc otherwise infers the level from the headings present,
  so adding a `#` would silently demote every `##` from a slide to a bullet.
- **Tables must be markdown**, built by hand and passed to `display(Markdown(...))`.
  `df.style` is HTML-only and renders as *nothing* in pptx; `df.to_markdown()` needs
  `tabulate`, which the container does not carry. Done right they become native
  PowerPoint tables (`graphicFrame`), not images.
- **Figures are sized to the slide.** The template is 16:9 at 10 × 5.625in, so figures
  are ~9 × 4.3. A tall stacked panel works in a scrolling page and is unreadable on a
  slide.
- **No prose after a table.** Trailing text gets pushed onto a slide of its own, titled
  with its first sentence fragment. The cohort total is a row in the table for exactly
  this reason.
- **`execute: error: true`** keeps the render going when a cell raises, so a section
  whose h5ad was not staged errors on its own slide instead of taking the deck down.
  The cost: a genuinely broken cell yields an error slide and **exit 0**. Read the
  slides, not the exit status.

`quarto render` resolves `reference-doc` relative to the qmd's own directory, so the
module stages the template beside the staged notebook rather than referencing its repo
path.

### Plotting traps already hit here

- **A log axis needs log-spaced bins.** `set_xscale("log")` with linear bins renders
  wildly uneven bar widths and misrepresents the distribution — and it looks plausible.
  The counts panel builds edges with `np.logspace`.
- **A colorbar inherits its mappable's alpha.** A semi-transparent scatter draws the
  ramp blended toward white — washed out and nothing like the viridis it is. Points are
  opaque now, which sidesteps it; drop them below 1 and the colorbar needs
  `cbar.solids.set_alpha(1.0)`.
- **Mito colour is unclipped 0–100 on purpose.** Clipping at 25 made a 25% cell and a
  90% cell the same yellow, collapsing the dead-cell population sitting at 75–90%. The
  cost is a near-uniform purple manifold, which is accepted.
- Fixes here have twice outlived their cause (a tick-count cap, the colorbar alpha
  workaround). When a panel changes, check whether its workarounds still apply.

### The QC report is a fan-in step, and stages its inputs

`qc_report` is the only step that is not per-sample: it `toSortedList()`s every h5ad
into one task. Two consequences worth keeping:

- **Staging is the input contract.** Every h5ad is staged flat into the work dir and
  the notebook globs `*.h5ad` from its own directory. There are no Quarto params, no
  `notebook_registry.json`, and no `quarto_params.nf` — that machinery exists in
  `xenium_nb` to validate params across many notebooks, and is not worth it for one
  report. If a second notebook ever appears, revisit; until then, do not add it.
- **The notebook reads `obs` only**, via `ad.read_h5ad(path, backed="r")`. The cohort's
  X is ~5 GB and none of it is needed to plot QC metrics. Keep it backed.

Sorting (`toSortedList`, not `collect`) makes the staged order — and so the report —
reproducible regardless of upstream task completion order.

### Sample ids carry the design

Sample ids are `<condition>_id_<study_id>` (e.g. `normal_id_1`, `obese_id_23`).
`bin/create_adata.py` gets the three obs columns from a plain
`args.sample.split("_")` in `main()` — no helper, no regex — so the samplesheet stays
`sample,path` and the design is not restated anywhere.

Keep it that simple. The three-way unpack is the validation: an id with the wrong
number of underscore-separated parts raises `ValueError` on its own, before the slow
matrix read, and that is the failure that matters. Do not add a helper or a pattern to
catch narrower typos (`normal_id_` parses to an empty id; `normal_ID_1` is accepted) —
those were considered and judged not worth the machinery.

`id` is stored as a string, not an int — it is a label to group and join on (e.g.
against clinical metadata keyed on StudyID), never a quantity to average.

This makes the sample id load-bearing: **renaming a sample silently changes its
`condition`.** If ids ever stop encoding the design, move `condition` to an explicit
samplesheet column rather than making the parse cleverer.

**Samples are clustered individually first — there is no concat step**, and one is not
planned until per-sample clustering says whether merging is warranted. `create_adata`
is strictly one-sample-in / one-h5ad-out; the per-sample h5ads are the current
endpoint.

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
```

The report lands at `<outdir>/qc_report/qc_report.pptx` (~193 KB with 2 samples), so it
is a single `scp` to view off the cluster.

Every artifact-producing step publishes a handoff samplesheet into `outdir`
(`<step>_samplesheet.csv`) listing its outputs as `sample,path`, so the next step's
`--samplesheet` can be pointed straight at it instead of hand-building a CSV.

### Profiles
Defined in `nextflow.config`:

| Profile | Executor | Container |
|---------|----------|-----------|
| (none)  | local, no container | requires an activated conda env (`environment.yml`) |
| `local` | local, Apptainer | `babiddy755/python_spatial:1.2.0`, 8 CPUs, 16 GB |
| `oscer` | SLURM on OSCER HPC, Apptainer | same image, 16 CPUs, memory retries 48→96→144 GB |

`oscer` is the real target — the Cell Ranger matrices in `assets/samplesheet.csv` are
on `/ourdisk/hpc/lilab/babiddy/dont_archive/sammy_r21/`. The container is reused from
`xenium_nb`; it carries the scanpy + session-info stack `bin/create_adata.py` needs.

A **local copy of all 15 raw samples** lives in `~/R21/`, laid out exactly like the
OSCER tree (same sample dirs, same `filtered_feature_bc_matrix_*` subdirs — only the
root differs), so a test sheet is the OSCER one with the prefix rewritten.
`assets/test_samplesheet.csv` is that: two samples, one per condition
(`normal_id_20`, `obese_id_23` — the smallest of each), and it is what the `local`
profile defaults to:

```bash
nextflow run main.nf --step create_adata -profile local   # uses assets/test_samplesheet.csv
```

Verified end-to-end on this data (2026-07-15): both samples produce h5ads with the
right obs, and the mito regex flags exactly the expected 13 protein-coding genes on
the real reference.

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

### Stub run (no script execution — verifies workflow wiring)
```bash
nextflow run main.nf --step create_adata -stub --samplesheet assets/samplesheet.csv
```

### Config parse check
```bash
nextflow config .
```

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
directly from its module's `script:` block. `bin/timer.py` is a shared helper providing
the `timer` context manager and `timing_summary`, which every script uses to emit
`<sample>_timing.tsv`.

### Process conventions
- Always use `script:` blocks, never `exec:` — processes must run through SLURM.
- Every process script sets `XDG_CACHE_HOME=$PWD/.cache` and `TMPDIR=$PWD/tmp` to avoid
  writing to a read-only compute-node `/tmp`.
- Keep named input variables; do not inline maps into process call arguments.
- Build command lines with optional arguments using a Groovy list + conditional append.

## Adding a step

1. Create `bin/<name>.py` with an `argparse` CLI (`parse_args()` function).
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
