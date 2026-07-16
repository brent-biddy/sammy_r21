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
  nothing and writes no h5ad. This is the step to look at **before** choosing clustering
  QC thresholds. Everything it plots comes off the **raw** count matrix.

  `MITO_CANDIDATE` is a discussion aid, not a decision — nothing downstream reads it,
  and the report filters nothing. It drives the summary table's `dropped ≥N%` column,
  which is the only place discard rates can be compared across samples side by side;
  the per-scatter tallies give the same numbers one sample at a time.

  **Mito-retention curves were built and then removed** (cells discarded vs cutoff, one
  line per sample, plus a per-sample grid). Do not rebuild them speculatively — and note
  they got the one call they were used for wrong, because they are marginal in mito
  while the question is joint (see the mito colour bullet). If the question ever comes
  back: the curve is the complement of the mito ECDF, so **its slope at a cutoff is the
  density of cells sitting there** — a flat stretch is a gap between populations, a
  steep one is a dense population being sliced.

- `cluster` — the standard scanpy pass over one `create_adata` h5ad: mito filter,
  scrublet (flag only), normalize/log1p, HVG, PCA, neighbours, UMAP, and a Leiden
  sweep writing one `leiden_<res>` column per resolution. One sample in, one
  `<sample>.h5ad` out — this is the step that finally **filters**, and the only
  threshold it applies is mito. See "The cluster step applies the cohort's one
  threshold" below; the short version is that `doublet_score` is trustworthy and
  `predicted_doublet` is not.

**Sample sections are hardcoded**, one `# <sample>` block per sample, and all 15 are now
present in samplesheet order (condition, then ascending id — the same order
`sample_order` produces). Hardcoding is deliberate: generating headings dynamically
needs `output: asis`, and under it an inline figure lands between the markdown blocks
and leaves the *next* heading not starting its own line, so pandoc emits a literal
`## sample` and drops the image. One cell per sample means `plt.show()` just works.
Regenerate them from `assets/samplesheet.csv` rather than typing ids by hand.

**A `local`-profile render now yields 26 error slides**, because only `normal_id_20` and
`obese_id_23` have h5ads built from the local `~/R21` copy and the other 13 sections
have nothing to plot. That is `execute: error: true` working — each missing section
fails on its own slide instead of taking the deck down — but it does mean the local test
render is no longer a clean deck to read. Check the two real samples' figures and that
the error count is exactly 26 (13 samples × 2 content slides); anything else is a real
break.

The **summary table is not hardcoded** — it iterates `sample_order`, so it picks up the
other 13 with no edit. Only the per-sample sections need adding.

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
- **Hidden subplots still reserve their space under `tight_layout`.** `set_visible(False)`
  hides an axes; it does not reclaim its grid cell. A small-multiples grid with a fixed
  column count and only the 2 test samples staged squeezed its panels into the left two
  fifths of the slide and stranded `supxlabel` out to the right of them. Any future grid
  wants `ncols = min(cap, len(sample_order))`.
- **A colorbar inherits its mappable's alpha.** A semi-transparent scatter draws the
  ramp blended toward white — washed out and nothing like the viridis it is. Points are
  opaque now, which sidesteps it; drop them below 1 and the colorbar needs
  `cbar.solids.set_alpha(1.0)`.
- **Mito colour steps across 30–70 and holds both ends flat.** The distribution is
  bottom-heavy (median ~4%), so a full 0–100 scale burns its range separating cells
  that are all equally fine and leaves the manifold a near-uniform purple.
  `MITO_BINS = np.arange(30, 71, 10)` puts the resolution in the band where the
  threshold call actually lives. Below 30 nothing is being decided; above 70 a cell is
  dead either way and 75% vs 90% changes nothing.
  - **The top end of 70 is generous, and known to be.** It came from the retention
    curves, which put the gap between the live and dead populations out at ~72; the full
    cohort then showed the 60–70 cells sitting in the debris cloud, so they are dead and
    70 is late. It was left at 70 deliberately — the bands are a reading aid, nothing
    downstream reads them, the report filters nothing, and the threshold question is
    settled (**50, fixed, all samples**). Moving it changes no decision, and by the same
    evidence 50–60 is probably debris too, so there is no principled place to stop.
  - **Why the curve was wrong and the scatter right: marginal vs joint.** The retention
    curve only saw mito, so a thin tail at 60–70 read as "still in the gap". The scatter
    sees those same cells at 100–300 genes — debris regardless of where the mito density
    thins. When the two disagree, the joint view wins.
  - **Fixed rather than per-sample, because the confound was checked and did not hold.**
    At n=1 per condition `obese_id_23` carried ~2× `normal_id_20`'s high-mito fraction,
    which would have made a fixed cutoff discard systematically more from one arm and
    handed downstream DE a technical covariate. Across all 15 that was sample variance.
    No condition-correlated discard rate, so no case for an adaptive per-sample
    threshold (miQC, median + 3×MAD) — **do not re-raise this without new evidence.**
  - **The steps are the point, not decoration.** This panel exists to pick a mito
    cutoff, and a continuous ramp makes you eyeball where one colour *becomes*
    another. Four 10-wide steps on round numbers put the candidate cutoffs in the
    legend, so a band you can see is a number you can name. Do not smooth it back out.
  - **The bin edges and the annotated cutoffs are the same list** (`MITO_CUTOFFS =
    list(MITO_BINS)`), and that correspondence is load-bearing: every colorbar tick is a
    row in `mito_tally`'s annotation, so the colour shows where those cells sit and the
    number says how many. Changing one without the other breaks the pairing that makes
    the slide readable.
  - **The tally is cumulative, not per-band** — a cut at 40 discards everything above
    it, not just 40–50 — so a `≥N` row covers every band from that tick up. It sits
    bottom-right, the one reliably empty corner, since the manifold runs bottom-left to
    top-right and the debris sits low and left of it.
  - **Collapsing the steps to a single band was tried and rejected.** Only the red end
    is spatially separated — the mid cells hug the underside of the manifold at every
    step width tried (5, 10) — which is an argument that the *data* has little structure
    there, not that the bands should go. Merging them breaks the tick-to-row pairing
    above and buys nothing.
  - **This is the inverse of the trap, not a repeat of it.** Clipping the *top* at 25
    once collapsed a 25% cell into a 90% one and hid the dead population entirely.
    Clipping at 70 collapses only cells already past saving. Direction is everything
    here — do not "fix" this back to a full-range scale.
  - **The steps are Okabe-Ito, not a sampled colormap** (`MITO_STEP_COLORS`). A
    sequential ramp is built so neighbouring steps blend, which is precisely wrong when
    the question is which band a point is in. The cost — losing the "higher = worse"
    read, so you consult the legend — is worth it at this few bands. Do not "restore" a
    sequential colormap here.
  - **Both flat ends must stay distinct from the step colours**, or the boundary is
    invisible and a 29% cell looks like a 31% one. Grey under makes the manifold recede
    (it is not the question); red over reads as dead against steps that never go red.
    That last constraint is why the orange-to-crimson end of Okabe-Ito is unusable for
    the steps — any new step colour has to clear both grey and red.
  - **`BoundaryNorm(..., extend="both")` maps the out-of-range regions to the
    colormap's own first and last entries**, so grey and red are ordinary members of
    the `ListedColormap` rather than `set_under`/`set_over`. Hence the list is
    `len(MITO_BINS) + 1` — four bins plus two ends — and the colorbar must *not* be
    passed `extend`, since it already takes it from the norm.
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

### The cluster step applies the cohort's one threshold

`cluster` is where the QC report's reading finally becomes a filter. It is per-sample
(one h5ad in, one out) and reads the same `create_adata` handoff sheet the report does.

- **`MITO_MAX = 50` is a module constant, not a CLI flag.** The threshold question is
  settled — 50, fixed, all samples — and the reasoning is under the mito colour bullet
  above. It is deliberately not exposed as a param because a flag invites the
  per-sample variation that was checked across all 15 samples and rejected. The
  direction matches the report's `dropped ≥N%` column: the report drops `>= N`, so
  `cluster` keeps `< N`. Keep those in step.
- **No `min_genes` / `min_cells` floor, deliberately.** The low-gene debris the QC
  scatter shows (100–300 genes) survives the mito cut. It is left in so it forms its
  own clusters and can be dropped by cluster identity after looking — a more informed
  call than a blanket per-cell floor, and reversible in a way a filter is not.
- **Use `doublet_score`; do not trust `predicted_doublet`.** Scrublet runs flag-only
  (it filters nothing, same spirit as the report), but it picks the threshold behind
  `predicted_doublet` per sample, and on this data that lands somewhere different every
  time: **0.68 for `normal_id_20` vs 0.30 for `obese_id_23`**, on samples whose
  *simulated* score distributions are near-identical (90th pct 0.405 vs 0.465). The
  0.68 sits above that sample's 99.9th observed percentile, so it flagged exactly 1
  cell out of 15,389. The score is comparable across samples; the thresholded call is
  not — it is the same per-sample-adaptive-threshold trap the mito cut avoids. Left
  unpinned rather than fixed cohort-wide because choosing that number wants the
  look-at-all-15 evidence the mito call got. **Revisit once the cohort has run.**
- **HVGs are selected but not subset out**, and there is no `sc.pp.scale`. Scaling
  densifies X across all ~38.6k genes (nothing is filtered upstream, so the var axis is
  full-width), and scanpy's current guidance is to run PCA on log-normalized data
  directly. PCA takes the HVG mask, so carrying every gene costs little and keeps DE's
  inputs intact. `flavor="seurat"` (the default) wants the log1p data it is handed here
  — `seurat_v3` would want raw counts plus `scikit-misc`, which the container lacks.
- **The Leiden sweep writes one obs column per resolution** (`leiden_0.2` … `leiden_1.2`).
  The right resolution is a judgement call made by looking, and a sweep makes that one
  render instead of one pipeline run per value. `--resolutions '0.4 0.8'` overrides it;
  quote it, it is one param with many values.
- **Input and output are both `<sample>.h5ad`**, and stage-in is a symlink into
  `create_adata`'s publish dir — so writing the output at the top level would clobber
  the upstream artifact *through the link*. `stageAs: 'input/*'` keeps the two names
  from ever meeting. Do not "simplify" that away.

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
planned until per-sample clustering says whether merging is warranted. Both
`create_adata` and `cluster` are strictly one-sample-in / one-h5ad-out; the per-sample
clustered h5ads are the current endpoint.

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

The report lands at `<outdir>/qc_report/qc_report.pptx` (~193 KB with 2 samples), so it
is a single `scp` to view off the cluster.

`cluster` reads the **same** `create_adata_samplesheet.csv` the report does — it consumes
the raw h5ads, not the report. The two are independent consumers of `create_adata` and
can run in either order, but read the report first: it is where the mito cutoff `cluster`
applies was chosen.

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
