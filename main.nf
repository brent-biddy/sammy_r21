// sammy_r21 — scRNA-seq analysis pipeline (normal vs obese).
//
// All steps run through this single entry point, selected with --step.
//
// Steps:
//   create_adata              samplesheet: sample, path

include { CREATE_ADATA } from './modules/create_adata'

// ── Entry workflow ────────────────────────────────────────────────────────────

workflow {
    if (!params.step) error "Please provide --step <name>. Valid steps: create_adata"

    if      (params.step == 'create_adata') create_adata()
    else error "Unknown --step '${params.step}'. Valid steps: create_adata"
}

// ── create_adata ──────────────────────────────────────────────────────────────

workflow create_adata {
    if (!params.samplesheet) error "Please provide --samplesheet"

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | CREATE_ADATA

    // Aggregate the per-sample rows the process emits into a ready-to-use handoff
    // samplesheet, so the clustering step can be pointed straight at it instead of
    // hand-building a sample,path CSV. The published path lives in the module (its
    // publishDir and the emitted row share one helper), so main.nf stays agnostic.
    CREATE_ADATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'create_adata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}
