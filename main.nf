// sammy_r21 — scRNA-seq analysis pipeline (normal vs obese).
//
// All steps run through this single entry point, selected with --step.
//
// Steps:
//   create_adata              samplesheet: sample, path
//   qc_report                 samplesheet: sample, path  (h5ads from create_adata)

include { CREATE_ADATA } from './modules/create_adata'
include { QC_REPORT }    from './modules/qc_report'

// ── Entry workflow ────────────────────────────────────────────────────────────

workflow {
    if (!params.step) error "Please provide --step <name>. Valid steps: create_adata, qc_report"

    if      (params.step == 'create_adata') create_adata()
    else if (params.step == 'qc_report')    qc_report()
    else error "Unknown --step '${params.step}'. Valid steps: create_adata, qc_report"
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

// ── qc_report ─────────────────────────────────────────────────────────────────

workflow qc_report {
    if (!params.samplesheet) error "Please provide --samplesheet"

    // Point --samplesheet at create_adata's published handoff sheet.
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.path) error "Samplesheet row missing 'path': ${row}"
            file(row.path)
        }                            // path(h5ad)
        // sort so the staged order — and therefore the report — is reproducible
        // regardless of the order tasks happen to finish upstream.
        .toSortedList()              // one list of every h5ad: fans in to a single task
        .set { qcReportH5ads }

    QC_REPORT(qcReportH5ads, file("${projectDir}/notebooks/qc_report.qmd"))
}
