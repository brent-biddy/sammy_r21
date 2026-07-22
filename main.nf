// sammy_r21 — scRNA-seq analysis pipeline (normal vs obese).
//
// All steps run through this single entry point, selected with --step.
//
// Steps:
//   create_adata              samplesheet: sample, path
//   qc_report                 samplesheet: sample, path  (h5ads from create_adata)
//   cluster                   samplesheet: sample, path  (h5ads from create_adata)
//   cluster_report            samplesheet: sample, path  (h5ads from cluster)
//   sample_summary            samplesheet: sample, path  (h5ads from cluster)

include { CREATE_ADATA }   from './modules/create_adata'
include { QC_REPORT }      from './modules/qc_report'
include { CLUSTER }        from './modules/cluster'
include { CLUSTER_REPORT } from './modules/cluster_report'
include { SAMPLE_SUMMARY } from './modules/sample_summary'

// ── Entry workflow ────────────────────────────────────────────────────────────

workflow {
    // Inline, not a script-level `def`: Nextflow's strict syntax rejects statements
    // mixed with script declarations, and `nextflow config .` does not compile main.nf
    // so it will not catch that — only a run or a -stub run will.
    def valid_steps = 'create_adata, qc_report, cluster, cluster_report, sample_summary'

    if (!params.step) error "Please provide --step <name>. Valid steps: ${valid_steps}"

    if      (params.step == 'create_adata')   create_adata()
    else if (params.step == 'qc_report')      qc_report()
    else if (params.step == 'cluster')        cluster()
    else if (params.step == 'cluster_report') cluster_report()
    else if (params.step == 'sample_summary') sample_summary()
    else error "Unknown --step '${params.step}'. Valid steps: ${valid_steps}"
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

    QC_REPORT(
        qcReportH5ads,
        file("${projectDir}/notebooks/qc_report.qmd"),
        file("${projectDir}/resources/ouhsc_ppt_template.pptx"),
    )
}

// ── cluster ───────────────────────────────────────────────────────────────────

workflow cluster {
    if (!params.samplesheet) error "Please provide --samplesheet"

    // Point --samplesheet at create_adata's published handoff sheet. Samples are
    // clustered individually — one h5ad in, one h5ad out, no concat — so this stays
    // a per-sample fan-out, unlike qc_report.
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | CLUSTER

    // Aggregate the per-sample rows into a handoff samplesheet, same contract as
    // create_adata's — the published path lives in the module.
    CLUSTER.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'cluster_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── cluster_report ────────────────────────────────────────────────────────────

workflow cluster_report {
    if (!params.samplesheet) error "Please provide --samplesheet"

    // Point --samplesheet at cluster's published handoff sheet — this report reads the
    // clustered h5ads (X_umap and the leiden columns), not create_adata's raw ones.
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
        .set { clusterReportH5ads }

    CLUSTER_REPORT(
        clusterReportH5ads,
        file("${projectDir}/notebooks/cluster_report.qmd"),
        file("${projectDir}/resources/ouhsc_ppt_template.pptx"),
    )
}

// ── sample_summary ────────────────────────────────────────────────────────────

workflow sample_summary {
    if (!params.samplesheet) error "Please provide --samplesheet"

    // Point --samplesheet at cluster's published handoff sheet — this reads the same
    // clustered h5ads as cluster_report, but tabulates them (cell accounting, per-cell
    // QC, clusters at each sample's chosen resolution) instead of plotting them.
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.path) error "Samplesheet row missing 'path': ${row}"
            file(row.path)
        }                            // path(h5ad)
        // sort so the staged order — and therefore the table row order fallback — is
        // reproducible regardless of the order tasks happen to finish upstream.
        .toSortedList()              // one list of every h5ad: fans in to a single task
        .set { sampleSummaryH5ads }

    SAMPLE_SUMMARY(
        sampleSummaryH5ads,
        file("${projectDir}/notebooks/sample_summary.qmd"),
        file("${projectDir}/resources/ouhsc_ppt_template.pptx"),
        file("${projectDir}/assets/chosen_resolutions.csv"),
    )
}
