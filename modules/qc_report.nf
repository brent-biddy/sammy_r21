// Published output directory for the cohort QC report. Unlike the per-sample steps
// there is no sample in the path — this is one fan-in task over the whole cohort.
def qcReportPublishDir() {
    "${params.outdir}/qc_report"
}

process QC_REPORT {
    tag "QC_REPORT"

    // 'copy' not 'link': the report is small and is the thing you scp off the
    // cluster, so it should survive the work dir being cleaned.
    publishDir { qcReportPublishDir() }, mode: 'copy'

    input:
    // Every sample's h5ad, staged flat into the work dir. The notebook globs *.h5ad
    // from its own directory, so staging IS the input contract — there are no params
    // to pass and nothing to keep in sync with the notebook's parameters cell.
    path h5ads
    path notebook

    output:
    path "qc_report.html", emit: report

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --output-dir .
    """

    stub:
    """
    touch qc_report.html
    """
}
