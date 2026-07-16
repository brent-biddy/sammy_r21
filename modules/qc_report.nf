process QC_REPORT {
    tag "QC_REPORT"

    // No sample in the path — this is one fan-in task over the whole cohort.
    // 'copy' not 'link': the report is small and is the thing you scp off the
    // cluster, so it should survive the work dir being cleaned.
    publishDir "${params.outdir}/qc_report", mode: 'copy'

    input:
    // Every sample's h5ad, staged flat into the work dir. The notebook globs *.h5ad
    // from its own directory, so staging IS the input contract — there are no params
    // to pass and nothing to keep in sync with the notebook's parameters cell.
    path h5ads
    path notebook
    // Quarto resolves reference-doc relative to the qmd's own directory, so the
    // template has to land beside the staged notebook, not at its repo path.
    path template

    output:
    path "qc_report.pptx", emit: report

    script:
    // Redirect caches and temp files into the task dir: an OSCER compute node's /tmp is
    // read-only, so anything defaulting there fails. Every process needs this preamble.
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --output-dir .
    """

    stub:
    """
    touch qc_report.pptx
    """
}
