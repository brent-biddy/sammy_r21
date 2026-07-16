process CLUSTER_REPORT {
    tag "CLUSTER_REPORT"

    // No sample in the path — this is one fan-in task over the whole cohort.
    // 'copy' not 'link': the report is small and is the thing you scp off the
    // cluster, so it should survive the work dir being cleaned.
    publishDir "${params.outdir}/cluster_report", mode: 'copy'

    input:
    // Every sample's clustered h5ad, staged flat into the work dir. The notebook globs
    // *.h5ad from its own directory, so staging IS the input contract — there are no
    // params to pass and nothing to keep in sync with the notebook.
    //
    // Only cluster's h5ads are staged, so the flat glob is safe. Staging create_adata's
    // alongside them would collide: both name their output <sample>.h5ad, and Nextflow
    // would silently rename one rather than fail. A future notebook wanting both wants
    // stageAs subdirs, as cluster.nf does.
    path h5ads
    path notebook
    // Quarto resolves reference-doc relative to the qmd's own directory, so the
    // template has to land beside the staged notebook, not at its repo path.
    path template

    output:
    path "cluster_report.pptx", emit: report

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
    touch cluster_report.pptx
    """
}
