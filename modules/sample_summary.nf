process SAMPLE_SUMMARY {
    tag "SAMPLE_SUMMARY"

    // No sample in the path — this is one fan-in task over the whole cohort.
    // 'copy' not 'link': the deck and its CSV are small and are what you scp off the
    // cluster, so they should survive the work dir being cleaned.
    publishDir "${params.outdir}/sample_summary", mode: 'copy'

    input:
    // Every sample's clustered h5ad, staged flat into the work dir. The notebook globs
    // *.h5ad from its own directory, so staging IS the input contract — same as
    // cluster_report. Only cluster's h5ads are staged, so the flat glob is safe.
    path h5ads
    path notebook
    // Quarto resolves reference-doc relative to the qmd's own directory, so the
    // template has to land beside the staged notebook, not at its repo path.
    path template
    // The curated per-sample resolution table, staged flat beside the notebook. The
    // notebook reads it by fixed name (chosen_resolutions.csv); n_clusters is reported
    // at the resolution chosen for each sample, so this is a required input, not a param.
    path resolutions

    output:
    path "sample_summary.pptx", emit: report
    // The full flat table as a side output — the deck splits it across slides for
    // legibility; the CSV keeps every column in one place for downstream use.
    path "sample_summary.csv", emit: table

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
    touch sample_summary.pptx sample_summary.csv
    """
}
