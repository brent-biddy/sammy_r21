// Published output directory for this step's per-sample artifacts. Single-sourced
// here so the publishDir directive, the emitted h5ad, and the handoff samplesheet
// row all reference the same location and cannot drift apart.
def clusterPublishDir(sample) {
    "${params.outdir}/${sample}/cluster"
}

process CLUSTER {
    tag "${sample}"

    // saveAs drops the per-sample row fragment from the published dir; it is only
    // needed on the channel for main.nf to collectFile into the aggregate sheet.
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the h5ad.
    publishDir { clusterPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    // stageAs an input/ subdir: the input h5ad and the output h5ad are both named
    // <sample>.h5ad, and stage-in is a symlink into create_adata's publish dir — so
    // writing the output at the top level would clobber the upstream artifact
    // through the link. The subdir keeps the two names from ever meeting.
    tuple val(sample), path(h5ad, stageAs: 'input/*')

    output:
    tuple val(sample), path("${sample}.h5ad"), emit: artifacts
    path "${sample}_timing.tsv", emit: timing
    path "${sample}_session_info.txt", emit: session_info
    // One `sample,path` line pointing at the published h5ad; main.nf collectFiles
    // these into a ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    // Optional args built as a list + conditional append, so an unset param simply
    // leaves the flag off and cluster.py's own default applies.
    def args = []
    if (params.resolutions) args << "--resolutions ${params.resolutions}"
    def args_str = args.join(' ')
    // Redirect caches and temp files into the task dir: an OSCER compute node's /tmp is
    // read-only, so anything defaulting there fails. Every process needs this preamble.
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cluster.py --sample ${sample} --path ${h5ad} ${args_str}

    printf '%s' '${sample},${clusterPublishDir(sample)}/${sample}.h5ad' > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    touch ${sample}.h5ad
    touch ${sample}_timing.tsv
    touch ${sample}_session_info.txt

    printf '%s' '${sample},${clusterPublishDir(sample)}/${sample}.h5ad' > ${sample}.samplesheet_row.csv
    """
}
