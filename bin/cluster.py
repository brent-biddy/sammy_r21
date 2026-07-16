#!/usr/bin/env python3
"""
cluster.py - Cluster one sample's h5ad and write the clustered store.

Consumes a create_adata h5ad (raw counts + QC metrics, nothing filtered) and runs
the standard scanpy pass: mito filter, doublet scoring, normalization, HVG, PCA,
neighbours, UMAP, and a Leiden sweep. Samples are clustered individually — there
is no concat step — so this is strictly one-sample-in / one-h5ad-out.

Writes <sample>.h5ad into the current working directory, alongside timing and
session info files.

Usage:
    cluster.py --sample normal_id_1 --path normal_id_1.h5ad
    cluster.py --sample normal_id_1 --path normal_id_1.h5ad --resolutions 0.4 0.8
"""

import argparse

import scanpy as sc
import session_info

from timer import timer, timing_summary

# Cells with >= this pct_counts_mt are discarded. Fixed cohort-wide rather than
# adaptive per-sample: the per-sample high-mito fraction was checked across all 15
# samples and is sample variance, not condition-correlated, so a fixed cut does not
# hand downstream DE a technical covariate. Deliberately not a CLI flag — exposing
# it would invite the per-sample variation that was considered and rejected.
# Mirrors the qc_report's `dropped >=N%` column, which is the number this cut was
# chosen against — the report drops >= N, so we keep < N. Keep those in step.
MITO_MAX = 50

DEFAULT_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

# Genes used for PCA. Selected but NOT subset out: DE later wants every gene, and
# with no gene filter upstream the full var axis costs little to carry.
N_TOP_GENES = 2000

# Counts per cell after normalization, before log1p — the conventional 10k.
TARGET_SUM = 1e4


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster one sample's h5ad and write the clustered store"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier, of the form <condition>_id_<study_id> (e.g. normal_id_1)",
    )
    parser.add_argument(
        "--path",
        required=True,
        help="Input h5ad from create_adata (raw counts, unfiltered)",
    )
    parser.add_argument(
        "--resolutions",
        nargs="+",
        type=float,
        default=DEFAULT_RESOLUTIONS,
        help="Leiden resolutions to sweep; each lands in its own obs column "
        "(default: %(default)s)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    output_path = f"{args.sample}.h5ad"

    print(f"Sample:      {args.sample}")
    print(f"Input:       {args.path}")
    print(f"Output:      {output_path}")
    print(f"Resolutions: {args.resolutions}")

    with timer("Read h5ad"):
        adata = sc.read_h5ad(args.path)
    print(f"Loaded:      {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # pct_counts_mt is already in obs — create_adata annotated it and filtered
    # nothing, so this is the first and only place the cohort's cutoff is applied.
    with timer("Mito filter"):
        n_before = adata.n_obs
        adata = adata[adata.obs["pct_counts_mt"] < MITO_MAX].copy()
        n_dropped = n_before - adata.n_obs
    print(
        f"Mito >={MITO_MAX}%:  dropped {n_dropped:,} of {n_before:,} cells "
        f"({n_dropped / n_before * 100:.1f}%), {adata.n_obs:,} remain"
    )

    # No min_genes / min_cells floor, deliberately. The QC scatter shows a low-gene
    # debris cloud (100-300 genes) that survives the mito cut; it is left in so it
    # forms its own clusters and can be dropped by cluster identity after looking,
    # which is a more informed call than a blanket per-cell floor.

    # Doublet scoring on raw counts, before any normalization — scrublet simulates
    # doublets by summing count vectors and needs the counts to do it. Flag only:
    # writes doublet_score and predicted_doublet to obs and filters nothing, same
    # read-only spirit as the QC report.
    #
    # Use doublet_score, NOT predicted_doublet. Scrublet picks the threshold behind
    # predicted_doublet per sample by looking for a bimodal split in the simulated
    # scores, and on this data it lands somewhere different every time: 0.68 for
    # normal_id_20 vs 0.30 for obese_id_23, on samples whose simulated distributions
    # are near-identical (90th pct 0.405 vs 0.465). The 0.68 sits above that sample's
    # 99.9th observed percentile, so it flagged exactly 1 cell. The score is
    # comparable across samples; the thresholded call is not. Left in rather than
    # pinned to a fixed threshold because choosing that number wants the same
    # look-at-all-15 evidence the mito cut got — revisit once the cohort has run.
    with timer("Scrublet"):
        sc.pp.scrublet(adata)
    n_doublet = int(adata.obs["predicted_doublet"].sum())
    print(
        f"Scrublet:    flagged {n_doublet:,} predicted doublets "
        f"({n_doublet / adata.n_obs * 100:.1f}%) — not filtered; "
        f"prefer doublet_score, this call is auto-thresholded per sample"
    )

    # Stash the counts before normalize_total overwrites X in place. DE and any
    # count-based method downstream need them back.
    adata.layers["counts"] = adata.X.copy()

    with timer("Normalize"):
        sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
        sc.pp.log1p(adata)

    # flavor="seurat" (the default) expects the log1p data it is being handed here.
    # seurat_v3 would want raw counts and the scikit-misc package, which the
    # container does not carry.
    with timer("HVG"):
        sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
    print(f"HVG:         {int(adata.var['highly_variable'].sum()):,} genes selected")

    # No sc.pp.scale: it densifies X (all ~36k genes, since nothing was filtered out
    # upstream), and scanpy's current guidance is to run PCA on log-normalized data
    # directly. PCA uses the HVG mask, so the unscaled full var axis costs nothing.
    with timer("PCA"):
        sc.tl.pca(adata, n_comps=50)

    with timer("Neighbors"):
        sc.pp.neighbors(adata)

    with timer("UMAP"):
        sc.tl.umap(adata)

    # One column per resolution rather than one run at a chosen resolution: the
    # right resolution is a judgement call made by looking, and a sweep makes that
    # one render instead of one pipeline run per value.
    with timer("Leiden sweep"):
        for resolution in args.resolutions:
            key = f"leiden_{resolution}"
            # flavor="igraph" with n_iterations=2 is scanpy's recommended pairing;
            # the legacy leidenalg default is deprecated and warns.
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            print(f"  {key}: {adata.obs[key].nunique()} clusters")

    with timer("Write h5ad"):
        adata.write_h5ad(output_path)
    print(f"Written to {output_path}")

    timing_summary(path=f"{args.sample}_timing.tsv")

    session_info_path = f"{args.sample}_session_info.txt"
    session_info.show(write_req_file=True, req_file_name=session_info_path)
    print(f"Session info written to {session_info_path}")


if __name__ == "__main__":
    main()
