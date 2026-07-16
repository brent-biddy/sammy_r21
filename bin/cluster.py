#!/usr/bin/env python3
"""
cluster.py - Cluster one sample's h5ad and write the clustered store.

Consumes a create_adata h5ad (raw counts + QC metrics, nothing filtered) and runs
the standard scanpy pass: mito filter, doublet scoring, normalization, HVG, scale,
PCA, neighbours, UMAP, and a Leiden sweep. Samples are clustered individually —
there is no concat step — so this is strictly one-sample-in / one-h5ad-out.

The embedding parameters follow scanpy's legacy 2017 clustering tutorial, which
exists to reproduce Seurat's results — see the Embedding parameters block below for
what is copied and what is deliberately not.

Writes <sample>.h5ad into the current working directory, alongside timing and
session info files. The written object keeps its full var axis and its counts
layer; only the embedding is computed on the HVG subset.

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

# ── Embedding parameters ──────────────────────────────────────────────────────
# These mirror scanpy's legacy 2017 clustering tutorial, which exists specifically to
# reproduce Seurat's PBMC3k results:
#   https://scanpy.scverse.org/en/stable/tutorials/basics/clustering-2017.html
# Its QC thresholds are deliberately NOT copied — they are tuned for 2.7k PBMCs
# (mito < 5, genes < 2500) and this cohort's mito cut is settled at 50 cohort-wide.
# Its regress_out(["total_counts", "pct_counts_mt"]) is also deliberately skipped.

N_TOP_GENES = 2000
TARGET_SUM = 1e4        # counts per cell after normalization, before log1p
SCALE_MAX = 10          # clip z-scores here, as the tutorial does
N_COMPS = 50            # PCs computed
N_PCS = 40              # PCs the neighbour graph actually uses
N_NEIGHBORS = 10
# Fixed so a rerun reproduces the same embedding and the same clusters.
RANDOM_STATE = 0


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
    # count-based method downstream need them back, and seurat_v3 HVG reads them here.
    adata.layers["counts"] = adata.X.copy()

    with timer("Normalize"):
        sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
        sc.pp.log1p(adata)

    # seurat_v3 selects on raw counts (hence layer="counts"), unlike the "seurat"
    # flavor which expects log data. It is what the tutorial uses.
    with timer("HVG"):
        sc.pp.highly_variable_genes(
            adata, layer="counts", n_top_genes=N_TOP_GENES, flavor="seurat_v3"
        )
    print(f"HVG:         {int(adata.var['highly_variable'].sum()):,} genes selected")

    # Embed on a HVG-only copy, then carry the results back onto the full object.
    #
    # The tutorial scales every gene into a layer and lets PCA mask to the HVGs. That
    # is numerically identical to this — scale z-scores each gene independently, so
    # subsetting before or after does not move a single PCA input — but it would
    # densify all ~38.6k genes (nothing is filtered upstream, so the var axis is
    # full-width) to compute 2,000 columns' worth of embedding. Scaling only what PCA
    # reads is also what Seurat's ScaleData does by default.
    #
    # Keeping the *published* h5ad full-width matters: cluster_report reads obs and
    # obsm, and any later pseudobulk DE needs layers["counts"] across every gene.
    with timer("Scale + PCA"):
        hvg = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(hvg, max_value=SCALE_MAX)
        sc.tl.pca(hvg, n_comps=N_COMPS, svd_solver="arpack", random_state=RANDOM_STATE)

    with timer("Neighbors"):
        sc.pp.neighbors(hvg, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_STATE)

    with timer("UMAP"):
        sc.tl.umap(hvg, random_state=RANDOM_STATE)

    # One column per resolution rather than one run at a chosen resolution: the
    # right resolution is a judgement call made by looking, and a sweep makes that
    # one render instead of one pipeline run per value.
    with timer("Leiden sweep"):
        for resolution in args.resolutions:
            key = f"leiden_{resolution}"
            # flavor="igraph" with n_iterations=2 is scanpy's recommended pairing and
            # what the tutorial uses; the legacy leidenalg default is deprecated.
            sc.tl.leiden(
                hvg,
                resolution=resolution,
                key_added=key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                random_state=RANDOM_STATE,
            )
            adata.obs[key] = hvg.obs[key]
            print(f"  {key}: {adata.obs[key].nunique()} clusters")

    # Carry the embedding back onto the full object. obs_names are untouched by the
    # var-axis subset, so these stay row-aligned.
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.obsm["X_umap"] = hvg.obsm["X_umap"]

    with timer("Write h5ad"):
        adata.write_h5ad(output_path)
    print(f"Written to {output_path}")

    timing_summary(path=f"{args.sample}_timing.tsv")

    session_info_path = f"{args.sample}_session_info.txt"
    session_info.show(write_req_file=True, req_file_name=session_info_path)
    print(f"Session info written to {session_info_path}")


if __name__ == "__main__":
    main()
