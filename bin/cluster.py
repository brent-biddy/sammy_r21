#!/usr/bin/env python3
"""
cluster.py - Cluster one sample's h5ad and write the clustered store.

Consumes a create_adata h5ad (raw counts + QC metrics, nothing filtered) and runs
the standard scanpy pass: mito filter, doublet scoring, normalization, HVG, PCA,
neighbours, UMAP, and a Leiden sweep. Samples are clustered individually — there
is no concat step — so this is strictly one-sample-in / one-h5ad-out.

The embedding parameters come from scanpy's legacy 2017 clustering tutorial, which
exists to reproduce Seurat's results — see the Embedding parameters block below for
what is taken from it and what is deliberately not.

Writes <sample>.h5ad into the current working directory, alongside timing and
session info files. The written object keeps its full var axis and its counts
layer; PCA masks to the HVGs, so only they reach the embedding.

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
MITO_MAX = 30

# 0.1–1.0 in 0.1 steps, so the sweep is the default and --resolutions rarely needs
# passing. Spelled out as literals rather than [i / 10 for i in range(1, 11)]: the value
# is formatted straight into the obs key (f"leiden_{resolution}"), and 0.3 from division
# is 0.30000000000000004, which would name a column leiden_0.30000000000000004.
DEFAULT_RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# ── Embedding parameters ──────────────────────────────────────────────────────
# These mirror scanpy's legacy 2017 clustering tutorial, which exists specifically to
# reproduce Seurat's PBMC3k results:
#   https://scanpy.scverse.org/en/stable/tutorials/basics/clustering-2017.html
# Its QC thresholds are deliberately NOT copied — they are tuned for 2.7k PBMCs
# (mito < 5, genes < 2500) and this cohort's mito cut is settled at 30 cohort-wide.
# Its regress_out(["total_counts", "pct_counts_mt"]) is also deliberately skipped.
#
# Its sc.pp.scale(max_value=10) is skipped too, following scanpy's current clustering
# tutorial, which goes normalize -> log1p -> HVG -> PCA with no z-scoring:
#   https://scanpy.scverse.org/en/stable/tutorials/basics/clustering.html
# Scaling flattens every gene to unit variance, which hands a rarely-expressed gene
# the same pull on the PCs as a highly variable one; leaving it out keeps PCA weighted
# by the expression variance actually in the data. The counts and log-normalized X are
# unaffected either way — scale only ever touched what PCA read.

N_TOP_GENES = 2000
TARGET_SUM = 1e4        # counts per cell after normalization, before log1p
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


def relabel_by_size(labels):
    """Renumber Leiden labels 1..k by descending cluster size (largest is "1").

    Leiden returns 0-indexed labels whose numbering carries no meaning — cluster "0" is
    not the most populous, and the ids are not comparable across resolutions. Relabel so
    the id encodes rank: cluster 1 is always the biggest, cluster 2 the next, and so on.
    That makes the sweep's columns read consistently and a reader's "cluster 1" mean the
    same thing everywhere. Ties in exact cell count are broken by the original label so
    the mapping is deterministic across reruns. Returns an ordered categorical (1..k) so
    downstream plots and tables sort numerically rather than lexically.
    """
    counts = labels.value_counts()                       # descending by count
    order = sorted(counts.index, key=lambda c: (-counts[c], int(c)))
    mapping = {old: str(rank) for rank, old in enumerate(order, start=1)}
    new_levels = [str(rank) for rank in range(1, len(order) + 1)]
    return (
        labels.map(mapping)
        .astype("category")
        .cat.set_categories(new_levels, ordered=True)
    )


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
    # Record the mito accounting in uns so the clustered h5ad is self-describing:
    # the dropped cells are physically gone from the object, so without this the only
    # way to recover "how many were removed" is to re-read the create_adata h5ad and
    # re-apply the cut. Downstream bookkeeping reads these scalars instead.
    adata.uns["mito_filter"] = {
        "mito_max": MITO_MAX,
        "n_before": n_before,
        "n_dropped": n_dropped,
        "n_after": adata.n_obs,
    }
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

    # Everything from here runs on the full object. pca defaults to
    # mask_var="highly_variable", so it reads only the HVGs without a subset copy and
    # without densifying the other ~36.6k genes. The published h5ad stays full-width,
    # which matters: cluster_report reads obs and obsm, and any later pseudobulk DE
    # needs layers["counts"] across every gene.
    with timer("PCA"):
        sc.tl.pca(adata, n_comps=N_COMPS, svd_solver="arpack", random_state=RANDOM_STATE)

    with timer("Neighbors"):
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_STATE)

    with timer("UMAP"):
        sc.tl.umap(adata, random_state=RANDOM_STATE)

    # One column per resolution rather than one run at a chosen resolution: the
    # right resolution is a judgement call made by looking, and a sweep makes that
    # one render instead of one pipeline run per value.
    with timer("Leiden sweep"):
        for resolution in args.resolutions:
            key = f"leiden_{resolution}"
            # flavor="igraph" with n_iterations=2 is scanpy's recommended pairing and
            # what the tutorial uses; the legacy leidenalg default is deprecated.
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                random_state=RANDOM_STATE,
            )
            # Renumber 1..k by descending size in place — see relabel_by_size. Done here,
            # not in a report, so every consumer of the h5ad (the deck, annotation, any
            # DE) sees the same size-ranked ids.
            adata.obs[key] = relabel_by_size(adata.obs[key])
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
