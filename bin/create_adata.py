#!/usr/bin/env python3
"""
create_adata.py - Convert a Cell Ranger count matrix to an AnnData h5ad store.

Reads a Cell Ranger filtered_feature_bc_matrix directory (matrix.mtx.gz,
barcodes.tsv.gz, features.tsv.gz) and writes an h5ad store carrying the raw
counts plus per-cell and per-gene QC metrics. Nothing is filtered — the h5ad is
the raw artifact consumed by downstream clustering, which owns its own QC
thresholds. Mitochondrial genes are detected from the gene symbols themselves,
so no species argument is needed.

The sample identifier and its experimental condition are attached to every cell
in obs, so each h5ad is self-describing and downstream code never has to recover
the design by parsing sample IDs or re-reading the samplesheet.

Writes <sample>.h5ad into the current working directory, alongside timing and
session info files.

Usage:
    create_adata.py --sample S1 --condition normal --path /data/S1/outs/filtered_feature_bc_matrix
"""

import argparse

import pandas as pd
import scanpy as sc
import session_info

from timer import timer, timing_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a Cell Ranger count matrix to an AnnData h5ad store"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument(
        "--condition",
        required=True,
        help="Experimental condition for this sample (e.g. normal, obese)",
    )
    parser.add_argument(
        "--path",
        required=True,
        help="Cell Ranger filtered_feature_bc_matrix directory",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    output_path = f"{args.sample}.h5ad"

    print(f"Sample:    {args.sample}")
    print(f"Condition: {args.condition}")
    print(f"Input:     {args.path}")
    print(f"Output:    {output_path}")

    # var_names="gene_symbols" makes downstream gene lookups readable; the Ensembl
    # IDs stay available in var["gene_ids"]. cache=False keeps scanpy from writing
    # an h5ad cache into the work dir, which would be a second copy of the matrix.
    with timer("Read matrix"):
        adata = sc.read_10x_mtx(args.path, var_names="gene_symbols", cache=False)

    # Gene symbols are not unique in a Cell Ranger reference (distinct Ensembl IDs
    # can share a symbol). Suffix the duplicates so var_names can index.
    adata.var_names_make_unique()

    # Categorical (not object) so scanpy's groupby/plotting treats these as discrete
    # and so a later concat unions the categories rather than falling back to object.
    # Each h5ad carries exactly one sample and one condition, hence one category each.
    adata.obs["sample"] = pd.Categorical([args.sample] * adata.n_obs)
    adata.obs["condition"] = pd.Categorical([args.condition] * adata.n_obs)

    print(f"Loaded:    {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    with timer("QC metrics"):
        # Case-insensitive so human (MT-ND1) and mouse (mt-Nd1) both match with no
        # species flag to pass. The trailing hyphen is what keeps this specific:
        # it excludes the metallothioneins (MT1A, Mt2), which carry no hyphen.
        adata.var["mt"] = adata.var_names.str.match(r"^[Mm][Tt]-")
        # qc_vars=["mt"] adds pct_counts_mt to obs. percent_top reports the share of
        # counts in the top-N genes per cell, a library-complexity signal.
        # log1p=False — the raw totals are what downstream thresholds are set on.
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=(10, 20, 50, 150),
            log1p=False,
            inplace=True,
        )

    n_mito = int(adata.var["mt"].sum())
    print(f"Flagged {n_mito:,} mitochondrial genes.")
    # A real reference has ~13 protein-coding mito genes. Zero means the symbols are
    # not what we expect (Ensembl IDs as var_names, a custom reference), and every
    # pct_counts_mt would silently read 0 — worth flagging rather than filtering on.
    if n_mito == 0:
        print("WARNING: no MT-/mt- genes found; pct_counts_mt will be 0 for all cells.")
    print(f"Median counts/cell: {adata.obs['total_counts'].median():,.0f}")
    print(f"Median genes/cell:  {adata.obs['n_genes_by_counts'].median():,.0f}")
    print(f"Median pct mito:    {adata.obs['pct_counts_mt'].median():.2f}%")

    with timer("Write h5ad"):
        adata.write_h5ad(output_path)
    print(f"Written to {output_path}")

    timing_summary(path=f"{args.sample}_timing.tsv")

    session_info_path = f"{args.sample}_session_info.txt"
    session_info.show(write_req_file=True, req_file_name=session_info_path)
    print(f"Session info written to {session_info_path}")


if __name__ == "__main__":
    main()
