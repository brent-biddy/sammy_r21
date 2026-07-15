#!/usr/bin/env python3

import anndata as ad
import scanpy as sc
from pathlib import Path


def create_anndata(input_path: Path) -> ad.AnnData:
    """Load the raw count matrix and return an AnnData object.

    Args:
        input_path: Path to the CellRanger filtered_feature_bc_matrix directory.

    Returns:
        AnnData with obs = cells, var = genes, and X = raw counts.
    """
    # TODO: confirm the input layout — this assumes CellRanger mtx output
    adata = sc.read_10x_mtx(input_path, var_names="gene_symbols")
    adata.var_names_make_unique()
    return adata


def add_qc_vars(adata: ad.AnnData) -> None:
    """Flag mitochondrial, ribosomal, and hemoglobin genes and compute QC metrics.

    Adds boolean columns to adata.var (mt, ribo, hb) and per-cell QC metrics
    to adata.obs via sc.pp.calculate_qc_metrics.

    Args:
        adata: AnnData object to annotate in place.
    """
    adata.var["mt"] = adata.var_names.str.match(r"^[Mm][Tt]-")
    adata.var["ribo"] = adata.var_names.str.match(r"^[Rr][Pp][SsLl]")
    adata.var["hb"] = adata.var_names.str.contains(r"^[Hh][Bb][AaBb]-", na=False)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True
    )


def main():
    """Build and write the processed AnnData object.

    Reads the raw counts, attaches per-cell sample metadata, computes QC
    metrics, and writes the result to data/processed/.
    """
    # paths are relative to the repo root — run from there
    input_path = Path("data/raw")  # TODO: point at the actual sample directory
    output_path = Path("data/processed/adata.h5ad")

    adata = create_anndata(input_path)

    # TODO: attach sample metadata (condition, timepoint, replicate, ...) as
    # Categorical obs columns, and set {key}_colors palettes in adata.uns

    add_qc_vars(adata)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {adata.n_obs} cells × {adata.n_vars} genes to {output_path}")
    adata.write_h5ad(output_path)


if __name__ == "__main__":
    main()
