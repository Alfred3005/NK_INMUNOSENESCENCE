"""
Preprocessing Module
====================

This module implements preprocessing steps for single-cell RNA-seq data:
- Normalization (library size, log transformation)
- Scaling and variance stabilization
- Highly variable gene (HVG) selection (batch-aware)
- Dimensionality reduction (PCA)
- Marker gene inclusion

Author: Alfred3005
"""

import logging
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from .utils import setup_logging, log_memory_usage, cleanup_memory


def normalize_total(
    adata: ad.AnnData,
    target_sum: float = 1e4,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Normalize counts per cell (library size normalization).

    Normalizes each cell to have `target_sum` total counts, then log-transforms.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with raw counts
    target_sum : float, default 10000
        Target sum for normalization (e.g., 10000 = CPM/100)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        Normalized AnnData

    Examples
    --------
    >>> adata = normalize_total(adata, target_sum=1e4, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    logger.info(f"Normalizing to {target_sum:.0f} counts per cell...")

    # Store raw counts in layer
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
        logger.info("Raw counts stored in adata.layers['counts']")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)

    logger.info("Normalization complete")

    return adata


def log_transform(
    adata: ad.AnnData,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Log-transform normalized counts: log(count + 1).

    Parameters
    ----------
    adata : AnnData
        Input AnnData with normalized counts
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        Log-transformed AnnData

    Examples
    --------
    >>> adata = log_transform(adata, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    logger.info("Applying log transformation: log(count + 1)...")

    # Check if already log-transformed
    if adata.X.max() < 20:  # Heuristic: log-transformed data is usually < 20
        logger.warning(
            "Data appears to already be log-transformed "
            "(max value < 20). Skipping log transformation."
        )
        return adata

    sc.pp.log1p(adata)

    logger.info("Log transformation complete")

    return adata


def scale_data(
    adata: ad.AnnData,
    max_value: float = 10.0,
    zero_center: bool = True,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Scale data to unit variance and optionally zero mean.

    Z-score transformation: (x - mean) / std
    Clips values to prevent extreme outliers.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    max_value : float, default 10.0
        Clip scaled values to [-max_value, max_value]
    zero_center : bool, default True
        Whether to center to zero mean
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        Scaled AnnData

    Examples
    --------
    >>> adata = scale_data(adata, max_value=10, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    logger.info(
        f"Scaling data (zero_center={zero_center}, max_value={max_value})..."
    )

    # Replace infinite values with zero (if any)
    if hasattr(adata.X, 'data'):  # Sparse matrix
        adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:  # Dense matrix
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    # Scale
    sc.pp.scale(adata, max_value=max_value, zero_center=zero_center)

    logger.info("Scaling complete")

    return adata


def select_highly_variable_genes(
    adata: ad.AnnData,
    config: Dict[str, Any],
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Select highly variable genes (HVGs) with batch-aware correction.

    Also forces inclusion of marker genes (NK markers and aging genes).

    Parameters
    ----------
    adata : AnnData
        Input AnnData (log-normalized)
    config : dict
        Configuration from analysis_params.yaml
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with 'highly_variable' in .var

    Examples
    --------
    >>> config = load_config("config/analysis_params.yaml")
    >>> adata = select_highly_variable_genes(adata, config, logger)
    >>> hvgs = adata.var_names[adata.var['highly_variable']]
    >>> print(f"Selected {len(hvgs)} HVGs")
    """
    if logger is None:
        logger = setup_logging()

    hvg_params = config['feature_selection']

    n_top_genes = hvg_params['n_top_genes']
    flavor = hvg_params['flavor']
    batch_key = hvg_params.get('batch_key', None)
    batch_aware = hvg_params.get('batch_aware', True)

    logger.info(
        f"Selecting {n_top_genes} highly variable genes "
        f"(flavor={flavor}, batch_aware={batch_aware})..."
    )

    # Select HVGs
    if batch_aware and batch_key and batch_key in adata.obs.columns:
        logger.info(f"Using batch key: {batch_key}")

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key,
            subset=False  # Keep all genes initially
        )

        # Log HVG statistics per batch
        if 'highly_variable_nbatches' in adata.var.columns:
            n_batches = adata.obs[batch_key].nunique()
            hvg_batch_dist = adata.var['highly_variable_nbatches'].value_counts().sort_index()

            logger.info(f"HVG distribution across {n_batches} batches:")
            for n_batch, count in hvg_batch_dist.items():
                logger.info(f"  HVG in {n_batch} batches: {count} genes")

            # Optionally filter genes that are HVG in multiple batches
            min_batches = hvg_params.get('min_batches', 1)
            if min_batches > 1:
                adata.var['highly_variable'] = (
                    adata.var['highly_variable'] &
                    (adata.var['highly_variable_nbatches'] >= min_batches)
                )
                logger.info(
                    f"Filtering HVGs present in ≥{min_batches} batches: "
                    f"{adata.var['highly_variable'].sum()} genes"
                )

    else:
        logger.info("Batch-aware HVG selection disabled or batch key not found")

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            subset=False
        )

    # Force include marker genes
    if hvg_params.get('force_include_markers', True):
        nk_markers = hvg_params.get('nk_marker_genes', [])
        aging_markers = hvg_params.get('aging_marker_genes', [])
        all_markers = nk_markers + aging_markers

        markers_in_data = [g for g in all_markers if g in adata.var_names]

        if markers_in_data:
            adata.var.loc[markers_in_data, 'highly_variable'] = True
            logger.info(
                f"Forced inclusion of {len(markers_in_data)} marker genes "
                f"({len(nk_markers)} NK + {len(aging_markers)} aging)"
            )

            markers_not_found = set(all_markers) - set(markers_in_data)
            if markers_not_found:
                logger.warning(
                    f"{len(markers_not_found)} marker genes not found in data: "
                    f"{', '.join(list(markers_not_found)[:10])}"
                    f"{'...' if len(markers_not_found) > 10 else ''}"
                )

    # Summary
    n_hvgs = adata.var['highly_variable'].sum()
    logger.info(
        f"Selected {n_hvgs} highly variable genes "
        f"({n_hvgs/adata.n_vars*100:.1f}% of {adata.n_vars} total genes)"
    )

    return adata


def run_pca(
    adata: ad.AnnData,
    n_comps: int = 50,
    use_highly_variable: bool = True,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Perform Principal Component Analysis (PCA).

    Parameters
    ----------
    adata : AnnData
        Input AnnData (scaled, with HVGs selected)
    n_comps : int, default 50
        Number of principal components to compute
    use_highly_variable : bool, default True
        Whether to use only HVGs for PCA
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with PCA in .obsm['X_pca'] and loadings in .varm['PCs']

    Examples
    --------
    >>> adata = run_pca(adata, n_comps=50, logger=logger)
    >>> print(adata.obsm['X_pca'].shape)  # (n_cells, 50)
    """
    if logger is None:
        logger = setup_logging()

    logger.info(f"Computing PCA ({n_comps} components)...")

    if use_highly_variable:
        if 'highly_variable' not in adata.var.columns:
            logger.warning("No HVGs selected. Using all genes for PCA.")
            use_highly_variable = False
        else:
            n_hvgs = adata.var['highly_variable'].sum()
            logger.info(f"Using {n_hvgs} highly variable genes for PCA")

    # Run PCA
    sc.tl.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable,
        svd_solver='arpack',
        random_state=42
    )

    # Log variance explained
    var_ratio = adata.uns['pca']['variance_ratio']
    cumsum = np.cumsum(var_ratio)

    logger.info(
        f"Variance explained by first {n_comps} PCs: {cumsum[-1]*100:.1f}%"
    )
    logger.info(f"  PC1-10: {cumsum[9]*100:.1f}%")
    logger.info(f"  PC1-20: {cumsum[19]*100:.1f}%")
    logger.info(f"  PC1-30: {cumsum[29]*100:.1f}%")

    return adata


def run_umap(
    adata: ad.AnnData,
    config: Dict[str, Any],
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Compute UMAP embedding for visualization.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with PCA computed
    config : dict
        Configuration from analysis_params.yaml
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with UMAP in .obsm['X_umap']

    Examples
    --------
    >>> config = load_config("config/analysis_params.yaml")
    >>> adata = run_umap(adata, config, logger)
    """
    if logger is None:
        logger = setup_logging()

    umap_params = config['dimensionality_reduction']['umap']

    logger.info("Computing UMAP embedding...")

    # Compute neighbors first (required for UMAP)
    sc.pp.neighbors(
        adata,
        n_neighbors=umap_params['n_neighbors'],
        metric=umap_params['metric'],
        random_state=umap_params['random_state']
    )

    # Compute UMAP
    sc.tl.umap(
        adata,
        min_dist=umap_params['min_dist'],
        spread=umap_params['spread'],
        random_state=umap_params['random_state']
    )

    logger.info("UMAP computation complete")

    return adata


def run_leiden_clustering(
    adata: ad.AnnData,
    resolution: float = 0.5,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Perform Leiden clustering for QC and visualization.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with neighbors computed
    resolution : float, default 0.5
        Clustering resolution (higher = more clusters)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with 'leiden' cluster assignments in .obs

    Examples
    --------
    >>> adata = run_leiden_clustering(adata, resolution=0.5, logger=logger)
    >>> print(adata.obs['leiden'].value_counts())
    """
    if logger is None:
        logger = setup_logging()

    logger.info(f"Performing Leiden clustering (resolution={resolution})...")

    # Run Leiden
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=42
    )

    n_clusters = adata.obs['leiden'].nunique()
    logger.info(f"Identified {n_clusters} clusters")

    # Log cluster sizes
    cluster_sizes = adata.obs['leiden'].value_counts().sort_index()
    logger.info("Cluster sizes:")
    for cluster, size in cluster_sizes.items():
        logger.info(f"  Cluster {cluster}: {size:,} cells ({size/adata.n_obs*100:.1f}%)")

    return adata


def preprocess_pipeline(
    adata: ad.AnnData,
    config: Dict[str, Any],
    skip_scaling: bool = False,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Run complete preprocessing pipeline.

    Steps:
    1. Normalize to target sum
    2. Log-transform
    3. Select HVGs (batch-aware)
    4. Scale data (optional)
    5. PCA
    6. UMAP
    7. Leiden clustering

    Parameters
    ----------
    adata : AnnData
        Input AnnData with raw counts (post-QC)
    config : dict
        Configuration from analysis_params.yaml
    skip_scaling : bool, default False
        Whether to skip scaling (useful before scVI)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        Preprocessed AnnData

    Examples
    --------
    >>> config = load_config("config/analysis_params.yaml")
    >>> adata = preprocess_pipeline(adata, config, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    logger.info("="*60)
    logger.info("Starting preprocessing pipeline")
    logger.info("="*60)

    # Get parameters
    norm_params = config['normalization']
    pca_params = config['dimensionality_reduction']['pca']
    cluster_params = config['clustering']

    # Step 1: Normalize
    adata = normalize_total(
        adata,
        target_sum=norm_params['target_sum'],
        logger=logger
    )

    # Step 2: Log-transform
    if norm_params['log_transform']:
        adata = log_transform(adata, logger=logger)

    # Step 3: Select HVGs
    adata = select_highly_variable_genes(adata, config, logger=logger)

    # Step 4: Scale (optional)
    if norm_params['scale'] and not skip_scaling:
        adata = scale_data(
            adata,
            max_value=norm_params['max_value'],
            logger=logger
        )
    elif skip_scaling:
        logger.info("Skipping scaling (will be done by integration method)")

    # Step 5: PCA
    adata = run_pca(
        adata,
        n_comps=pca_params['n_comps'],
        logger=logger
    )

    # Step 6: UMAP
    adata = run_umap(adata, config, logger=logger)

    # Step 7: Leiden clustering
    adata = run_leiden_clustering(
        adata,
        resolution=cluster_params['resolution'],
        logger=logger
    )

    logger.info("="*60)
    logger.info("Preprocessing pipeline complete")
    logger.info("="*60)

    log_memory_usage(logger)

    return adata
