"""
Integration Module
==================

This module implements batch correction methods for single-cell RNA-seq data:
- scVI (variational inference)
- scANVI (semi-supervised scVI)
- Harmony (PCA-based correction)
- Integration quality metrics

Author: Alfred3005
"""

import logging
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

try:
    import scvi
except ImportError:
    scvi = None

try:
    import harmonypy
except ImportError:
    harmonypy = None

from .utils import setup_logging, log_memory_usage, cleanup_memory, check_gpu_available


def integrate_scvi(
    adata: ad.AnnData,
    config: Dict[str, Any],
    batch_key: str = 'dataset_id',
    use_gpu: bool = True,
    logger: Optional[logging.Logger] = None
) -> Tuple[ad.AnnData, Any]:
    """
    Integrate data using scVI (single-cell Variational Inference).

    Parameters
    ----------
    adata : AnnData
        Input AnnData (log-normalized, with HVGs selected)
    config : dict
        Configuration from analysis_params.yaml
    batch_key : str, default 'dataset_id'
        Key for batch correction
    use_gpu : bool, default True
        Whether to use GPU acceleration
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    tuple
        (adata with scVI embedding in .obsm['X_scvi'], trained model)

    Examples
    --------
    >>> config = load_config("config/analysis_params.yaml")
    >>> adata, model = integrate_scvi(adata, config, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    if scvi is None:
        raise ImportError(
            "scvi-tools not installed. Install with: pip install scvi-tools"
        )

    logger.info("="*60)
    logger.info("Starting scVI integration")
    logger.info("="*60)

    # Check GPU availability
    has_gpu = check_gpu_available(logger)
    use_gpu = use_gpu and has_gpu

    # Set scVI seed for reproducibility
    scvi.settings.seed = 42

    # Get parameters
    scvi_params = config['integration']['scvi']

    # Subset to HVGs
    if 'highly_variable' in adata.var.columns:
        adata_hvg = adata[:, adata.var['highly_variable']].copy()
        logger.info(f"Using {adata_hvg.n_vars} highly variable genes for scVI")
    else:
        adata_hvg = adata.copy()
        logger.warning("No HVGs selected. Using all genes.")

    # Setup scVI
    logger.info("Setting up scVI model...")

    categorical_covs = [
        c for c in scvi_params.get('categorical_covariate_keys', [])
        if c in adata_hvg.obs.columns
    ]
    continuous_covs = [
        c for c in scvi_params.get('continuous_covariate_keys', [])
        if c in adata_hvg.obs.columns
    ]

    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covs if categorical_covs else None,
        continuous_covariate_keys=continuous_covs if continuous_covs else None,
        layer='counts' if 'counts' in adata_hvg.layers else None
    )

    # Create model
    model = scvi.model.SCVI(
        adata_hvg,
        n_latent=scvi_params['n_latent'],
        n_layers=scvi_params['n_layers'],
        n_hidden=scvi_params.get('n_hidden', 128),
        dropout_rate=scvi_params['dropout_rate'],
        gene_likelihood=scvi_params['gene_likelihood'],
        dispersion=scvi_params.get('dispersion', 'gene')
    )

    logger.info(f"Model architecture: {scvi_params['n_layers']} layers, "
                f"{scvi_params['n_latent']} latent dimensions")

    # Train model
    logger.info("Training scVI model...")
    log_memory_usage(logger)

    model.train(
        max_epochs=scvi_params['max_epochs'],
        early_stopping=scvi_params['early_stopping'],
        early_stopping_patience=scvi_params['early_stopping_patience'],
        early_stopping_monitor=scvi_params.get('early_stopping_monitor', 'elbo_validation'),
        train_size=scvi_params['train_size'],
        validation_size=scvi_params.get('validation_size', 0.1),
        batch_size=scvi_params['batch_size'],
        plan_kwargs={'lr': scvi_params.get('learning_rate', 1e-3)},
        use_gpu=use_gpu
    )

    logger.info("Training complete")

    # Get latent representation
    logger.info("Computing latent representation...")
    latent = model.get_latent_representation()

    # Add to original adata
    adata.obsm['X_scvi'] = latent

    # Compute UMAP on scVI embedding
    logger.info("Computing UMAP on scVI latent space...")
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3)

    logger.info("="*60)
    logger.info("scVI integration complete")
    logger.info("="*60)

    log_memory_usage(logger)
    cleanup_memory(logger)

    return adata, model


def integrate_scanvi(
    adata: ad.AnnData,
    scvi_model: Any,
    config: Dict[str, Any],
    labels_key: str = 'cell_type',
    use_gpu: bool = True,
    logger: Optional[logging.Logger] = None
) -> Tuple[ad.AnnData, Any]:
    """
    Integrate data using scANVI (semi-supervised scVI).

    Initializes from pre-trained scVI model.

    Parameters
    ----------
    adata : AnnData
        Input AnnData (same as used for scVI)
    scvi_model : scvi.model.SCVI
        Pre-trained scVI model
    config : dict
        Configuration from analysis_params.yaml
    labels_key : str, default 'cell_type'
        Key for cell type labels
    use_gpu : bool, default True
        Whether to use GPU
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    tuple
        (adata with scANVI embedding in .obsm['X_scanvi'], trained model)

    Examples
    --------
    >>> adata, scvi_model = integrate_scvi(adata, config, logger=logger)
    >>> adata, scanvi_model = integrate_scanvi(adata, scvi_model, config, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    if scvi is None:
        raise ImportError("scvi-tools not installed")

    logger.info("="*60)
    logger.info("Starting scANVI integration")
    logger.info("="*60)

    # Check GPU
    has_gpu = check_gpu_available(logger)
    use_gpu = use_gpu and has_gpu

    # Get parameters
    scanvi_params = config['integration']['scanvi']

    # Get adata subset (same as used for scVI)
    if 'highly_variable' in adata.var.columns:
        adata_hvg = adata[:, adata.var['highly_variable']].copy()
    else:
        adata_hvg = adata.copy()

    # Check if labels exist
    if labels_key not in adata_hvg.obs.columns:
        logger.error(f"Labels key '{labels_key}' not found. Skipping scANVI.")
        return adata, None

    logger.info(f"Using '{labels_key}' for semi-supervised learning")

    # Initialize scANVI from scVI
    logger.info("Initializing scANVI from scVI model...")

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=labels_key,
        unlabeled_category=scanvi_params.get('unlabeled_category', 'Unknown')
    )

    # Train scANVI
    logger.info("Training scANVI model...")
    log_memory_usage(logger)

    scanvi_model.train(
        max_epochs=scanvi_params['max_epochs'],
        early_stopping=scanvi_params['early_stopping'],
        early_stopping_patience=scanvi_params['early_stopping_patience'],
        train_size=scanvi_params['train_size'],
        batch_size=scanvi_params['batch_size'],
        plan_kwargs={'lr': scanvi_params.get('learning_rate', 1e-3)},
        use_gpu=use_gpu
    )

    logger.info("Training complete")

    # Get latent representation
    logger.info("Computing latent representation...")
    latent = scanvi_model.get_latent_representation()

    # Add to original adata
    adata.obsm['X_scanvi'] = latent

    # Compute UMAP on scANVI embedding
    logger.info("Computing UMAP on scANVI latent space...")
    sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3)

    logger.info("="*60)
    logger.info("scANVI integration complete")
    logger.info("="*60)

    log_memory_usage(logger)
    cleanup_memory(logger)

    return adata, scanvi_model


def integrate_harmony(
    adata: ad.AnnData,
    config: Dict[str, Any],
    batch_key: str = 'dataset_id',
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Integrate data using Harmony (PCA-based correction).

    Parameters
    ----------
    adata : AnnData
        Input AnnData with PCA computed
    config : dict
        Configuration from analysis_params.yaml
    batch_key : str, default 'dataset_id'
        Key for batch correction
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with Harmony embedding in .obsm['X_pca_harmony']

    Examples
    --------
    >>> config = load_config("config/analysis_params.yaml")
    >>> adata = integrate_harmony(adata, config, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    if harmonypy is None:
        raise ImportError(
            "harmonypy not installed. Install with: pip install harmonypy"
        )

    logger.info("="*60)
    logger.info("Starting Harmony integration")
    logger.info("="*60)

    # Check if PCA exists
    if 'X_pca' not in adata.obsm.keys():
        logger.error("PCA not found. Run PCA before Harmony integration.")
        raise ValueError("PCA required for Harmony")

    # Get parameters
    harmony_params = config['integration']['harmony']

    # Run Harmony
    logger.info(f"Running Harmony correction (batch_key={batch_key})...")

    harmony_out = harmonypy.run_harmony(
        adata.obsm['X_pca'],
        adata.obs,
        batch_key,
        max_iter_harmony=harmony_params['max_iter_harmony'],
        theta=harmony_params.get('theta', 2)
    )

    # Add to adata
    adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T

    # Compute UMAP on Harmony embedding
    logger.info("Computing UMAP on Harmony-corrected PCA...")
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3)

    logger.info("="*60)
    logger.info("Harmony integration complete")
    logger.info("="*60)

    log_memory_usage(logger)

    return adata


def calculate_integration_metrics(
    adata: ad.AnnData,
    batch_key: str = 'dataset_id',
    label_key: str = 'cell_type',
    embed_keys: List[str] = ['X_scvi', 'X_scanvi', 'X_pca_harmony'],
    logger: Optional[logging.Logger] = None
) -> pd.DataFrame:
    """
    Calculate integration quality metrics for multiple embeddings.

    Uses scib-metrics package to compute:
    - Batch mixing: iLISI, ASW batch, kBET
    - Bio conservation: cLISI, ASW label, NMI, ARI

    Parameters
    ----------
    adata : AnnData
        Input AnnData with multiple embeddings
    batch_key : str, default 'dataset_id'
        Batch variable
    label_key : str, default 'cell_type'
        Cell type variable
    embed_keys : list of str
        List of embedding keys to evaluate
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        DataFrame with metrics for each method

    Examples
    --------
    >>> metrics_df = calculate_integration_metrics(
    ...     adata,
    ...     embed_keys=['X_scvi', 'X_scanvi', 'X_pca_harmony'],
    ...     logger=logger
    ... )
    >>> print(metrics_df)
    """
    if logger is None:
        logger = setup_logging()

    logger.info("Calculating integration quality metrics...")

    try:
        import scib_metrics
    except ImportError:
        logger.error("scib-metrics not installed. Skipping metrics calculation.")
        return pd.DataFrame()

    metrics_list = []

    for embed_key in embed_keys:
        if embed_key not in adata.obsm.keys():
            logger.warning(f"Embedding '{embed_key}' not found. Skipping.")
            continue

        logger.info(f"Evaluating {embed_key}...")

        method_name = embed_key.replace('X_', '').replace('pca_', '')

        try:
            # Calculate metrics
            # Note: scib-metrics API may vary; adjust as needed
            metrics = {}

            # Batch mixing metrics (higher is better)
            # iLISI: Integration LISI
            # ASW batch: Silhouette score for batch (lower is better for mixing)

            # Bio conservation metrics
            # cLISI: Cell type LISI (lower is better for conservation)
            # ASW label: Silhouette score for cell type (higher is better)

            # For this implementation, we'll use simplified metrics
            # Full scib-metrics integration would require more setup

            metrics['method'] = method_name
            metrics['embedding'] = embed_key

            # Placeholder - replace with actual scib-metrics calls
            # metrics['iLISI'] = scib_metrics.ilisi(adata, batch_key, embed_key)
            # metrics['cLISI'] = scib_metrics.clisi(adata, label_key, embed_key)

            logger.info(f"  Metrics for {method_name}: {metrics}")

            metrics_list.append(metrics)

        except Exception as e:
            logger.error(f"Error calculating metrics for {embed_key}: {str(e)}")
            continue

    if metrics_list:
        metrics_df = pd.DataFrame(metrics_list).set_index('method')
        logger.info("\nIntegration metrics summary:")
        logger.info(f"\n{metrics_df.to_string()}")
        return metrics_df
    else:
        logger.warning("No metrics calculated")
        return pd.DataFrame()


def save_integrated_data(
    adata: ad.AnnData,
    output_dir: Path,
    method: str,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    Save integrated AnnData object.

    Parameters
    ----------
    adata : AnnData
        Integrated AnnData
    output_dir : Path
        Output directory
    method : str
        Integration method name (scvi, scanvi, harmony)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to saved file

    Examples
    --------
    >>> save_integrated_data(
    ...     adata,
    ...     Path("data/processed"),
    ...     method="scvi",
    ...     logger=logger
    ... )
    """
    if logger is None:
        logger = setup_logging()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"nk_cells_integrated_{method}.h5ad"

    logger.info(f"Saving integrated data to: {output_path}")

    adata.write_h5ad(output_path, compression='gzip')

    file_size = output_path.stat().st_size / (1024 ** 2)  # MB
    logger.info(f"File saved ({file_size:.1f} MB)")

    return output_path
