"""
Quality Control Module
======================

This module implements rigorous quality control for single-cell RNA-seq data:
- Basic QC metrics calculation (genes, counts, MT%, ribo%, HB%)
- MAD-based outlier detection
- Doublet detection (Scrublet + DoubletDetection)
- Cell and gene filtering
- QC visualization

Author: Alfred3005
"""

import logging
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import median_abs_deviation

try:
    import scrublet as scr
except ImportError:
    scr = None

try:
    import doubletdetection
except ImportError:
    doubletdetection = None

from .utils import setup_logging, log_memory_usage, cleanup_memory


def calculate_qc_metrics(
    adata: ad.AnnData,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Calculate comprehensive QC metrics for single-cell data.

    Calculates:
    - n_genes_by_counts: Number of genes detected per cell
    - total_counts: Total UMI counts per cell
    - pct_counts_mt: Percentage mitochondrial gene expression
    - pct_counts_ribo: Percentage ribosomal gene expression
    - pct_counts_hb: Percentage hemoglobin gene expression

    Parameters
    ----------
    adata : AnnData
        Input AnnData object with raw counts
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with QC metrics in .obs

    Examples
    --------
    >>> adata = calculate_qc_metrics(adata, logger)
    >>> print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']])
    """
    if logger is None:
        logger = setup_logging()

    logger.info("Calculating QC metrics...")

    # Identify gene types
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.match(r'^RP[SL]')
    adata.var['hb'] = adata.var_names.str.contains(r'^HB[^(P)]')

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo', 'hb'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Log summary statistics
    logger.info("QC metrics summary:")
    logger.info(f"  Median genes/cell: {np.median(adata.obs['n_genes_by_counts']):.0f}")
    logger.info(f"  Median counts/cell: {np.median(adata.obs['total_counts']):.0f}")
    logger.info(f"  Median MT%: {np.median(adata.obs['pct_counts_mt']):.2f}")
    logger.info(f"  Median Ribo%: {np.median(adata.obs['pct_counts_ribo']):.2f}")
    logger.info(f"  Median HB%: {np.median(adata.obs['pct_counts_hb']):.2f}")

    # Check for potential issues
    high_mt = (adata.obs['pct_counts_mt'] > 20).sum()
    if high_mt > 0:
        logger.warning(
            f"{high_mt} cells ({high_mt/adata.n_obs*100:.1f}%) have >20% MT genes"
        )

    high_hb = (adata.obs['pct_counts_hb'] > 10).sum()
    if high_hb > 0:
        logger.warning(
            f"{high_hb} cells ({high_hb/adata.n_obs*100:.1f}%) have >10% HB genes "
            f"(possible RBC contamination)"
        )

    return adata


def detect_outliers_mad(
    values: np.ndarray,
    n_mads: float = 5.0
) -> np.ndarray:
    """
    Detect outliers using Median Absolute Deviation (MAD) method.

    MAD is robust to outliers and works better than standard deviation
    for non-normal distributions.

    Parameters
    ----------
    values : np.ndarray
        Array of values to check for outliers
    n_mads : float, default 5.0
        Number of MADs from median to consider outlier

    Returns
    -------
    np.ndarray
        Boolean array where True indicates outlier

    Examples
    --------
    >>> outliers = detect_outliers_mad(adata.obs['total_counts'], n_mads=5)
    >>> print(f"Detected {outliers.sum()} outliers")
    """
    median = np.median(values)
    mad = median_abs_deviation(values)

    # Handle case where MAD is zero (all values identical)
    if mad == 0:
        return np.zeros(len(values), dtype=bool)

    lower_bound = median - n_mads * mad
    upper_bound = median + n_mads * mad

    outliers = (values < lower_bound) | (values > upper_bound)

    return outliers


def filter_cells_qc(
    adata: ad.AnnData,
    config: Dict[str, Any],
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Filter cells based on QC thresholds.

    Applies:
    1. MAD-based outlier detection for n_genes, total_counts, pct_mt
    2. Hard thresholds for min_genes, min_counts, max_pct_mt

    Parameters
    ----------
    adata : AnnData
        Input AnnData with QC metrics
    config : dict
        Configuration from qc_thresholds.yaml
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        Filtered AnnData

    Examples
    --------
    >>> config = load_config("config/qc_thresholds.yaml")
    >>> adata = filter_cells_qc(adata, config, logger)
    """
    if logger is None:
        logger = setup_logging()

    n_cells_before = adata.n_obs
    logger.info(f"Starting cell filtering: {n_cells_before:,} cells")

    # Extract thresholds from config
    qc_params = config['qc_metrics']
    outlier_params = config['outlier_detection']

    # Initialize filtering flags
    adata.obs['pass_qc'] = True

    # 1. Hard thresholds
    adata.obs['pass_min_genes'] = adata.obs['n_genes_by_counts'] >= qc_params['min_genes']
    adata.obs['pass_min_counts'] = adata.obs['total_counts'] >= qc_params['min_counts']
    adata.obs['pass_max_mt'] = adata.obs['pct_counts_mt'] <= qc_params['max_pct_mt']
    adata.obs['pass_max_ribo'] = adata.obs['pct_counts_ribo'] <= qc_params['max_pct_ribo']
    adata.obs['pass_max_hb'] = adata.obs['pct_counts_hb'] <= qc_params['max_pct_hb']

    # 2. MAD-based outlier detection
    n_mads = outlier_params['n_mads']

    adata.obs['outlier_genes'] = detect_outliers_mad(
        adata.obs['n_genes_by_counts'].values,
        n_mads=n_mads
    )

    adata.obs['outlier_counts'] = detect_outliers_mad(
        adata.obs['total_counts'].values,
        n_mads=n_mads
    )

    adata.obs['outlier_mt'] = detect_outliers_mad(
        adata.obs['pct_counts_mt'].values,
        n_mads=n_mads
    )

    # Combine all filters
    adata.obs['pass_qc'] = (
        adata.obs['pass_min_genes'] &
        adata.obs['pass_min_counts'] &
        adata.obs['pass_max_mt'] &
        adata.obs['pass_max_ribo'] &
        adata.obs['pass_max_hb'] &
        ~adata.obs['outlier_genes'] &
        ~adata.obs['outlier_counts'] &
        ~adata.obs['outlier_mt']
    )

    # Log filtering statistics
    logger.info("Cell filtering results:")
    logger.info(f"  Failing min_genes ({qc_params['min_genes']}): {(~adata.obs['pass_min_genes']).sum():,}")
    logger.info(f"  Failing min_counts ({qc_params['min_counts']}): {(~adata.obs['pass_min_counts']).sum():,}")
    logger.info(f"  Failing max_mt ({qc_params['max_pct_mt']}%): {(~adata.obs['pass_max_mt']).sum():,}")
    logger.info(f"  Outliers (genes): {adata.obs['outlier_genes'].sum():,}")
    logger.info(f"  Outliers (counts): {adata.obs['outlier_counts'].sum():,}")
    logger.info(f"  Outliers (MT%): {adata.obs['outlier_mt'].sum():,}")

    # Filter
    adata = adata[adata.obs['pass_qc']].copy()

    n_cells_after = adata.n_obs
    logger.info(
        f"Cells retained: {n_cells_after:,} / {n_cells_before:,} "
        f"({n_cells_after/n_cells_before*100:.1f}%)"
    )

    return adata


def detect_doublets_scrublet(
    adata: ad.AnnData,
    config: Dict[str, Any],
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Detect doublets using Scrublet (simulation-based method).

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    config : dict
        Configuration from qc_thresholds.yaml
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with 'doublet_score_scrublet' and 'predicted_doublet_scrublet' in .obs

    Examples
    --------
    >>> adata = detect_doublets_scrublet(adata, config, logger)
    """
    if logger is None:
        logger = setup_logging()

    if scr is None:
        logger.error("Scrublet not installed. Skipping scrublet doublet detection.")
        adata.obs['doublet_score_scrublet'] = 0.0
        adata.obs['predicted_doublet_scrublet'] = False
        return adata

    logger.info("Detecting doublets with Scrublet...")

    scrub_params = config['doublet_detection']['scrublet']

    # Run Scrublet
    scrub = scr.Scrublet(
        adata.X,
        expected_doublet_rate=scrub_params['expected_doublet_rate']
    )

    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=scrub_params['min_counts'],
        min_cells=scrub_params['min_cells'],
        min_gene_variability_pctile=scrub_params['min_gene_variability_pctile'],
        n_prin_comps=scrub_params['n_prin_comps']
    )

    # Add to adata
    adata.obs['doublet_score_scrublet'] = doublet_scores
    adata.obs['predicted_doublet_scrublet'] = predicted_doublets

    doublet_rate = predicted_doublets.sum() / len(predicted_doublets) * 100

    logger.info(
        f"Scrublet detected {predicted_doublets.sum():,} doublets "
        f"({doublet_rate:.2f}%)"
    )

    return adata


def detect_doublets_doubletdetection(
    adata: ad.AnnData,
    config: Dict[str, Any],
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Detect doublets using DoubletDetection (machine learning method).

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    config : dict
        Configuration from qc_thresholds.yaml
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with 'doublet_score_dd' and 'predicted_doublet_dd' in .obs

    Examples
    --------
    >>> adata = detect_doublets_doubletdetection(adata, config, logger)
    """
    if logger is None:
        logger = setup_logging()

    if doubletdetection is None:
        logger.error("DoubletDetection not installed. Skipping DD doublet detection.")
        adata.obs['doublet_score_dd'] = 0.0
        adata.obs['predicted_doublet_dd'] = False
        return adata

    logger.info("Detecting doublets with DoubletDetection...")

    dd_params = config['doublet_detection']['doubletdetection']

    # Run DoubletDetection
    clf = doubletdetection.BoostClassifier(
        boost_rate=dd_params['boost_rate'],
        n_iters=dd_params['n_iters'],
        use_phenograph=dd_params['use_phenograph'],
        standard_scaling=dd_params['standard_scaling']
    )

    doublet_scores = clf.fit(adata.X).predict(p_thresh=1e-7, voter_thresh=0.5)
    predicted_doublets = doublet_scores == 1

    # Add to adata
    adata.obs['doublet_score_dd'] = doublet_scores
    adata.obs['predicted_doublet_dd'] = predicted_doublets

    doublet_rate = predicted_doublets.sum() / len(predicted_doublets) * 100

    logger.info(
        f"DoubletDetection detected {predicted_doublets.sum():,} doublets "
        f"({doublet_rate:.2f}%)"
    )

    return adata


def remove_doublets(
    adata: ad.AnnData,
    config: Dict[str, Any],
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Remove doublets using consensus approach.

    If consensus=True in config, only removes cells predicted as doublets
    by both methods. Otherwise, removes cells predicted by either method.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with doublet predictions
    config : dict
        Configuration from qc_thresholds.yaml
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with doublets removed

    Examples
    --------
    >>> config = load_config("config/qc_thresholds.yaml")
    >>> adata = detect_doublets_scrublet(adata, config, logger)
    >>> adata = detect_doublets_doubletdetection(adata, config, logger)
    >>> adata = remove_doublets(adata, config, logger)
    """
    if logger is None:
        logger = setup_logging()

    n_cells_before = adata.n_obs

    consensus = config['doublet_detection']['consensus']

    # Check if doublet predictions exist
    has_scrublet = 'predicted_doublet_scrublet' in adata.obs.columns
    has_dd = 'predicted_doublet_dd' in adata.obs.columns

    if not (has_scrublet or has_dd):
        logger.warning("No doublet predictions found. Skipping doublet removal.")
        return adata

    # Determine doublets based on consensus setting
    if consensus and has_scrublet and has_dd:
        # Both methods must agree
        adata.obs['is_doublet'] = (
            adata.obs['predicted_doublet_scrublet'] &
            adata.obs['predicted_doublet_dd']
        )
        logger.info("Using consensus approach (both methods must agree)")

    elif has_scrublet and has_dd:
        # Either method
        adata.obs['is_doublet'] = (
            adata.obs['predicted_doublet_scrublet'] |
            adata.obs['predicted_doublet_dd']
        )
        logger.info("Using union approach (either method)")

    elif has_scrublet:
        adata.obs['is_doublet'] = adata.obs['predicted_doublet_scrublet']
        logger.info("Using Scrublet only")

    else:
        adata.obs['is_doublet'] = adata.obs['predicted_doublet_dd']
        logger.info("Using DoubletDetection only")

    # Remove doublets
    n_doublets = adata.obs['is_doublet'].sum()
    adata = adata[~adata.obs['is_doublet']].copy()

    n_cells_after = adata.n_obs

    logger.info(
        f"Doublets removed: {n_doublets:,} "
        f"({n_doublets/n_cells_before*100:.1f}%)"
    )
    logger.info(
        f"Cells retained: {n_cells_after:,} / {n_cells_before:,} "
        f"({n_cells_after/n_cells_before*100:.1f}%)"
    )

    return adata


def filter_genes(
    adata: ad.AnnData,
    min_cells: int = 3,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Filter genes expressed in fewer than min_cells.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    min_cells : int, default 3
        Minimum number of cells expressing a gene
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        Filtered AnnData

    Examples
    --------
    >>> adata = filter_genes(adata, min_cells=3, logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    n_genes_before = adata.n_vars

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_genes_after = adata.n_vars

    logger.info(
        f"Genes retained: {n_genes_after:,} / {n_genes_before:,} "
        f"({n_genes_after/n_genes_before*100:.1f}%) "
        f"[min_cells={min_cells}]"
    )

    return adata


def plot_qc_metrics(
    adata: ad.AnnData,
    output_dir: Path,
    groupby: Optional[str] = None,
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Generate QC metric visualization plots.

    Creates violin plots for:
    - n_genes_by_counts
    - total_counts
    - pct_counts_mt
    - pct_counts_ribo
    - pct_counts_hb

    Parameters
    ----------
    adata : AnnData
        Input AnnData with QC metrics
    output_dir : Path
        Output directory for plots
    groupby : str, optional
        Variable to group by (e.g., 'dataset_id')
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> plot_qc_metrics(adata, Path("results/figures"), groupby='dataset_id', logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    logger.info("Generating QC plots...")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = [
        'n_genes_by_counts',
        'total_counts',
        'pct_counts_mt',
        'pct_counts_ribo',
        'pct_counts_hb'
    ]

    # Violin plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, metric in enumerate(metrics):
        if metric not in adata.obs.columns:
            continue

        sc.pl.violin(
            adata,
            keys=metric,
            groupby=groupby,
            ax=axes[idx],
            show=False,
            rotation=45
        )
        axes[idx].set_title(metric.replace('_', ' ').title())

    # Remove empty subplot
    fig.delaxes(axes[5])

    plt.tight_layout()
    output_path = output_dir / "qc_violin_plots.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"QC plots saved to: {output_path}")

    # Scatter plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sc.pl.scatter(
        adata,
        x='total_counts',
        y='n_genes_by_counts',
        color='pct_counts_mt',
        ax=axes[0],
        show=False
    )
    axes[0].set_title('Counts vs Genes (colored by MT%)')

    sc.pl.scatter(
        adata,
        x='total_counts',
        y='pct_counts_mt',
        ax=axes[1],
        show=False
    )
    axes[1].set_title('Counts vs MT%')

    plt.tight_layout()
    output_path = output_dir / "qc_scatter_plots.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Scatter plots saved to: {output_path}")
