"""
Visualization Module
====================

This module provides publication-quality plotting functions for scRNA-seq analysis:
- QC plots
- UMAP embeddings
- Heatmaps
- Volcano plots
- Gene expression plots

Author: Alfred3005
"""

import logging
from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad

from .utils import setup_logging, ensure_dir


# Set publication-quality defaults
sc.set_figure_params(
    dpi=300,
    dpi_save=300,
    frameon=False,
    vector_friendly=True,
    fontsize=10,
    figsize=(6, 6),
    format='pdf'
)


def plot_umap_grid(
    adata: ad.AnnData,
    color_by: List[str],
    output_path: Path,
    ncols: int = 3,
    figsize: Tuple[int, int] = (18, 12),
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Plot UMAP embeddings colored by multiple variables in a grid.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with UMAP computed
    color_by : list of str
        List of variables to color by (from .obs or .var_names)
    output_path : Path
        Output file path
    ncols : int, default 3
        Number of columns in grid
    figsize : tuple, default (18, 12)
        Figure size
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> plot_umap_grid(
    ...     adata,
    ...     color_by=['dataset_id', 'age_group', 'leiden', 'NCAM1'],
    ...     output_path=Path("results/figures/umap_grid.pdf")
    ... )
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_path.parent)

    logger.info(f"Generating UMAP grid plot: {output_path}")

    # Calculate number of rows
    nrows = int(np.ceil(len(color_by) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)

    if nrows == 1:
        axes = axes.reshape(1, -1)

    axes = axes.flatten()

    for idx, var in enumerate(color_by):
        sc.pl.umap(
            adata,
            color=var,
            ax=axes[idx],
            show=False,
            frameon=False,
            legend_loc='right margin' if var in adata.obs.columns else 'on data',
            legend_fontsize=8
        )
        axes[idx].set_title(var.replace('_', ' ').title())

    # Remove empty subplots
    for idx in range(len(color_by), len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"UMAP grid saved to: {output_path}")


def plot_volcano(
    de_results: pd.DataFrame,
    lfc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
    output_path: Path,
    top_n_labels: int = 20,
    figsize: Tuple[int, int] = (10, 8),
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Generate volcano plot for differential expression results.

    Parameters
    ----------
    de_results : pd.DataFrame
        DE results with columns: log2FoldChange, padj, gene
    lfc_threshold : float, default 0.5
        Log2 fold change threshold for significance
    fdr_threshold : float, default 0.05
        FDR threshold for significance
    output_path : Path
        Output file path
    top_n_labels : int, default 20
        Number of top genes to label
    figsize : tuple, default (10, 8)
        Figure size
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> plot_volcano(
    ...     de_results,
    ...     output_path=Path("results/figures/volcano_young_vs_aged.pdf")
    ... )
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_path.parent)

    logger.info(f"Generating volcano plot: {output_path}")

    # Prepare data
    df = de_results.copy()

    # Add -log10(padj)
    df['-log10(padj)'] = -np.log10(df['padj'].clip(lower=1e-300))

    # Classify genes
    df['significant'] = (
        (df['padj'] < fdr_threshold) &
        (np.abs(df['log2FoldChange']) > lfc_threshold)
    )

    df['direction'] = 'Not significant'
    df.loc[
        df['significant'] & (df['log2FoldChange'] > 0),
        'direction'
    ] = f'Up (LFC > {lfc_threshold})'
    df.loc[
        df['significant'] & (df['log2FoldChange'] < 0),
        'direction'
    ] = f'Down (LFC < -{lfc_threshold})'

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    colors = {
        'Not significant': '#CCCCCC',
        f'Up (LFC > {lfc_threshold})': '#E74C3C',
        f'Down (LFC < -{lfc_threshold})': '#3498DB'
    }

    for direction, color in colors.items():
        subset = df[df['direction'] == direction]
        ax.scatter(
            subset['log2FoldChange'],
            subset['-log10(padj)'],
            c=color,
            label=f"{direction} (n={len(subset)})",
            alpha=0.6,
            s=10
        )

    # Add threshold lines
    ax.axhline(
        -np.log10(fdr_threshold),
        color='black',
        linestyle='--',
        linewidth=1,
        alpha=0.5,
        label=f'FDR = {fdr_threshold}'
    )
    ax.axvline(lfc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(-lfc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)

    # Label top genes
    if top_n_labels > 0:
        top_up = df[df['direction'].str.contains('Up')].nlargest(
            top_n_labels // 2,
            'log2FoldChange'
        )
        top_down = df[df['direction'].str.contains('Down')].nsmallest(
            top_n_labels // 2,
            'log2FoldChange'
        )

        for _, row in pd.concat([top_up, top_down]).iterrows():
            ax.text(
                row['log2FoldChange'],
                row['-log10(padj)'],
                row.name,  # Gene name (assuming index is gene names)
                fontsize=8,
                alpha=0.7
            )

    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10(Adjusted P-value)', fontsize=12)
    ax.set_title('Differential Expression Volcano Plot', fontsize=14)
    ax.legend(loc='upper right', frameon=True, fontsize=10)
    ax.grid(alpha=0.3, linestyle=':')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Volcano plot saved to: {output_path}")


def plot_heatmap_markers(
    adata: ad.AnnData,
    marker_genes: List[str],
    groupby: str,
    output_path: Path,
    figsize: Tuple[int, int] = (10, 12),
    cmap: str = 'RdBu_r',
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Plot heatmap of marker gene expression across groups.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    marker_genes : list of str
        List of genes to plot
    groupby : str
        Variable to group by (e.g., 'age_group')
    output_path : Path
        Output file path
    figsize : tuple, default (10, 12)
        Figure size
    cmap : str, default 'RdBu_r'
        Colormap
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> plot_heatmap_markers(
    ...     adata,
    ...     marker_genes=['NCAM1', 'FCGR3A', 'NKG7'],
    ...     groupby='age_group',
    ...     output_path=Path("results/figures/heatmap_markers.pdf")
    ... )
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_path.parent)

    logger.info(f"Generating heatmap: {output_path}")

    # Filter genes that exist in data
    genes_in_data = [g for g in marker_genes if g in adata.var_names]

    if len(genes_in_data) < len(marker_genes):
        missing = set(marker_genes) - set(genes_in_data)
        logger.warning(
            f"{len(missing)} genes not found in data: {', '.join(list(missing)[:5])}"
        )

    if len(genes_in_data) == 0:
        logger.error("No marker genes found in data. Skipping heatmap.")
        return

    # Create heatmap
    sc.pl.heatmap(
        adata,
        var_names=genes_in_data,
        groupby=groupby,
        cmap=cmap,
        dendrogram=True,
        figsize=figsize,
        show=False,
        save=None
    )

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Heatmap saved to: {output_path}")


def plot_integration_metrics(
    metrics_df: pd.DataFrame,
    output_path: Path,
    figsize: Tuple[int, int] = (12, 6),
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Plot integration quality metrics comparison.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        DataFrame with integration metrics (rows=methods, cols=metrics)
    output_path : Path
        Output file path
    figsize : tuple, default (12, 6)
        Figure size
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> metrics_df = pd.DataFrame({
    ...     'method': ['scVI', 'scANVI', 'Harmony'],
    ...     'iLISI': [2.5, 2.8, 2.3],
    ...     'cLISI': [1.2, 1.1, 1.4],
    ...     'ASW_batch': [0.1, 0.05, 0.15]
    ... }).set_index('method')
    >>> plot_integration_metrics(
    ...     metrics_df,
    ...     output_path=Path("results/figures/integration_metrics.pdf")
    ... )
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_path.parent)

    logger.info(f"Generating integration metrics plot: {output_path}")

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Batch correction metrics (higher is better)
    batch_metrics = [c for c in metrics_df.columns if c in ['iLISI', 'kBET']]
    if batch_metrics:
        metrics_df[batch_metrics].plot(
            kind='bar',
            ax=axes[0],
            color=['#3498DB', '#E74C3C']
        )
        axes[0].set_title('Batch Mixing Metrics (higher is better)', fontsize=12)
        axes[0].set_ylabel('Score', fontsize=10)
        axes[0].set_xlabel('Integration Method', fontsize=10)
        axes[0].legend(loc='best')
        axes[0].grid(alpha=0.3, axis='y')

    # Bio conservation metrics (varies by metric)
    bio_metrics = [c for c in metrics_df.columns if c in ['cLISI', 'ASW_label', 'NMI', 'ARI']]
    if bio_metrics:
        metrics_df[bio_metrics].plot(
            kind='bar',
            ax=axes[1],
            color=['#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C']
        )
        axes[1].set_title('Biological Conservation Metrics', fontsize=12)
        axes[1].set_ylabel('Score', fontsize=10)
        axes[1].set_xlabel('Integration Method', fontsize=10)
        axes[1].legend(loc='best')
        axes[1].grid(alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Integration metrics plot saved to: {output_path}")


def plot_age_distribution(
    adata: ad.AnnData,
    age_key: str = 'donor_age',
    age_group_key: str = 'age_group',
    output_path: Path,
    figsize: Tuple[int, int] = (12, 5),
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Plot age distribution across samples.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    age_key : str, default 'donor_age'
        Column with continuous age values
    age_group_key : str, default 'age_group'
        Column with age group labels
    output_path : Path
        Output file path
    figsize : tuple, default (12, 5)
        Figure size
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> plot_age_distribution(
    ...     adata,
    ...     output_path=Path("results/figures/age_distribution.pdf")
    ... )
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_path.parent)

    logger.info(f"Generating age distribution plot: {output_path}")

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Histogram of age
    if age_key in adata.obs.columns:
        ages = adata.obs[age_key].dropna()
        axes[0].hist(ages, bins=30, color='#3498DB', alpha=0.7, edgecolor='black')
        axes[0].set_xlabel('Age (years)', fontsize=12)
        axes[0].set_ylabel('Number of Cells', fontsize=12)
        axes[0].set_title('Age Distribution', fontsize=14)
        axes[0].grid(alpha=0.3, axis='y')

    # Bar plot of age groups
    if age_group_key in adata.obs.columns:
        age_group_counts = adata.obs[age_group_key].value_counts().sort_index()
        age_group_counts.plot(
            kind='bar',
            ax=axes[1],
            color=['#3498DB', '#E74C3C', '#F39C12'],
            edgecolor='black'
        )
        axes[1].set_xlabel('Age Group', fontsize=12)
        axes[1].set_ylabel('Number of Cells', fontsize=12)
        axes[1].set_title('Age Group Distribution', fontsize=14)
        axes[1].grid(alpha=0.3, axis='y')
        axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Age distribution plot saved to: {output_path}")
