"""
Data Download Module - CELLxGENE Census Integration
====================================================

This module handles dataset discovery and download from the CELLxGENE Census.

Key functions:
- discover_datasets: Find datasets matching inclusion criteria
- download_dataset: Download individual H5AD files
- batch_download: Download multiple datasets with progress tracking
- extract_metadata: Parse dataset metadata

Author: Alfred3005
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
import numpy as np
from tqdm import tqdm

try:
    import cellxgene_census
except ImportError:
    raise ImportError(
        "cellxgene-census not installed. "
        "Install with: pip install cellxgene-census"
    )

import anndata as ad

from .utils import setup_logging, ensure_dir, log_memory_usage


def discover_datasets(
    config: Dict[str, Any],
    logger: Optional[logging.Logger] = None
) -> pd.DataFrame:
    """
    Discover datasets from CELLxGENE Census matching inclusion criteria.

    Parameters
    ----------
    config : dict
        Configuration dictionary from datasets.yaml with keys:
        - census_version: Census version to use
        - inclusion_criteria: Dict with cell_types, disease_status, etc.
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    pd.DataFrame
        DataFrame with metadata for matching datasets:
        - dataset_id
        - title
        - cell_count
        - tissue_general
        - disease
        - assay
        - nk_cell_count (estimated)

    Examples
    --------
    >>> config = load_config("config/datasets.yaml")
    >>> datasets_df = discover_datasets(config, logger)
    >>> print(f"Found {len(datasets_df)} matching datasets")
    """
    if logger is None:
        logger = setup_logging()

    logger.info("Starting dataset discovery from CELLxGENE Census")

    census_version = config.get('census_version', 'stable')
    criteria = config['inclusion_criteria']

    # Open census
    logger.info(f"Opening CELLxGENE Census (version: {census_version})")
    census = cellxgene_census.open_soma(census_version=census_version)

    try:
        # Get dataset metadata
        logger.info("Retrieving dataset metadata...")
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()

        logger.info(f"Total datasets in census: {len(datasets)}")

        # Filter by organism
        organism = criteria.get('organism', 'Homo sapiens')
        datasets = datasets[datasets['organism'] == organism]
        logger.info(f"After organism filter ({organism}): {len(datasets)}")

        # Get cell metadata for filtering
        logger.info("Querying cell metadata (this may take a few minutes)...")

        # Query expression data to get cell-level metadata
        query = census["census_data"]["homo_sapiens"]

        # Get obs (cell metadata)
        obs = query.obs.read(
            column_names=[
                "dataset_id",
                "cell_type",
                "disease",
                "tissue_general",
                "assay",
                "donor_id",
                "sex",
                "is_primary_data"
            ],
            value_filter="is_primary_data == True"
        ).concat().to_pandas()

        logger.info(f"Retrieved metadata for {len(obs):,} cells")

        # Filter by cell type (NK cells)
        cell_types = criteria['cell_types']
        nk_mask = obs['cell_type'].str.lower().str.contains(
            '|'.join([ct.lower() for ct in cell_types]),
            na=False
        )
        nk_cells = obs[nk_mask]

        logger.info(
            f"Found {len(nk_cells):,} NK cells "
            f"({len(nk_cells)/len(obs)*100:.2f}% of total)"
        )

        # Filter by disease status (healthy/control)
        disease_status = [ds.lower() for ds in criteria['disease_status']]
        disease_mask = nk_cells['disease'].str.lower().isin(disease_status)
        nk_cells = nk_cells[disease_mask]

        logger.info(
            f"After disease filter (healthy/control): {len(nk_cells):,} cells"
        )

        # Filter by tissue
        tissues = [t.lower() for t in criteria['tissue_general']]
        tissue_mask = nk_cells['tissue_general'].str.lower().isin(tissues)
        nk_cells = nk_cells[tissue_mask]

        logger.info(
            f"After tissue filter (blood): {len(nk_cells):,} cells"
        )

        # Filter by assay
        assays = [a.lower() for a in criteria['assay']]
        assay_mask = nk_cells['assay'].str.lower().isin(assays)
        nk_cells = nk_cells[assay_mask]

        logger.info(
            f"After assay filter (10X): {len(nk_cells):,} cells"
        )

        # Count NK cells per dataset
        nk_counts = nk_cells.groupby('dataset_id').size().reset_index(name='nk_cell_count')

        # Filter datasets with minimum NK cells
        min_nk_cells = criteria.get('min_nk_cells', 50)
        nk_counts = nk_counts[nk_counts['nk_cell_count'] >= min_nk_cells]

        logger.info(
            f"Datasets with ≥{min_nk_cells} NK cells: {len(nk_counts)}"
        )

        # Merge with dataset metadata
        matching_datasets = datasets.merge(nk_counts, on='dataset_id', how='inner')

        # Sort by NK cell count (descending)
        matching_datasets = matching_datasets.sort_values(
            'nk_cell_count',
            ascending=False
        ).reset_index(drop=True)

        # Select relevant columns
        output_cols = [
            'dataset_id',
            'collection_id',
            'dataset_title',
            'dataset_h5ad_path',
            'dataset_total_cell_count',
            'nk_cell_count'
        ]

        matching_datasets = matching_datasets[output_cols]

        logger.info(
            f"\n{'='*60}\n"
            f"Dataset Discovery Summary:\n"
            f"  Total matching datasets: {len(matching_datasets)}\n"
            f"  Total NK cells: {matching_datasets['nk_cell_count'].sum():,}\n"
            f"  Mean NK cells/dataset: {matching_datasets['nk_cell_count'].mean():.0f}\n"
            f"  Median NK cells/dataset: {matching_datasets['nk_cell_count'].median():.0f}\n"
            f"{'='*60}"
        )

        return matching_datasets

    finally:
        census.close()
        logger.info("Census connection closed")


def download_dataset(
    dataset_id: str,
    output_dir: Path,
    census_version: str = "stable",
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    Download a single dataset from CELLxGENE Census.

    Parameters
    ----------
    dataset_id : str
        Dataset ID to download
    output_dir : Path
        Output directory for H5AD file
    census_version : str, default "stable"
        Census version
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    Path
        Path to downloaded H5AD file

    Examples
    --------
    >>> filepath = download_dataset(
    ...     "abc123",
    ...     Path("data/raw"),
    ...     logger=logger
    ... )
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_dir)

    output_path = output_dir / f"{dataset_id}.h5ad"

    # Check if already downloaded
    if output_path.exists():
        logger.info(f"Dataset already exists: {output_path}")
        return output_path

    logger.info(f"Downloading dataset: {dataset_id}")

    try:
        # Download using CELLxGENE Census API
        cellxgene_census.download_source_h5ad(
            dataset_id,
            to_path=str(output_path),
            census_version=census_version
        )

        # Verify download
        if output_path.exists():
            file_size = output_path.stat().st_size / (1024 ** 2)  # MB
            logger.info(
                f"Successfully downloaded {dataset_id} "
                f"({file_size:.1f} MB) to {output_path}"
            )
        else:
            raise FileNotFoundError(f"Download failed: {output_path} not created")

        return output_path

    except Exception as e:
        logger.error(f"Failed to download {dataset_id}: {str(e)}")
        # Clean up partial download
        if output_path.exists():
            output_path.unlink()
        raise


def batch_download(
    datasets_df: pd.DataFrame,
    output_dir: Path,
    census_version: str = "stable",
    max_datasets: Optional[int] = None,
    logger: Optional[logging.Logger] = None
) -> Dict[str, Path]:
    """
    Download multiple datasets with progress tracking.

    Parameters
    ----------
    datasets_df : pd.DataFrame
        DataFrame with dataset_id column
    output_dir : Path
        Output directory
    census_version : str, default "stable"
        Census version
    max_datasets : int, optional
        Maximum number of datasets to download (for testing)
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    dict
        Dictionary mapping dataset_id to filepath

    Examples
    --------
    >>> datasets_df = discover_datasets(config, logger)
    >>> filepaths = batch_download(datasets_df, Path("data/raw"), logger=logger)
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_dir)

    # Limit datasets if specified
    if max_datasets is not None:
        datasets_df = datasets_df.head(max_datasets)

    logger.info(f"Starting batch download of {len(datasets_df)} datasets")

    downloaded = {}
    failed = []

    for idx, row in tqdm(
        datasets_df.iterrows(),
        total=len(datasets_df),
        desc="Downloading datasets"
    ):
        dataset_id = row['dataset_id']

        try:
            filepath = download_dataset(
                dataset_id,
                output_dir,
                census_version,
                logger
            )
            downloaded[dataset_id] = filepath

        except Exception as e:
            logger.error(f"Failed to download {dataset_id}: {str(e)}")
            failed.append(dataset_id)
            continue

    logger.info(
        f"\nDownload complete:\n"
        f"  Successfully downloaded: {len(downloaded)}\n"
        f"  Failed: {len(failed)}"
    )

    if failed:
        logger.warning(f"Failed dataset IDs: {', '.join(failed)}")

    return downloaded


def extract_nk_cells(
    adata: ad.AnnData,
    cell_type_key: str = "cell_type",
    nk_keywords: Optional[List[str]] = None,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Extract NK cells from AnnData object using cell type annotations.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    cell_type_key : str, default "cell_type"
        Column name in adata.obs containing cell type annotations
    nk_keywords : list of str, optional
        Keywords to identify NK cells. If None, uses default list.
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        Filtered AnnData containing only NK cells

    Examples
    --------
    >>> adata = ad.read_h5ad("data/raw/dataset123.h5ad")
    >>> nk_adata = extract_nk_cells(adata, logger=logger)
    >>> print(f"Extracted {nk_adata.n_obs} NK cells")
    """
    if logger is None:
        logger = setup_logging()

    if nk_keywords is None:
        nk_keywords = [
            "natural killer",
            "NK cell",
            "CD56",
            "CD16"
        ]

    logger.info(f"Extracting NK cells from {adata.n_obs:,} total cells")

    # Check if cell_type_key exists
    if cell_type_key not in adata.obs.columns:
        raise KeyError(
            f"Cell type key '{cell_type_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Create mask for NK cells
    nk_mask = adata.obs[cell_type_key].str.lower().str.contains(
        '|'.join([kw.lower() for kw in nk_keywords]),
        na=False
    )

    nk_adata = adata[nk_mask].copy()

    logger.info(
        f"Extracted {nk_adata.n_obs:,} NK cells "
        f"({nk_adata.n_obs/adata.n_obs*100:.2f}% of total)"
    )

    # Log NK cell subtypes if available
    if nk_adata.n_obs > 0:
        subtype_counts = nk_adata.obs[cell_type_key].value_counts()
        logger.info("NK cell subtypes:")
        for subtype, count in subtype_counts.items():
            logger.info(f"  {subtype}: {count:,} cells")

    return nk_adata


def add_age_groups(
    adata: ad.AnnData,
    age_key: str = "donor_age",
    age_group_config: Dict[str, Dict[str, int]] = None,
    logger: Optional[logging.Logger] = None
) -> ad.AnnData:
    """
    Add age group annotations to AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    age_key : str, default "donor_age"
        Column name in adata.obs containing age values
    age_group_config : dict, optional
        Age group configuration. If None, uses default:
        Young: <35, Adult: 35-59, Aged: ≥60
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    AnnData
        AnnData with 'age_group' column in obs

    Examples
    --------
    >>> adata = add_age_groups(adata, age_key="donor_age", logger=logger)
    >>> print(adata.obs['age_group'].value_counts())
    """
    if logger is None:
        logger = setup_logging()

    if age_group_config is None:
        age_group_config = {
            'young': {'min': 18, 'max': 34, 'label': 'Young'},
            'adult': {'min': 35, 'max': 59, 'label': 'Adult'},
            'aged': {'min': 60, 'max': 100, 'label': 'Aged'}
        }

    # Check if age_key exists
    if age_key not in adata.obs.columns:
        logger.warning(
            f"Age key '{age_key}' not found in adata.obs. "
            f"Trying alternative column names..."
        )

        # Try alternative age column names
        alt_names = ['age', 'Age', 'donor_age', 'development_stage']
        for alt_name in alt_names:
            if alt_name in adata.obs.columns:
                age_key = alt_name
                logger.info(f"Using '{age_key}' as age column")
                break
        else:
            raise KeyError(
                f"No age column found. Available columns: {list(adata.obs.columns)}"
            )

    # Extract age values (handle different formats)
    ages = adata.obs[age_key].copy()

    # Convert to numeric if needed
    if ages.dtype == 'object':
        # Try to extract numeric values from strings like "50-60" or "50 year"
        ages = ages.str.extract(r'(\d+)')[0].astype(float)

    # Add age group column
    adata.obs['age_group'] = 'Unknown'

    for group_name, group_config in age_group_config.items():
        mask = (ages >= group_config['min']) & (ages <= group_config['max'])
        adata.obs.loc[mask, 'age_group'] = group_config['label']

    # Log distribution
    group_counts = adata.obs['age_group'].value_counts()
    logger.info("Age group distribution:")
    for group, count in group_counts.items():
        logger.info(f"  {group}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")

    # Warn about unknown ages
    unknown_count = (adata.obs['age_group'] == 'Unknown').sum()
    if unknown_count > 0:
        logger.warning(
            f"{unknown_count} cells ({unknown_count/adata.n_obs*100:.1f}%) "
            f"have unknown/unassigned age group"
        )

    return adata


def save_dataset_catalog(
    datasets_df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Save dataset catalog to CSV file.

    Parameters
    ----------
    datasets_df : pd.DataFrame
        Dataset metadata DataFrame
    output_path : Path
        Output CSV file path
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> save_dataset_catalog(datasets_df, Path("data/metadata/catalog.csv"))
    """
    if logger is None:
        logger = setup_logging()

    ensure_dir(output_path.parent)

    datasets_df.to_csv(output_path, index=False)

    logger.info(f"Dataset catalog saved to: {output_path}")
    logger.info(f"  {len(datasets_df)} datasets")
    logger.info(f"  {datasets_df['nk_cell_count'].sum():,} total NK cells")
