"""
Utility Functions and Logging Configuration
============================================

This module provides helper functions for:
- Logging setup and management
- Configuration file loading
- Memory monitoring
- File I/O operations
- Random seed setting for reproducibility

Author: Alfred3005
"""

import logging
import os
import sys
import gc
import yaml
import psutil
import numpy as np
import random
from pathlib import Path
from typing import Dict, Any, Optional, Union
from datetime import datetime

import torch


def setup_logging(
    log_file: Optional[str] = None,
    log_level: str = "INFO",
    console_output: bool = True
) -> logging.Logger:
    """
    Configure logging for the pipeline.

    Parameters
    ----------
    log_file : str, optional
        Path to log file. If None, logs only to console.
    log_level : str, default "INFO"
        Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL
    console_output : bool, default True
        Whether to output logs to console

    Returns
    -------
    logging.Logger
        Configured logger instance

    Examples
    --------
    >>> logger = setup_logging(log_file="results/reports/pipeline.log")
    >>> logger.info("Pipeline started")
    """
    # Create logger
    logger = logging.getLogger("NK_Immunosenescence")
    logger.setLevel(getattr(logging, log_level.upper()))

    # Remove existing handlers
    logger.handlers = []

    # Create formatter
    formatter = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console handler
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, log_level.upper()))
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # File handler
    if log_file is not None:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(log_file), exist_ok=True)

        file_handler = logging.FileHandler(log_file, mode='a')
        file_handler.setLevel(getattr(logging, log_level.upper()))
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load YAML configuration file.

    Parameters
    ----------
    config_path : str or Path
        Path to YAML configuration file

    Returns
    -------
    dict
        Configuration dictionary

    Examples
    --------
    >>> config = load_config("config/datasets.yaml")
    >>> print(config['inclusion_criteria']['cell_types'])
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    return config


def set_random_seeds(seed: int = 42) -> None:
    """
    Set random seeds for reproducibility.

    Sets seeds for:
    - Python random module
    - NumPy
    - PyTorch (CPU and GPU)

    Parameters
    ----------
    seed : int, default 42
        Random seed value

    Examples
    --------
    >>> set_random_seeds(42)
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    # Set environment variable for Python hash seed
    os.environ['PYTHONHASHSEED'] = str(seed)


def get_memory_usage() -> Dict[str, float]:
    """
    Get current memory usage statistics.

    Returns
    -------
    dict
        Dictionary with memory statistics:
        - ram_used_gb: RAM used in GB
        - ram_available_gb: Available RAM in GB
        - ram_percent: RAM usage percentage
        - gpu_used_gb: GPU memory used in GB (if available)
        - gpu_total_gb: Total GPU memory in GB (if available)

    Examples
    --------
    >>> mem = get_memory_usage()
    >>> print(f"RAM usage: {mem['ram_used_gb']:.2f} GB ({mem['ram_percent']:.1f}%)")
    """
    # RAM statistics
    memory = psutil.virtual_memory()
    stats = {
        'ram_used_gb': memory.used / (1024 ** 3),
        'ram_available_gb': memory.available / (1024 ** 3),
        'ram_percent': memory.percent
    }

    # GPU statistics (if available)
    if torch.cuda.is_available():
        try:
            gpu_memory = torch.cuda.mem_get_info()
            stats['gpu_available_gb'] = gpu_memory[0] / (1024 ** 3)
            stats['gpu_total_gb'] = gpu_memory[1] / (1024 ** 3)
            stats['gpu_used_gb'] = (gpu_memory[1] - gpu_memory[0]) / (1024 ** 3)
        except:
            pass

    return stats


def log_memory_usage(logger: logging.Logger) -> None:
    """
    Log current memory usage.

    Parameters
    ----------
    logger : logging.Logger
        Logger instance

    Examples
    --------
    >>> logger = setup_logging()
    >>> log_memory_usage(logger)
    """
    mem = get_memory_usage()

    logger.info(
        f"Memory usage - RAM: {mem['ram_used_gb']:.2f} GB "
        f"({mem['ram_percent']:.1f}%), "
        f"Available: {mem['ram_available_gb']:.2f} GB"
    )

    if 'gpu_used_gb' in mem:
        logger.info(
            f"GPU memory - Used: {mem['gpu_used_gb']:.2f} GB, "
            f"Total: {mem['gpu_total_gb']:.2f} GB"
        )


def cleanup_memory(logger: Optional[logging.Logger] = None) -> None:
    """
    Perform garbage collection and clear GPU cache.

    Parameters
    ----------
    logger : logging.Logger, optional
        Logger instance for logging cleanup operations

    Examples
    --------
    >>> cleanup_memory(logger)
    """
    # Python garbage collection
    gc.collect()

    # PyTorch GPU cache
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    if logger is not None:
        logger.debug("Memory cleanup performed (gc + GPU cache clear)")


def ensure_dir(path: Union[str, Path]) -> Path:
    """
    Create directory if it doesn't exist.

    Parameters
    ----------
    path : str or Path
        Directory path

    Returns
    -------
    Path
        Path object of the directory

    Examples
    --------
    >>> ensure_dir("results/figures")
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_timestamp(format: str = "%Y%m%d_%H%M%S") -> str:
    """
    Get current timestamp as string.

    Parameters
    ----------
    format : str, default "%Y%m%d_%H%M%S"
        Datetime format string

    Returns
    -------
    str
        Formatted timestamp

    Examples
    --------
    >>> timestamp = get_timestamp()
    >>> print(f"Run started at: {timestamp}")
    """
    return datetime.now().strftime(format)


def check_gpu_available(logger: Optional[logging.Logger] = None) -> bool:
    """
    Check if GPU is available and log information.

    Parameters
    ----------
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    bool
        True if GPU is available, False otherwise

    Examples
    --------
    >>> has_gpu = check_gpu_available(logger)
    >>> if has_gpu:
    >>>     print("GPU acceleration enabled")
    """
    available = torch.cuda.is_available()

    if logger is not None:
        if available:
            gpu_name = torch.cuda.get_device_name(0)
            gpu_memory = torch.cuda.get_device_properties(0).total_memory / (1024 ** 3)
            logger.info(f"GPU available: {gpu_name} ({gpu_memory:.1f} GB)")
        else:
            logger.warning("GPU not available. Using CPU (slower).")

    return available


def validate_age_groups(age_series, age_group_config: Dict[str, Dict]) -> bool:
    """
    Validate that age groups are properly defined and non-overlapping.

    Parameters
    ----------
    age_series : pd.Series
        Series containing age values
    age_group_config : dict
        Configuration dictionary with age group definitions

    Returns
    -------
    bool
        True if age groups are valid

    Raises
    ------
    ValueError
        If age groups are overlapping or invalid

    Examples
    --------
    >>> age_config = {
    ...     'young': {'min': 18, 'max': 34},
    ...     'adult': {'min': 35, 'max': 59},
    ...     'aged': {'min': 60, 'max': 100}
    ... }
    >>> validate_age_groups(df['age'], age_config)
    """
    # Check for overlaps
    groups = list(age_group_config.values())

    for i in range(len(groups) - 1):
        if groups[i]['max'] >= groups[i + 1]['min']:
            raise ValueError(
                f"Age groups overlap: {groups[i]} and {groups[i + 1]}"
            )

    # Check that all ages fall into a group
    min_age = min(g['min'] for g in groups)
    max_age = max(g['max'] for g in groups)

    ages_out_of_range = age_series[
        (age_series < min_age) | (age_series > max_age)
    ]

    if len(ages_out_of_range) > 0:
        raise ValueError(
            f"{len(ages_out_of_range)} ages outside defined groups: "
            f"{ages_out_of_range.min()}-{ages_out_of_range.max()}"
        )

    return True


def save_checkpoint(
    obj: Any,
    filepath: Union[str, Path],
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Save checkpoint file (for recovery from interruptions).

    Parameters
    ----------
    obj : any
        Object to save (typically AnnData or dict)
    filepath : str or Path
        Path to save checkpoint
    logger : logging.Logger, optional
        Logger instance

    Examples
    --------
    >>> save_checkpoint(adata, "data/processed/checkpoint_qc.h5ad", logger)
    """
    filepath = Path(filepath)
    ensure_dir(filepath.parent)

    # For AnnData objects
    if hasattr(obj, 'write_h5ad'):
        obj.write_h5ad(filepath, compression='gzip')
    # For other objects (pickle)
    else:
        import pickle
        with open(filepath, 'wb') as f:
            pickle.dump(obj, f)

    if logger is not None:
        logger.info(f"Checkpoint saved: {filepath}")


def load_checkpoint(
    filepath: Union[str, Path],
    logger: Optional[logging.Logger] = None
) -> Any:
    """
    Load checkpoint file.

    Parameters
    ----------
    filepath : str or Path
        Path to checkpoint file
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    any
        Loaded object

    Examples
    --------
    >>> adata = load_checkpoint("data/processed/checkpoint_qc.h5ad", logger)
    """
    import anndata
    import pickle

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Checkpoint not found: {filepath}")

    # Try loading as H5AD first
    if filepath.suffix == '.h5ad':
        obj = anndata.read_h5ad(filepath)
    # Otherwise try pickle
    else:
        with open(filepath, 'rb') as f:
            obj = pickle.load(f)

    if logger is not None:
        logger.info(f"Checkpoint loaded: {filepath}")

    return obj
