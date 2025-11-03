"""
NK Immunosenescence Analysis Pipeline
======================================

A comprehensive bioinformatics pipeline for identifying aging biomarkers
in Natural Killer cells using single-cell RNA sequencing data.

Modules:
--------
- download: Data retrieval from CELLxGENE Census
- qc: Quality control and filtering
- preprocessing: Normalization and feature selection
- integration: Batch correction methods
- differential_expression: Pseudobulk DE analysis
- visualization: Publication-quality plotting
- utils: Helper functions and logging

Author: Alfred3005
License: MIT
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Alfred3005"

from . import utils
from . import download
from . import qc
from . import preprocessing
from . import integration
from . import visualization

__all__ = [
    "utils",
    "download",
    "qc",
    "preprocessing",
    "integration",
    "visualization",
]

# Version info
__version__ = "1.0.0"
__author__ = "Alfred3005"
__description__ = "NK Immunosenescence Analysis Pipeline"
