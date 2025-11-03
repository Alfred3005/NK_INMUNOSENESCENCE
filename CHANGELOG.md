# Changelog

All notable changes to the NK_INMUNOSENESCENCE project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-01-03

### Added
- Initial project structure with modular architecture
- Comprehensive README with scientific background and usage instructions
- YAML-based configuration system for reproducibility
  - `datasets.yaml`: Dataset selection criteria
  - `qc_thresholds.yaml`: Quality control parameters
  - `analysis_params.yaml`: All analysis hyperparameters
- Docker environment optimized for NVIDIA RTX 4060 (8GB VRAM)
- Core Python modules:
  - `download.py`: CELLxGENE Census data retrieval
  - `qc.py`: Quality control and filtering
  - `preprocessing.py`: Normalization and feature selection
  - `integration.py`: Batch correction (scVI, scANVI, Harmony)
  - `differential_expression.py`: Pseudobulk DESeq2 analysis
  - `visualization.py`: Publication-quality plotting
  - `utils.py`: Helper functions and logging
- Jupyter notebooks for interactive analysis:
  - `01_dataset_discovery.ipynb`: Dataset search and download
  - `02_exploratory_analysis.ipynb`: Initial data exploration
  - `03_quality_control.ipynb`: QC metrics and filtering
  - `04_preprocessing.ipynb`: Normalization and HVG selection
  - `05_integration.ipynb`: Batch correction comparison
  - `06_differential_expression.ipynb`: Pseudobulk DE
  - `07_biomarker_validation.ipynb`: GSEA and validation
  - `08_final_figures.ipynb`: Publication figures
- Unit tests for critical functions
- Comprehensive .gitignore for data and results
- MIT License
- This CHANGELOG

### Project Specifications
- **Age Groups**: Young (<35), Adult (35-59), Aged (≥60)
- **Data Source**: CELLxGENE Census (stable release)
- **Inclusion**: Healthy controls only, peripheral blood, 10X platforms
- **Cell Type**: NK cells (using CELLxGENE annotations)
- **Integration**: scVI + scANVI + Harmony comparison
- **DE Method**: Pseudobulk with DESeq2
- **Validation**: Biological (GSEA, literature comparison)

### Technical Details
- Python 3.9+ with pinned dependencies
- GPU acceleration (CUDA 12.x support)
- Memory optimization for 64GB RAM, 8GB VRAM constraints
- Reproducible analysis with fixed random seeds
- Comprehensive logging and checkpointing

## [1.1.0] - 2025-01-03

### Updated
- **CRITICAL**: Updated to CELLxGENE Census v1.17.0+ API
- `src/download.py`: Now uses new `get_obs()` method (more efficient)
- Context manager pattern for census connections (best practice)
- Support for categorical columns (Census v1.13.0+ requirement)
- Census version now points to 2025-01-30 LTS (109M human cells)
- `requirements.txt`: cellxgene-census>=1.17.0

### Added
- `UPDATES_2025.md`: Comprehensive migration guide for Census API changes
- Fallback method in `discover_datasets()` for API compatibility
- Census version info logging (release date, total cells)
- Better filter construction for SOMA value filters

### Changed
- More efficient cell metadata querying
- Improved logging with census statistics
- Better error handling for API failures

### Fixed
- Compatibility with latest CELLxGENE Census schema (2.1.0)
- Categorical column handling in metadata

## [Unreleased]

### Planned Features
- Automated dataset update checker
- NK cell subtype analysis (CD56bright vs CD56dim)
- RNA velocity and pseudotime analysis
- Machine learning age predictor
- Interactive Shiny/Dash dashboard
- Integration with additional data repositories (GEO, ArrayExpress)
- Benchmarking against published aging signatures

---

**Note**: This project is under active development. Please report issues on [GitHub](https://github.com/Alfred3005/NK_INMUNOSENESCENCE/issues).
