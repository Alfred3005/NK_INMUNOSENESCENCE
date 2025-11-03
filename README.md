# NK_INMUNOSENESCENCE

**Identification of Aging Biomarkers in Natural Killer Cells Through Integrated Single-Cell Transcriptomic Analysis**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-available-blue.svg)](https://www.docker.com/)

---

## Overview

This project implements a robust, reproducible bioinformatics pipeline to identify aging-related biomarkers in Natural Killer (NK) cells using publicly available single-cell RNA sequencing (scRNA-seq) data from the CELLxGENE Census database.

NK cells are critical components of the innate immune system, and their functionality changes with age (immunosenescence). This pipeline systematically analyzes transcriptomic profiles of NK cells across different age groups to identify molecular signatures associated with aging.

### Key Features

- **Automated dataset discovery** from CELLxGENE Census with transparent inclusion criteria
- **Rigorous quality control** with MAD-based outlier detection
- **Batch effect correction** using state-of-the-art integration methods (scVI, scANVI, Harmony)
- **Pseudobulk differential expression** analysis with DESeq2
- **Comprehensive documentation** and logging for full reproducibility
- **GPU-accelerated** processing for efficient large-scale analysis
- **Modular codebase** with reusable functions and clear separation of concerns

---

## Scientific Rationale

### Background

Immunosenescence refers to the gradual deterioration of the immune system with aging, characterized by:
- Reduced immune surveillance capacity
- Increased susceptibility to infections
- Impaired vaccine responses
- Elevated chronic inflammation (inflammaging)

NK cells exhibit significant age-related functional changes:
- Altered cytotoxic capacity
- Changes in receptor expression patterns
- Modified cytokine production profiles
- Shifts in subset composition (CD56bright vs CD56dim)

### Objectives

1. **Identify transcriptomic biomarkers** associated with NK cell aging
2. **Characterize molecular pathways** altered during immunosenescence
3. **Establish reproducible methodology** for large-scale scRNA-seq meta-analysis
4. **Provide processed datasets** for downstream hypothesis-driven research

### Age Group Stratification

- **Young**: < 35 years
- **Adult**: 35-59 years
- **Aged**: ≥ 60 years

This stratification captures early adulthood (immune maturity), middle age (initial senescence), and older age (advanced immunosenescence).

---

## Project Structure

```
NK_INMUNOSENESCENCE/
├── README.md                          # This file
├── CHANGELOG.md                       # Version history and updates
├── LICENSE                            # MIT License
├── .gitignore                         # Git exclusions
├── requirements.txt                   # Python dependencies (pinned versions)
├── Dockerfile                         # Reproducible computing environment
│
├── config/                            # Configuration files (YAML)
│   ├── datasets.yaml                  # Dataset selection criteria
│   ├── qc_thresholds.yaml            # Quality control parameters
│   └── analysis_params.yaml          # Analysis hyperparameters
│
├── data/                              # Data directory (gitignored)
│   ├── raw/                          # Downloaded H5AD files from CELLxGENE
│   ├── processed/                    # Intermediate processing outputs
│   └── metadata/                     # Dataset metadata and provenance
│
├── src/                               # Source code modules
│   ├── __init__.py
│   ├── download.py                   # CELLxGENE data retrieval functions
│   ├── qc.py                         # Quality control and filtering
│   ├── preprocessing.py              # Normalization and feature selection
│   ├── integration.py                # Batch correction methods
│   ├── differential_expression.py    # Pseudobulk DE analysis
│   ├── visualization.py              # Publication-quality plotting
│   └── utils.py                      # Helper functions and logging
│
├── notebooks/                         # Analysis notebooks (numbered workflow)
│   ├── 01_dataset_discovery.ipynb    # Dataset search and download
│   ├── 02_exploratory_analysis.ipynb # Initial data exploration
│   ├── 03_quality_control.ipynb      # QC metrics and filtering
│   ├── 04_preprocessing.ipynb        # Normalization and HVG selection
│   ├── 05_integration.ipynb          # Batch correction comparison
│   ├── 06_differential_expression.ipynb  # Pseudobulk DE analysis
│   ├── 07_biomarker_validation.ipynb     # GSEA and literature validation
│   └── 08_final_figures.ipynb        # Publication-ready visualizations
│
├── scripts/                           # Standalone execution scripts
│   ├── run_pipeline.sh               # Master pipeline execution script
│   └── download_all_datasets.py      # Batch dataset download
│
├── results/                           # Analysis outputs (gitignored)
│   ├── figures/                      # PNG/PDF figures
│   ├── tables/                       # CSV/TSV result tables
│   └── reports/                      # HTML reports and logs
│
└── tests/                             # Unit tests
    ├── test_qc.py
    ├── test_preprocessing.py
    └── test_integration.py
```

---

## Installation

### Prerequisites

- **Operating System**: Linux (Ubuntu 20.04+ recommended) or WSL2
- **GPU**: NVIDIA GPU with CUDA support (recommended for scVI/scANVI)
  - Minimum: 8GB VRAM (tested on RTX 4060)
  - Driver: NVIDIA drivers 470.x or newer
- **RAM**: 32GB minimum, 64GB recommended
- **Storage**: 100GB+ free space for datasets and intermediate files
- **Software**: Docker with NVIDIA Container Toolkit

### Quick Start with Docker (Recommended)

```bash
# Clone the repository
git clone https://github.com/Alfred3005/NK_INMUNOSENESCENCE.git
cd NK_INMUNOSENESCENCE

# Build Docker image (includes CUDA support)
docker build -t nk-immunosenescence:latest .

# Run container with GPU support
docker run --gpus all -v $(pwd):/workspace -p 8888:8888 nk-immunosenescence:latest

# Inside container, start Jupyter
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```

### Manual Installation (Alternative)

```bash
# Create Python 3.9 environment
# Note: Conda environments have shown driver conflicts; Docker is preferred

pip install -r requirements.txt

# Verify GPU availability
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
```

---

## Usage

### 1. Dataset Discovery and Download

The first step identifies all CELLxGENE datasets matching our inclusion criteria:

**Inclusion Criteria:**
- Contains NK cells (annotated as "natural killer cell" or subsets)
- Age metadata available and documented
- Disease status: "normal", "healthy", or "control"
- Tissue: Primarily peripheral blood (PBMC, whole blood)
- Sequencing platform: 10X Chromium (for consistency)

```bash
# Open notebook 01
jupyter notebook notebooks/01_dataset_discovery.ipynb

# Or run automated download script
python scripts/download_all_datasets.py --output data/raw --log results/reports/download_log.txt
```

**Expected Output:**
- `data/metadata/dataset_catalog.csv` - Metadata for all discovered datasets
- `data/raw/*.h5ad` - Downloaded H5AD files
- `results/reports/download_summary.html` - Summary report

### 2. Quality Control

Applies rigorous per-dataset QC before integration:

**QC Metrics Calculated:**
- Total UMI counts per cell
- Number of genes detected per cell
- Percentage mitochondrial gene expression
- Percentage ribosomal gene expression
- Percentage hemoglobin contamination
- Doublet scores (Scrublet + DoubletDetection)

**Filtering Strategy:**
- MAD-based outlier detection (5 MAD threshold)
- Mitochondrial content < 20%
- Minimum 200 genes/cell, 500 UMI/cell
- Doublet removal (consensus approach)

```bash
jupyter notebook notebooks/03_quality_control.ipynb
```

### 3. Preprocessing

Normalization and feature selection:

**Normalization:**
- Library size normalization (target: 10,000 counts/cell)
- Log-transformation: log(count + 1)
- Scaling: z-score with clipping at ±10 SD

**Feature Selection:**
- Highly variable genes (HVGs): 2,000 genes
- Batch-aware selection using Seurat v3 method
- Force inclusion of NK marker genes (NCAM1, FCGR3A, NKG7, etc.)
- Force inclusion of aging-related genes (CDKN2A, TP53, IL6, etc.)

```bash
jupyter notebook notebooks/04_preprocessing.ipynb
```

### 4. Batch Correction / Integration

Comparison of integration methods:

**Methods Implemented:**
1. **scVI** (single-cell Variational Inference)
   - Deep generative model for UMI count data
   - Accounts for batch effects via latent variables
   - GPU-accelerated training

2. **scANVI** (semi-supervised scVI)
   - Leverages cell type annotations
   - Improves biological signal preservation

3. **Harmony** (Fast integration baseline)
   - PCA-based correction
   - Computationally efficient

**Integration Quality Metrics:**
- Batch mixing: iLISI, ASW batch, kBET
- Biological conservation: cLISI, ASW cell type, NMI
- Overall: Silhouette score

```bash
jupyter notebook notebooks/05_integration.ipynb
```

### 5. Differential Expression Analysis

Pseudobulk aggregation + DESeq2 for robust DE testing:

**Workflow:**
1. Aggregate counts by donor × age group
2. Create pseudobulk samples (biological replicates)
3. DESeq2 modeling with covariates (sex, dataset, sequencing depth)
4. Multiple testing correction (Benjamini-Hochberg FDR < 0.05)
5. Effect size filtering (|log2FC| > 0.5)

**Comparisons:**
- Young vs Adult
- Young vs Aged
- Adult vs Aged

```bash
jupyter notebook notebooks/06_differential_expression.ipynb
```

### 6. Biomarker Validation

Biological validation of identified markers:

**Validation Strategies:**
1. **Gene Set Enrichment Analysis (GSEA)**
   - KEGG pathways
   - GO Biological Processes
   - Reactome pathways
   - Hallmark gene sets

2. **Literature Comparison**
   - Cross-reference with known aging markers
   - PubMed literature mining
   - Immune aging databases

3. **Functional Scoring**
   - Calculate aging signature scores per cell
   - Correlate with chronological age
   - Test predictive capacity

```bash
jupyter notebook notebooks/07_biomarker_validation.ipynb
```

---

## Configuration

All analysis parameters are externalized to YAML configuration files for transparency and reproducibility.

### `config/datasets.yaml`

```yaml
inclusion_criteria:
  cell_types:
    - "natural killer cell"
    - "CD56-positive, CD16-positive NK cell"
    - "CD56-positive, CD16-negative NK cell"

  disease_status:
    - "normal"
    - "healthy"
    - "control"

  tissue_general:
    - "blood"
    - "peripheral blood"
    - "PBMC"

  assay:
    - "10x 3' v2"
    - "10x 3' v3"
    - "10x 5' v2"

  age_required: true
  min_nk_cells: 50  # Minimum NK cells per dataset
```

### `config/qc_thresholds.yaml`

```yaml
qc_metrics:
  min_genes: 200
  min_counts: 500
  max_pct_mt: 20
  max_pct_ribo: 50

doublet_detection:
  methods:
    - scrublet
    - doubletdetection
  consensus: true  # Both methods must agree

outlier_detection:
  method: "MAD"  # Median Absolute Deviation
  n_mads: 5
```

### `config/analysis_params.yaml`

```yaml
normalization:
  target_sum: 10000
  log_transform: true
  scale: true
  max_value: 10  # Clip scaled values

feature_selection:
  n_top_genes: 2000
  flavor: "seurat_v3"
  batch_key: "dataset_id"

integration:
  scvi:
    n_latent: 30
    n_layers: 2
    dropout_rate: 0.1
    gene_likelihood: "nb"
    max_epochs: 400
    early_stopping: true
    early_stopping_patience: 15

  scanvi:
    max_epochs: 150
    early_stopping_patience: 10

  harmony:
    max_iter_harmony: 10

differential_expression:
  method: "pseudobulk_deseq2"
  min_cells_per_sample: 10
  fdr_threshold: 0.05
  lfc_threshold: 0.5
```

---

## Computational Requirements

### Memory Optimization Strategies

Given the 64GB RAM constraint, the pipeline implements:

1. **Chunked Processing**: Process datasets in batches
2. **Sparse Matrices**: Use CSR format throughout
3. **Garbage Collection**: Explicit memory cleanup between steps
4. **Checkpointing**: Save intermediate results to disk
5. **GPU Offloading**: Move scVI computations to GPU (8GB VRAM)

### Expected Resource Usage

| Step | RAM Usage | GPU VRAM | Time (estimated) |
|------|-----------|----------|------------------|
| Download | < 5GB | - | 30-60 min |
| QC | 15-25GB | - | 10-20 min |
| Preprocessing | 20-35GB | - | 15-30 min |
| scVI Integration | 25-45GB | 6-8GB | 1-3 hours |
| Pseudobulk DE | 10-20GB | - | 20-40 min |

**Note**: Times depend on total number of cells downloaded.

---

## Output Files

### Key Result Files

1. **Processed Data Objects**
   - `data/processed/nk_cells_qc.h5ad` - Post-QC NK cells
   - `data/processed/nk_cells_integrated_scvi.h5ad` - scVI integrated
   - `data/processed/nk_cells_integrated_scanvi.h5ad` - scANVI integrated
   - `data/processed/nk_cells_integrated_harmony.h5ad` - Harmony integrated

2. **Differential Expression Results**
   - `results/tables/de_young_vs_aged.csv` - DE results
   - `results/tables/biomarkers_fdr05.csv` - Significant biomarkers

3. **Metadata and Logs**
   - `data/metadata/dataset_catalog.csv` - Dataset provenance
   - `data/metadata/sample_metadata.csv` - Per-sample metadata
   - `results/reports/qc_report.html` - QC summary report
   - `results/reports/integration_metrics.csv` - Integration quality

4. **Figures**
   - `results/figures/qc_violin_plots.pdf`
   - `results/figures/umap_integration.pdf`
   - `results/figures/de_volcano_plot.pdf`
   - `results/figures/biomarker_heatmap.pdf`

---

## Reproducibility

### Version Control

All software versions are pinned in `requirements.txt`:

```
scanpy==1.10.2
scvi-tools==1.1.5
anndata==0.10.8
...
```

### Random Seeds

All stochastic operations use fixed random seeds:
- NumPy: `np.random.seed(42)`
- PyTorch: `torch.manual_seed(42)`
- scVI: `scvi.settings.seed = 42`

### Logging

Every analysis step generates detailed logs:
- Timestamps
- Parameters used
- Input/output file paths
- QC statistics
- Warnings and errors

Logs are saved to `results/reports/pipeline_log.txt`.

### Docker Image

The Dockerfile specifies exact base image and package versions:

```dockerfile
FROM nvcr.io/nvidia/pytorch:24.06-py3
# All dependencies installed with pinned versions
```

Build and tag images for version tracking:
```bash
docker build -t nk-immunosenescence:v1.0.0 .
```

---

## Troubleshooting

### Common Issues

**1. CUDA Out of Memory Error**

If scVI training fails with OOM:
```python
# In notebooks, reduce batch size
model.train(batch_size=256)  # Default is 512
```

Or process datasets in smaller batches.

**2. Docker Permission Errors**

```bash
# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker
```

**3. Kernel Dies During Integration**

Reduce memory footprint:
```python
# In integration notebook
adata.X = adata.X.astype('float32')  # Use float32 instead of float64
```

**4. Slow Download Speeds**

CELLxGENE Census uses AWS S3. If slow:
```bash
# Configure boto3 for multipart downloads
export AWS_METADATA_SERVICE_TIMEOUT=5
export AWS_METADATA_SERVICE_NUM_ATTEMPTS=0
```

---

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{nk_immunosenescence2025,
  author = {Alfred3005},
  title = {NK_INMUNOSENESCENCE: Identification of Aging Biomarkers in Natural Killer Cells},
  year = {2025},
  url = {https://github.com/Alfred3005/NK_INMUNOSENESCENCE}
}
```

### Data Sources

This project uses data from the CELLxGENE Census:

> CZ CELLxGENE Discover: A single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *bioRxiv* (2023). https://doi.org/10.1101/2023.10.30.563174

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -m 'Add feature'`)
4. Push to branch (`git push origin feature/improvement`)
5. Open a Pull Request

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Project Maintainer**: Alfred3005

- GitHub: [@Alfred3005](https://github.com/Alfred3005)
- Repository: [NK_INMUNOSENESCENCE](https://github.com/Alfred3005/NK_INMUNOSENESCENCE)

For questions, issues, or suggestions, please open an [issue](https://github.com/Alfred3005/NK_INMUNOSENESCENCE/issues) on GitHub.

---

## Acknowledgments

- **CELLxGENE Team** for providing curated single-cell data
- **scVI-tools developers** for deep learning integration methods
- **Scanpy community** for single-cell analysis ecosystem
- **Bioinformatics community** for open-source tools and best practices

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and updates.

---

**Last Updated**: 2025-01-03
**Version**: 1.0.0
