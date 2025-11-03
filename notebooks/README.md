# Analysis Notebooks

This directory contains Jupyter notebooks for the NK Immunosenescence analysis pipeline.

## Execution Order

Run notebooks in numerical order:

1. **01_dataset_discovery.ipynb** - Dataset discovery and download from CELLxGENE
2. **02_exploratory_analysis.ipynb** - Initial data exploration
3. **03_quality_control.ipynb** - QC metrics and cell/gene filtering
4. **04_preprocessing.ipynb** - Normalization, HVG selection, PCA
5. **05_integration.ipynb** - Batch correction (scVI, scANVI, Harmony)
6. **06_differential_expression.ipynb** - Pseudobulk DE analysis
7. **07_biomarker_validation.ipynb** - GSEA and biological validation
8. **08_final_figures.ipynb** - Publication-ready visualizations

## Running Notebooks

### Individual Notebook

```bash
jupyter notebook notebooks/01_dataset_discovery.ipynb
```

### Automated Execution

```bash
# Run all notebooks sequentially
bash scripts/run_pipeline.sh

# Test mode (limited datasets)
bash scripts/run_pipeline.sh --test
```

### In Docker

```bash
# Start Jupyter Lab in Docker
docker run --gpus all -v $(pwd):/workspace -p 8888:8888 nk-immunosenescence:latest

# Access at: http://localhost:8888
```

## Expected Outputs

Each notebook generates:
- **Processed data**: Saved to `data/processed/`
- **Figures**: Saved to `results/figures/`
- **Tables**: Saved to `results/tables/`
- **Logs**: Saved to `results/reports/`

## Notes

- Notebooks use relative paths (`../data`, `../results`)
- All notebooks are designed to be run from the `notebooks/` directory
- Checkpoint files are saved to allow resuming from failures
- Memory usage is logged at each step

## Troubleshooting

**Kernel dies during integration:**
- Reduce batch size in scVI parameters (config/analysis_params.yaml)
- Process datasets in smaller batches

**Out of disk space:**
- Raw data (~10-50 GB) can be deleted after processing
- Keep only processed H5AD files for downstream analysis

**GPU out of memory:**
- Reduce `batch_size` in integration config
- Set `use_gpu=False` (slower but works)
