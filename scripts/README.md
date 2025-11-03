# Scripts

Utility scripts for pipeline execution and data management.

## Available Scripts

### `run_pipeline.sh`

Master script to execute the complete analysis pipeline.

**Usage:**
```bash
# Full pipeline
bash scripts/run_pipeline.sh

# Test mode (5 datasets only)
bash scripts/run_pipeline.sh --test
```

**Features:**
- Runs all notebooks in sequence
- Checks dependencies
- Creates output directories
- Error handling and logging
- Colored terminal output

**Requirements:**
- Jupyter installed (`pip install jupyter`)
- All Python dependencies (`pip install -r requirements.txt`)

## Creating Additional Scripts

Place new scripts in this directory and make them executable:

```bash
chmod +x scripts/your_script.sh
```

## Future Scripts (Planned)

- `download_all_datasets.py` - Standalone dataset downloader
- `merge_datasets.py` - Merge processed datasets
- `generate_report.py` - Generate HTML summary report
- `benchmark_integration.py` - Compare integration methods
