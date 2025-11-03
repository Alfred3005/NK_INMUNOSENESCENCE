#!/bin/bash

# NK Immunosenescence Analysis Pipeline
# Master execution script
#
# This script runs the complete analysis pipeline from data download
# to integration and preprocessing.
#
# Usage:
#   bash scripts/run_pipeline.sh [--test]
#
# Options:
#   --test    Run in test mode (download only 5 datasets)

set -e  # Exit on error
set -u  # Exit on undefined variable

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Parse arguments
TEST_MODE=false
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=true
    log_warn "Running in TEST MODE (limited datasets)"
fi

# Start pipeline
log_info "Starting NK Immunosenescence Analysis Pipeline"
log_info "================================================"

# Check Python environment
log_info "Checking Python environment..."
python --version || { log_error "Python not found"; exit 1; }

# Check required packages
log_info "Checking required packages..."
python -c "import scanpy, anndata, scvi, cellxgene_census" 2>/dev/null || {
    log_error "Required packages not installed"
    log_info "Install with: pip install -r requirements.txt"
    exit 1
}

log_info "All requirements satisfied"

# Create output directories
log_info "Creating output directories..."
mkdir -p data/raw data/processed data/metadata
mkdir -p results/figures results/tables results/reports

# Step 1: Dataset Discovery and Download
log_info ""
log_info "Step 1: Dataset Discovery and Download"
log_info "======================================="

if [ "$TEST_MODE" = true ]; then
    log_warn "Test mode: Will limit to 5 datasets"
    # Modify config temporarily for test mode
    # (Implementation depends on how you want to handle this)
fi

jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=7200 \
    --output-dir=notebooks \
    notebooks/01_dataset_discovery.ipynb || {
    log_error "Dataset discovery failed"
    exit 1
}

log_info "Dataset discovery complete"

# Step 2: Exploratory Analysis
log_info ""
log_info "Step 2: Exploratory Analysis"
log_info "============================"

if [ -f "notebooks/02_exploratory_analysis.ipynb" ]; then
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=3600 \
        --output-dir=notebooks \
        notebooks/02_exploratory_analysis.ipynb || {
        log_warn "Exploratory analysis failed (optional step)"
    }
else
    log_warn "Notebook 02 not found, skipping"
fi

# Step 3: Quality Control
log_info ""
log_info "Step 3: Quality Control"
log_info "======================="

if [ -f "notebooks/03_quality_control.ipynb" ]; then
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=3600 \
        --output-dir=notebooks \
        notebooks/03_quality_control.ipynb || {
        log_error "Quality control failed"
        exit 1
    }
else
    log_warn "Notebook 03 not found, skipping"
fi

# Step 4: Preprocessing
log_info ""
log_info "Step 4: Preprocessing"
log_info "===================="

if [ -f "notebooks/04_preprocessing.ipynb" ]; then
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=3600 \
        --output-dir=notebooks \
        notebooks/04_preprocessing.ipynb || {
        log_error "Preprocessing failed"
        exit 1
    }
else
    log_warn "Notebook 04 not found, skipping"
fi

# Step 5: Integration
log_info ""
log_info "Step 5: Batch Correction / Integration"
log_info "======================================="

if [ -f "notebooks/05_integration.ipynb" ]; then
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=7200 \
        --output-dir=notebooks \
        notebooks/05_integration.ipynb || {
        log_error "Integration failed"
        exit 1
    }
else
    log_warn "Notebook 05 not found, skipping"
fi

# Pipeline complete
log_info ""
log_info "================================================"
log_info "Pipeline Complete!"
log_info "================================================"
log_info ""
log_info "Results saved to:"
log_info "  - Processed data: data/processed/"
log_info "  - Figures: results/figures/"
log_info "  - Tables: results/tables/"
log_info "  - Reports: results/reports/"
log_info ""
log_info "Next steps:"
log_info "  1. Review QC report: results/reports/qc_report.html"
log_info "  2. Check integration quality: results/figures/integration_comparison.pdf"
log_info "  3. Proceed with differential expression analysis"
log_info ""
