# NK Immunosenescence Analysis Pipeline
# Optimized for NVIDIA RTX 4060 (8GB VRAM) with CUDA 12.x support
# Base: PyTorch with CUDA support

FROM nvcr.io/nvidia/pytorch:24.06-py3

# Metadata
LABEL maintainer="Alfred3005"
LABEL description="Single-cell RNA-seq analysis pipeline for NK cell immunosenescence"
LABEL version="1.0.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV CUDA_VISIBLE_DEVICES=0
ENV TORCH_CUDA_ARCH_LIST="8.9"

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    wget \
    vim \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libncursesw5-dev \
    xz-utils \
    tk-dev \
    libxml2-dev \
    libxmlsec1-dev \
    libffi-dev \
    liblzma-dev \
    libhdf5-dev \
    libcairo2-dev \
    pkg-config \
    # R dependencies for DESeq2 via rpy2
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install build tools
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install Python dependencies
# Copy requirements first for better Docker layer caching
COPY requirements.txt /tmp/requirements.txt

# Install JAX with CUDA 12 support explicitly
RUN pip install --no-cache-dir "jax[cuda12]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Install other Python packages
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Install R packages for DESeq2 (pseudobulk differential expression)
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'ComplexHeatmap'))"

# Install additional single-cell specific tools
RUN pip install --no-cache-dir \
    doubletdetection \
    harmonypy \
    bbknn

# Create necessary directories
RUN mkdir -p /workspace/data/raw \
    /workspace/data/processed \
    /workspace/data/metadata \
    /workspace/results/figures \
    /workspace/results/tables \
    /workspace/results/reports \
    /workspace/src \
    /workspace/notebooks \
    /workspace/scripts \
    /workspace/config \
    /workspace/tests

# Copy project files
COPY . /workspace/

# Set Python path
ENV PYTHONPATH="/workspace:${PYTHONPATH}"

# Configure Jupyter
RUN jupyter lab --generate-config && \
    echo "c.NotebookApp.ip = '0.0.0.0'" >> ~/.jupyter/jupyter_lab_config.py && \
    echo "c.NotebookApp.open_browser = False" >> ~/.jupyter/jupyter_lab_config.py && \
    echo "c.NotebookApp.port = 8888" >> ~/.jupyter/jupyter_lab_config.py && \
    echo "c.NotebookApp.allow_root = True" >> ~/.jupyter/jupyter_lab_config.py

# Configure matplotlib for non-interactive backend
RUN mkdir -p ~/.config/matplotlib && \
    echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

# Memory optimization for limited VRAM (8GB)
ENV PYTORCH_CUDA_ALLOC_CONF="max_split_size_mb:512"
ENV XLA_PYTHON_CLIENT_PREALLOCATE=false
ENV XLA_PYTHON_CLIENT_ALLOCATOR=platform

# Expose Jupyter port
EXPOSE 8888

# Default command: Start Jupyter Lab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
