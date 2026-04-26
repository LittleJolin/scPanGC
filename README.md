# scPanGC
A computational framework for extracting Gene Clusters from large-scale single-cell autoimmune atlases

## Key Features
* **Strict Background Control**: Performs Differential Expression Analysis (DEA) explicitly within homologous tissue and cell-type constraints to eliminate tissue-specific confounders.
* **MetaCell Aggregation**: Effectively mitigates scRNA-seq sparsity and technical noise while preserving true biological heterogeneity.
* **Automated GC Extraction**: Calculates background-corrected $log_2(FC)$ matrices and executes robust hierarchical clustering to define co-expressed gene modules.

## Installation

You can install `scPanGC` directly from the source code.

```bash
# Clone the repository
git clone https://github.com/LittleJolin/scPanGC.git
cd scPanGC
# Install the package and its dependencies
pip install .
```

# Quick Start
```python
import scanpy as sc
import scPanGC as gc

## 1. Load your global raw AnnData object (with metadata)
### Ensure metadata includes 'Celltype_new', 'Disease2', and 'Tissue1'
adata_raw = sc.read_h5ad("path/to/your_raw_data.h5ad")

## 2. Extract Consensus DEGs
### Strictly compares Disease vs HC within the same tissue and cell type
deg_set = gc.get_consensus_degs_from_raw(
    adata_raw, 
    out_dir="./results", 
    min_cells=3, 
    logfc_thresh=1.0, 
    pval_thresh=0.05
)

## 3. Construct MetaCells for Noise Reduction
### Divides data into granular homologous subsets before building metacells
adata_mc = gc.run_metacell_pipeline(
    adata_raw, 
    output_path="./results/metacell_global.h5ad"
)

## 4. Compute Gene Clusters (GCs)
### Calculates expression shifts and performs Ward's hierarchical clustering
logfc_matrix, gc_modules = gc.compute_gc_modules(
    adata_mc, 
    deg_list=deg_set, 
    output_path="./results/GC_Modules.csv",
    distance_threshold=70.0
)

print("Pipeline completed successfully!")
```

## One-step API

```python
import scanpy as sc
import scPanGC as gc

adata_raw = sc.read_h5ad("path/to/your_raw_data.h5ad")
result = gc.run_pipeline(adata_raw, out_dir="./results")

deg_list = result["deg_list"]
gc_modules = result["gc_modules"]
```

## Command line

After installation, run the full workflow from a terminal:

```bash
scpangc path/to/your_raw_data.h5ad --out-dir results
```

## Development with uv

Create and sync a local development environment:

```bash
uv python install 3.11 --install-dir .uv-python
uv sync
```

Run lightweight checks:

```bash
uv run python -m py_compile scPanGC/*.py
uv run pytest
```

Run the CLI from the uv environment:

```bash
uv run scpangc path/to/your_raw_data.h5ad --out-dir results
```

## Required metadata

The raw AnnData object must contain these `obs` columns:

* `Celltype_new`
* `Disease2`
* `Tissue1`

The metacell AnnData object used for GC extraction must contain:

* `Celltype3`
* `Disease`
* `Tissue`
