# Usage

`scPanGC` can be used as three separate steps or as a one-step pipeline.

## Python API

```python
import scanpy as sc
import scPanGC as gc

adata_raw = sc.read_h5ad("path/to/your_raw_data.h5ad")
result = gc.run_pipeline(adata_raw, out_dir="./results")
```

## CLI

```bash
scpangc path/to/your_raw_data.h5ad --out-dir results
```

## Input metadata

Raw input requires `Celltype`, `Disease`, and `Tissue` in `adata.obs`.
Gene-cluster extraction from a metacell object requires `Celltype`, `Disease`, and `Tissue`.
