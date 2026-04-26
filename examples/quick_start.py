import scanpy as sc

import scPanGC as gc


adata_raw = sc.read_h5ad("path/to/your_raw_data.h5ad")
result = gc.run_pipeline(adata_raw, out_dir="./results")

print(result["gc_modules"].head())
