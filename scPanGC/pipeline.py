from pathlib import Path

from .config import PipelineConfig


def run_pipeline(adata_raw, out_dir, config=None):
    """
    Run the full scPanGC workflow: DEG discovery, metacell construction, and GC extraction.
    """
    from .clustering import compute_gc_modules
    from .dea import get_consensus_degs_from_raw
    from .metacell import run_metacell_pipeline

    config = config or PipelineConfig()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    deg_list = get_consensus_degs_from_raw(
        adata_raw,
        out_dir=out_dir,
        min_cells=config.dea.min_cells,
        logfc_thresh=config.dea.logfc_thresh,
        pval_thresh=config.dea.pval_thresh,
        control_label=config.dea.control_label,
        exclude_gene_prefixes=config.dea.exclude_gene_prefixes,
    )
    adata_mc = run_metacell_pipeline(
        adata_raw,
        output_path=out_dir / "metacell_global.h5ad",
        target_metacell_size=config.metacell.target_metacell_size,
        min_cells=config.metacell.min_cells,
        random_seed=config.metacell.random_seed,
        compute_umap=config.metacell.compute_umap,
    )
    logfc_matrix, gc_modules = compute_gc_modules(
        adata_mc,
        deg_list=deg_list,
        output_path=out_dir / "GC_Modules.csv",
        distance_threshold=config.clustering.distance_threshold,
        pseudocount=config.clustering.pseudocount,
        control_label=config.clustering.control_label,
    )

    return {
        "deg_list": deg_list,
        "adata_mc": adata_mc,
        "logfc_matrix": logfc_matrix,
        "gc_modules": gc_modules,
    }
