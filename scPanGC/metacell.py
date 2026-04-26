import logging

import anndata as ad
import metacells as mc

from .config import MetacellConfig
from .io import ensure_parent_dir
from .validation import RAW_OBS_COLUMNS, validate_obs_columns

logger = logging.getLogger(__name__)


def run_metacell_pipeline(
    adata,
    output_path,
    target_metacell_size=50,
    min_cells=50,
    random_seed=123456,
    compute_umap=True,
):
    """
    Build metacells within homologous cell-type, disease, and tissue subsets.
    """
    config = MetacellConfig(
        target_metacell_size=target_metacell_size,
        min_cells=min_cells,
        random_seed=random_seed,
        compute_umap=compute_umap,
    )
    validate_obs_columns(adata, RAW_OBS_COLUMNS, context="adata")

    logger.info("Starting metacell construction.")
    celltypes = adata.obs["Celltype_new"].dropna().unique()
    diseases = adata.obs["Disease2"].dropna().unique()
    tissues = adata.obs["Tissue1"].dropna().unique()

    metacell_obj_list = []

    for celltype in celltypes:
        adata_celltype = adata[adata.obs["Celltype_new"] == celltype]
        for disease in diseases:
            adata_disease = adata_celltype[adata_celltype.obs["Disease2"] == disease]
            valid_tissues = [tissue for tissue in tissues if tissue in adata_disease.obs["Tissue1"].values]

            for tissue in valid_tissues:
                adata_subset = adata_disease[adata_disease.obs["Tissue1"] == tissue].copy()
                subset_name = f"{celltype}-{disease}-{tissue}"

                if adata_subset.n_obs <= config.min_cells:
                    logger.info("%s skipped for metacell construction (n_obs=%s).", subset_name, adata_subset.n_obs)
                    continue

                logger.info("Processing metacells for %s (n_obs=%s).", subset_name, adata_subset.n_obs)
                for required_mask in config.required_gene_masks:
                    if required_mask not in adata_subset.var.columns:
                        adata_subset.var[required_mask] = False

                mc.pl.divide_and_conquer_pipeline(
                    adata_subset,
                    target_metacell_size=config.target_metacell_size,
                    random_seed=config.random_seed,
                )
                metacells = mc.pl.collect_metacells(
                    adata_subset,
                    name=f"{subset_name}.metacells",
                    random_seed=config.random_seed,
                )

                if config.compute_umap:
                    try:
                        mc.pl.compute_umap_by_features(
                            metacells,
                            max_top_feature_genes=config.max_top_feature_genes,
                            min_dist=config.min_dist,
                            random_seed=config.random_seed,
                        )
                    except Exception as error:
                        logger.warning("UMAP computation skipped for %s: %s", subset_name, error)

                if "__name__" in metacells.uns:
                    del metacells.uns["__name__"]

                metacells.obs["Celltype_new"] = celltype
                metacells.obs["Disease2"] = disease
                metacells.obs["Tissue1"] = tissue
                metacells.obs["Celltype3"] = celltype
                metacells.obs["Disease"] = disease
                metacells.obs["Tissue"] = tissue
                metacell_obj_list.append(metacells)

    if not metacell_obj_list:
        raise ValueError("No metacells were generated. Please check input data quality and min_cells.")

    logger.info("Concatenating %s metacell chunks into a global atlas.", len(metacell_obj_list))
    adata_metacell_global = ad.concat(metacell_obj_list, join="inner")
    adata_metacell_global.obs_names_make_unique()

    ensure_parent_dir(output_path)
    adata_metacell_global.write_h5ad(output_path)
    logger.info("Global metacell construction completed. Saved to %s", output_path)

    return adata_metacell_global
