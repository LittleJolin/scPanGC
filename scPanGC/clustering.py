import logging

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.cluster.hierarchy as sch
import scipy.stats as stats
from scipy.spatial.distance import pdist

from .config import ClusteringConfig
from .io import ensure_parent_dir
from .validation import METACELL_OBS_COLUMNS, validate_nonempty_gene_list, validate_obs_columns

logger = logging.getLogger(__name__)


def compute_gc_modules(
    adata_mc,
    deg_list,
    output_path,
    distance_threshold=70.0,
    pseudocount=1.0,
    control_label="HC",
):
    """
    Compute a log2FC matrix and derive gene clusters using Ward clustering.
    """
    config = ClusteringConfig(
        distance_threshold=distance_threshold,
        pseudocount=pseudocount,
        control_label=control_label,
    )
    validate_obs_columns(adata_mc, METACELL_OBS_COLUMNS, context="adata_mc")
    validate_nonempty_gene_list(deg_list, context="deg_list")

    logger.info("Starting GC module computation.")
    adata_mc.obs_names_make_unique()

    adata_mc.obs["C-D-T"] = (
        adata_mc.obs["Celltype3"].astype(str)
        + "-"
        + adata_mc.obs["Disease"].astype(str)
        + "-"
        + adata_mc.obs["Tissue"].astype(str)
    )

    overlap_genes = [gene for gene in deg_list if gene in adata_mc.var_names]
    validate_nonempty_gene_list(overlap_genes, context="overlap between deg_list and adata_mc.var_names")

    adata_sub = adata_mc[:, overlap_genes].copy()
    obs_df = sc.get.obs_df(adata_sub, keys=overlap_genes)
    obs_df["C-D-T"] = adata_sub.obs["C-D-T"].values
    avg_score = obs_df.groupby("C-D-T")[overlap_genes].mean() + config.pseudocount

    cdt_list = avg_score.index.tolist()
    df_logfc = pd.DataFrame(index=cdt_list, columns=overlap_genes, dtype="float32")

    for cdt in cdt_list:
        parts = str(cdt).split("-")
        if len(parts) < 3:
            df_logfc.loc[cdt] = 0.0
            continue

        celltype, _, tissue = parts[0], parts[1], parts[2]
        control_name = f"{celltype}-{config.control_label}-{tissue}"
        if control_name in avg_score.index:
            df_logfc.loc[cdt] = np.log2(avg_score.loc[cdt] / avg_score.loc[control_name])
        else:
            df_logfc.loc[cdt] = 0.0

    df_logfc = df_logfc.dropna(axis=1, how="all").fillna(0)
    validate_nonempty_gene_list(df_logfc.columns.tolist(), context="valid clustered genes")

    z_matrix = stats.zscore(df_logfc.values, axis=0, nan_policy="omit")
    z_matrix = np.nan_to_num(z_matrix)

    if df_logfc.shape[1] == 1:
        cluster_labels = np.array([1])
    else:
        distance_matrix = pdist(z_matrix.T, metric="euclidean")
        linkage_matrix = sch.linkage(distance_matrix, method="ward")
        cluster_labels = sch.fcluster(
            linkage_matrix,
            t=config.distance_threshold,
            criterion="distance",
        )

    gc_result = pd.DataFrame({"Gene": df_logfc.columns, "GC_Cluster": cluster_labels})
    ensure_parent_dir(output_path)
    gc_result.to_csv(output_path, index=False)

    logger.info("Generated %s gene clusters. Saved to %s", len(set(cluster_labels)), output_path)
    return df_logfc, gc_result
