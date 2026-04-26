import logging
from pathlib import Path

import scanpy as sc

from .io import ensure_dir, write_gene_list
from .validation import RAW_OBS_COLUMNS, validate_obs_columns

logger = logging.getLogger(__name__)


def get_consensus_degs_from_raw(
    adata_raw,
    out_dir,
    min_cells=3,
    logfc_thresh=1.0,
    pval_thresh=0.05,
    control_label="HC",
    exclude_gene_prefixes=("MT-", "RPS", "RPL"),
):
    """
    Run disease-versus-control DEA within matched cell-type and tissue subsets.
    """
    validate_obs_columns(adata_raw, RAW_OBS_COLUMNS, context="adata_raw")
    logger.info("Starting differential expression analysis on raw single-cell matrix.")

    celltypes = adata_raw.obs["Celltype_new"].dropna().unique()
    tissues = adata_raw.obs["Tissue1"].dropna().unique()
    diseases = [disease for disease in adata_raw.obs["Disease2"].dropna().unique() if disease != control_label]

    all_deg_list = []

    for celltype in celltypes:
        adata_ct = adata_raw[adata_raw.obs["Celltype_new"] == celltype]
        for tissue in tissues:
            adata_ct_tissue = adata_ct[adata_ct.obs["Tissue1"] == tissue]
            cells_control = adata_ct_tissue[adata_ct_tissue.obs["Disease2"] == control_label].n_obs

            if cells_control < min_cells:
                continue

            for disease in diseases:
                cells_disease = adata_ct_tissue[adata_ct_tissue.obs["Disease2"] == disease].n_obs
                if cells_disease >= min_cells:
                    adata_test = adata_ct_tissue[
                        adata_ct_tissue.obs["Disease2"].isin([disease, control_label])
                    ].copy()
                    try:
                        sc.tl.rank_genes_groups(
                            adata_test,
                            groupby="Disease2",
                            groups=[disease],
                            reference=control_label,
                            method="wilcoxon",
                            use_raw=False,
                        )
                        result = sc.get.rank_genes_groups_df(adata_test, group=disease)
                        sig_genes = result[
                            (result["pvals_adj"] < pval_thresh)
                            & (result["logfoldchanges"] >= logfc_thresh)
                        ]["names"].tolist()
                        if sig_genes:
                            all_deg_list.extend(sig_genes)
                            logger.info(
                                "[%s | %s | %s vs %s] %s up-regulated DEGs.",
                                celltype,
                                tissue,
                                disease,
                                control_label,
                                len(sig_genes),
                            )
                    except Exception as error:
                        logger.warning(
                            "Skipped DEA for %s-%s-%s due to error: %s",
                            celltype,
                            tissue,
                            disease,
                            error,
                        )

    final_deg_set = sorted(set(all_deg_list))
    final_deg_set = [
        gene for gene in final_deg_set if not str(gene).startswith(tuple(exclude_gene_prefixes))
    ]
    logger.info("DEA completed. Total unique consensus up-regulated DEGs found: %s", len(final_deg_set))

    ensure_dir(out_dir)
    write_gene_list(final_deg_set, Path(out_dir) / "consensus_DEGs.csv")

    return final_deg_set
