RAW_OBS_COLUMNS = ("Celltype_new", "Disease2", "Tissue1")
METACELL_OBS_COLUMNS = ("Celltype3", "Disease", "Tissue")


def validate_obs_columns(adata, columns, context="adata"):
    """Raise a clear error if required AnnData obs columns are missing."""
    missing = [column for column in columns if column not in adata.obs.columns]
    if missing:
        missing_text = ", ".join(missing)
        required_text = ", ".join(columns)
        raise ValueError(
            f"{context} is missing required obs column(s): {missing_text}. "
            f"Required columns are: {required_text}."
        )


def validate_nonempty_gene_list(genes, context="gene list"):
    """Raise a clear error if a gene list is empty."""
    if not genes:
        raise ValueError(f"{context} is empty.")
