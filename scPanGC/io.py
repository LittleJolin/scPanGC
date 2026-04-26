from pathlib import Path

import pandas as pd


def ensure_dir(path):
    """Create a directory if it does not already exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def ensure_parent_dir(path):
    """Create the parent directory for an output file path."""
    parent = Path(path).expanduser().parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)


def write_gene_list(genes, output_path, column="Gene"):
    """Write a one-column gene list CSV."""
    ensure_parent_dir(output_path)
    pd.DataFrame(list(genes), columns=[column]).to_csv(output_path, index=False)
