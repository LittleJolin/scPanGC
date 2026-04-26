import pandas as pd
import pytest

from scPanGC.validation import validate_obs_columns


class DummyAnnData:
    def __init__(self, columns):
        self.obs = pd.DataFrame(columns=columns)


def test_validate_obs_columns_accepts_present_columns():
    adata = DummyAnnData(["Celltype_new", "Disease2", "Tissue1"])

    validate_obs_columns(adata, ["Celltype_new", "Disease2", "Tissue1"])


def test_validate_obs_columns_reports_missing_columns():
    adata = DummyAnnData(["Celltype_new"])

    with pytest.raises(ValueError, match="Disease2"):
        validate_obs_columns(adata, ["Celltype_new", "Disease2"])
