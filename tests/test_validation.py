import pandas as pd
import pytest

from scPanGC.validation import validate_obs_columns


class DummyAnnData:
    def __init__(self, columns):
        self.obs = pd.DataFrame(columns=columns)


def test_validate_obs_columns_accepts_present_columns():
    adata = DummyAnnData(["Celltype", "Disease", "Tissue"])

    validate_obs_columns(adata, ["Celltype", "Disease", "Tissue"])


def test_validate_obs_columns_reports_missing_columns():
    adata = DummyAnnData(["Celltype"])

    with pytest.raises(ValueError, match="Disease"):
        validate_obs_columns(adata, ["Celltype", "Disease"])
