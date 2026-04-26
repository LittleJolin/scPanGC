from dataclasses import dataclass, field
from typing import Tuple


@dataclass(frozen=True)
class DeaConfig:
    """Configuration for disease-versus-control differential expression."""

    min_cells: int = 3
    logfc_thresh: float = 1.0
    pval_thresh: float = 0.05
    control_label: str = "HC"
    exclude_gene_prefixes: Tuple[str, ...] = ("MT-", "RPS", "RPL")


@dataclass(frozen=True)
class MetacellConfig:
    """Configuration for metacell construction."""

    target_metacell_size: int = 50
    min_cells: int = 50
    random_seed: int = 123456
    compute_umap: bool = True
    max_top_feature_genes: int = 1000
    min_dist: float = 2.0
    required_gene_masks: Tuple[str, ...] = ("lateral_gene", "noisy_gene")


@dataclass(frozen=True)
class ClusteringConfig:
    """Configuration for gene-cluster extraction."""

    distance_threshold: float = 70.0
    pseudocount: float = 1.0
    control_label: str = "HC"


@dataclass(frozen=True)
class PipelineConfig:
    """Configuration for the complete scPanGC workflow."""

    dea: DeaConfig = field(default_factory=DeaConfig)
    metacell: MetacellConfig = field(default_factory=MetacellConfig)
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
