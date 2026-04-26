from importlib import import_module

__version__ = "0.1.0"
__author__ = "Jolin"

__all__ = [
    "get_consensus_degs_from_raw",
    "run_metacell_pipeline",
    "compute_gc_modules",
    "run_pipeline",
    "DeaConfig",
    "MetacellConfig",
    "ClusteringConfig",
    "PipelineConfig",
]

_API_MODULES = {
    "get_consensus_degs_from_raw": ".dea",
    "run_metacell_pipeline": ".metacell",
    "compute_gc_modules": ".clustering",
    "run_pipeline": ".pipeline",
    "DeaConfig": ".config",
    "MetacellConfig": ".config",
    "ClusteringConfig": ".config",
    "PipelineConfig": ".config",
}


def __getattr__(name):
    if name in _API_MODULES:
        module = import_module(_API_MODULES[name], __name__)
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
