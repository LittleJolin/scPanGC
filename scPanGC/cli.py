import argparse
import logging


def build_parser():
    parser = argparse.ArgumentParser(prog="scpangc", description="Run the scPanGC workflow.")
    parser.add_argument("input", help="Input raw AnnData .h5ad file.")
    parser.add_argument("-o", "--out-dir", default="results", help="Output directory.")
    parser.add_argument("--min-cells", type=int, default=3, help="Minimum cells for DEA comparisons.")
    parser.add_argument("--logfc-thresh", type=float, default=1.0, help="Minimum log fold-change for DEGs.")
    parser.add_argument("--pval-thresh", type=float, default=0.05, help="Adjusted p-value threshold for DEGs.")
    parser.add_argument("--distance-threshold", type=float, default=70.0, help="Hierarchical clustering cut threshold.")
    return parser


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
    args = build_parser().parse_args(argv)

    import scanpy as sc

    from .pipeline import run_pipeline

    adata_raw = sc.read_h5ad(args.input)
    from .config import ClusteringConfig, DeaConfig, PipelineConfig

    config = PipelineConfig(
        dea=DeaConfig(
            min_cells=args.min_cells,
            logfc_thresh=args.logfc_thresh,
            pval_thresh=args.pval_thresh,
        ),
        clustering=ClusteringConfig(distance_threshold=args.distance_threshold),
    )
    run_pipeline(adata_raw, args.out_dir, config=config)


if __name__ == "__main__":
    main()
