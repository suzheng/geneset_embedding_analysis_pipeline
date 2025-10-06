from __future__ import annotations

import argparse
from typing import List, Optional

from .pipeline import geneset_embedding


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Geneset embedding analysis pipeline (production-ready CLI)",
    )
    p.add_argument("--embeddings-path", required=True, help="Path to embeddings .txt file")
    p.add_argument("--entity-id", required=True, help="Unique identifier for this run")
    p.add_argument("--output-folder", required=False, default=None, help="Override output folder name")

    # Inputs (single sample)
    p.add_argument(
        "--input-genes",
        help="Comma-separated list of gene symbols for single-sample mode",
    )
    p.add_argument(
        "--input-genes-file",
        help="Path to text file with genes (one per line, or comma/space-separated)",
    )

    # Cohort metadata + precomputed terms
    p.add_argument("--metadata-csv", help="Metadata CSV/TSV with columns: sample_id, filepath, ...")
    p.add_argument("--precomputed-terms-jsonl", help="Precomputed term signatures JSONL path")
    p.add_argument("--grouping-column", help="Primary grouping column for cohort comparisons")
    p.add_argument("--subgroup-column", help="Optional secondary subgroup column (not used yet)")
    p.add_argument(
        "--comparison-prefix",
        default=None,
        help="Prefix to distinguish multiple comparison strategies",
    )

    # Steps selection
    p.add_argument(
        "--steps-to-run",
        default=None,
        help="Comma-separated steps (default: all). Options: preprocess,qc,signature,term_matching,cohort_integration,group_comparison",
    )

    # Signature and matching options
    p.add_argument("--dims-filter", default="quantile", help="Dimension filter: quantile|cv|zscore|topk|none")
    p.add_argument("--q-low", type=float, default=0.3, help="Low quantile threshold for std filtering")
    p.add_argument("--cv-threshold", type=float, default=1.0, help="CV threshold for cv filter")
    p.add_argument("--z-threshold", type=float, default=0.5, help="Z-score threshold for zscore filter (negative tail)")
    p.add_argument("--top-k", type=int, default=25, help="Top K results to keep for matches")

    p.add_argument("--enable-clustering", action="store_true", default=True, help="Enable optional clustering")
    p.add_argument("--no-clustering", dest="enable_clustering", action="store_false", help="Disable clustering")
    p.add_argument("--cluster-min-genes", type=int, default=20, help="Min genes to attempt clustering")
    p.add_argument("--random-state", type=int, default=17, help="Random seed for deterministic behavior")

    return p.parse_args()


def main() -> None:
    args = _parse_args()

    input_genes: Optional[List[str]] = None
    if args.input_genes:
        input_genes = [g.strip() for g in args.input_genes.split(",") if g.strip()]

    geneset_embedding(
        embeddings_path=args.embeddings_path,
        entity_id=args.entity_id,
        output_folder=args.output_folder,
        input_genes=input_genes,
        input_genes_file=args.input_genes_file,
        precomputed_terms_jsonl=args.precomputed_terms_jsonl,
        metadata_csv=args.metadata_csv,
        grouping_column=args.grouping_column,
        subgroup_column=args.subgroup_column,
        comparison_prefix=args.comparison_prefix,
        steps_to_run=args.steps_to_run,
        dims_filter=args.dims_filter,
        q_low=args.q_low,
        cv_threshold=args.cv_threshold,
        z_threshold=args.z_threshold,
        top_k=args.top_k,
        enable_clustering=args.enable_clustering,
        cluster_min_genes=args.cluster_min_genes,
        random_state=args.random_state,
    )


if __name__ == "__main__":
    main()


