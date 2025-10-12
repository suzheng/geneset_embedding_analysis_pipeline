from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import List, Optional

from .pipeline import geneset_embedding


def _read_library_groups(selected_libraries_file: str) -> dict[str, List[str]]:
    """
    Read selected_libraries.txt and return a mapping of Group_Clean -> list of library filenames.
    """
    group_to_libraries: dict[str, List[str]] = {}
    with open(selected_libraries_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            library_filename = row["Library_Filename"]
            group_clean = row["Group_Clean"]
            if group_clean not in group_to_libraries:
                group_to_libraries[group_clean] = []
            group_to_libraries[group_clean].append(library_filename)
    return group_to_libraries


def _get_libraries_for_group(
    group_name: str, 
    selected_libraries_file: str
) -> List[str]:
    """
    Get list of library filenames for a specific group or all groups.
    
    Args:
        group_name: Either a specific group name or "all"
        selected_libraries_file: Path to selected_libraries.txt
        
    Returns:
        List of library filenames
    """
    group_to_libraries = _read_library_groups(selected_libraries_file)
    
    if group_name.lower() == "all":
        # Return all libraries from all groups
        all_libraries = []
        for libraries in group_to_libraries.values():
            all_libraries.extend(libraries)
        return all_libraries
    else:
        # Return libraries for specific group
        if group_name not in group_to_libraries:
            available = ", ".join(sorted(group_to_libraries.keys()))
            raise ValueError(
                f"Group '{group_name}' not found in {selected_libraries_file}. "
                f"Available groups: {available}"
            )
        return group_to_libraries[group_name]


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
    
    # Reference genes (optional, for differential embedding)
    p.add_argument(
        "--reference-genes",
        help="Comma-separated list of reference gene symbols (computes input - reference embedding)",
    )
    p.add_argument(
        "--reference-genes-file",
        help="Path to text file with reference genes (one per line, or comma/space-separated)",
    )

    # Precomputed terms
    p.add_argument("--precomputed-terms-jsonl", help="Precomputed term signatures JSONL path (single file)")
    p.add_argument(
        "--library-group", 
        help="Library group name from Group_Clean column in selected_libraries.txt, or 'all' for all groups"
    )
    p.add_argument(
        "--selected-libraries-file",
        default=None,
        help="Path to selected_libraries.txt file (default: scripts/selected_libraries.txt)"
    )
    p.add_argument(
        "--precomputed-dir",
        default=None,
        help="Directory containing precomputed signatures (default: precomputed/)"
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
    
    reference_genes: Optional[List[str]] = None
    if args.reference_genes:
        reference_genes = [g.strip() for g in args.reference_genes.split(",") if g.strip()]

    # Determine if we're running with library groups or single file
    if args.library_group:
        # Running with library groups
        # Set up default paths
        project_root = Path(__file__).resolve().parent.parent
        
        if args.selected_libraries_file:
            selected_libraries_file = args.selected_libraries_file
        else:
            selected_libraries_file = str(project_root / "scripts" / "selected_libraries.txt")
        
        if args.precomputed_dir:
            precomputed_dir = Path(args.precomputed_dir)
        else:
            precomputed_dir = project_root / "precomputed"
        
        # Get list of libraries for the specified group
        libraries = _get_libraries_for_group(args.library_group, selected_libraries_file)
        
        print(f"Running analysis for library group: {args.library_group}")
        print(f"Found {len(libraries)} libraries to process")
        print("=" * 80)
        
        # Run analysis for each library separately
        results = []
        for idx, library_filename in enumerate(libraries, 1):
            library_name = library_filename.replace(".tsv", "")
            jsonl_filename = f"{library_name}.signatures.jsonl"
            jsonl_path = precomputed_dir / jsonl_filename
            
            print(f"\n[{idx}/{len(libraries)}] Processing library: {library_filename}")
            print(f"  Precomputed signatures: {jsonl_path}")
            
            if not jsonl_path.exists():
                print(f"  WARNING: Precomputed signatures not found, skipping.")
                continue
            
            # Create a unique entity_id for this library
            library_entity_id = f"{args.entity_id}_{library_name}"
            
            try:
                result = geneset_embedding(
                    embeddings_path=args.embeddings_path,
                    entity_id=library_entity_id,
                    output_folder=args.output_folder,
                    input_genes=input_genes,
                    input_genes_file=args.input_genes_file,
                    reference_genes=reference_genes,
                    reference_genes_file=args.reference_genes_file,
                    precomputed_terms_jsonl=str(jsonl_path),
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
                results.append({
                    "library": library_filename,
                    "success": result.get("success", False),
                    "entity_id": library_entity_id
                })
                print(f"  ✓ Completed: {library_filename}")
            except Exception as e:
                print(f"  ✗ ERROR processing {library_filename}: {str(e)}")
                results.append({
                    "library": library_filename,
                    "success": False,
                    "error": str(e)
                })
        
        # Print summary
        print("\n" + "=" * 80)
        print("ANALYSIS SUMMARY")
        print("=" * 80)
        successful = sum(1 for r in results if r.get("success", False))
        print(f"Total libraries:  {len(libraries)}")
        print(f"Successful:       {successful}")
        print(f"Failed:           {len(results) - successful}")
        print("=" * 80)
        
    else:
        # Single file mode (original behavior)
        if not args.precomputed_terms_jsonl:
            raise ValueError(
                "Either --library-group or --precomputed-terms-jsonl must be specified"
            )
        
        geneset_embedding(
            embeddings_path=args.embeddings_path,
            entity_id=args.entity_id,
            output_folder=args.output_folder,
            input_genes=input_genes,
            input_genes_file=args.input_genes_file,
            reference_genes=reference_genes,
            reference_genes_file=args.reference_genes_file,
            precomputed_terms_jsonl=args.precomputed_terms_jsonl,
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



