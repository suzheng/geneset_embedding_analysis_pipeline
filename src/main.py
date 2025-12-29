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
        group_name: Either a specific group name, "all", or "all-but-disease"
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
    elif group_name.lower() == "all-but-disease":
        # Return all libraries except from disease-related groups
        excluded_groups = {"Disease_and_Phenotype_Associations", "Rare_Disease_Associations"}
        all_libraries = []
        for group, libraries in group_to_libraries.items():
            if group not in excluded_groups:
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
    p.add_argument("--associated-item", required=True, help="What the input genes are associated with (e.g., disease, phenotype, molecular event)")
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
        help="Library group name from Group_Clean column in selected_libraries.txt, or 'all' for all groups, or 'all-but-disease' to exclude disease-related groups"
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
    p.add_argument(
        "--find-opposite-terms",
        default="auto",
        help="Find most opposite/different terms: true|false|auto (auto enables for Drug_and_Chemical_Perturbations, Gene_Perturbations, Protein_Domains_and_Complexes, Gene_Regulation, Miscellaneous groups)"
    )

    p.add_argument("--enable-clustering", action="store_true", default=True, help="Enable optional clustering")
    p.add_argument("--no-clustering", dest="enable_clustering", action="store_false", help="Disable clustering")
    p.add_argument("--cluster-min-genes", type=int, default=20, help="Min genes to attempt clustering")
    p.add_argument("--random-state", type=int, default=17, help="Random seed for deterministic behavior")

    # Figure generation options
    p.add_argument("--generate-figures", action="store_true", default=True, help="Generate publication-quality figures")
    p.add_argument("--figure-formats", nargs='+', default=["svg"], help="Figure formats (default: svg png)")
    
    # Output options
    p.add_argument("--show-top5-matches", action="store_true", default=True, help="Show top 5 matches with term names and similarity/distance scores in stdout for each library")

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
        print()
        print("Output files for each library (naming pattern: {entity_id}_{library_name}_geneset_embed_*):")
        print()
        print("  Main Result Files:")
        print("    - *_qc.csv                    : Quality control metrics (gene coverage, clustering)")
        print("    - *_cohort_signatures.csv     : Sample signature matrix (rows=samples, cols=dimensions)")
        print("    - *_signature.json            : Sample-specific signature (mean, std vectors)")
        print("    - *_term_matches.csv          : Top matching terms from the library")
        print("    - *_run_info.json             : Run metadata and parameters")
        print()
        print("  Documentation Files:")
        print("    - *_preprocess.PARAMS.md      : Preprocessing parameters")
        print("    - *_preprocess.RESULTS.md     : Preprocessing numerical results")
        print("    - *_qc.PARAMS.md              : QC parameters")
        print("    - *_qc.RESULTS.md             : QC numerical results")
        print("    - *_signature.PARAMS.md       : Signature computation parameters")
        print("    - *_signature.RESULTS.md      : Signature numerical results")
        print("    - *_term_matching.PARAMS.md   : Term matching parameters")
        print("    - *_term_matching.RESULTS.md  : Term matching numerical results")
        print("    - *_cohort_integration.PARAMS.md   : Cohort integration parameters")
        print("    - *_cohort_integration.RESULTS.md  : Cohort integration numerical results")
        print()
        print("  Shared Documentation (in shared_doc/ subfolder):")
        print("    - SHARED_*.METHOD.md          : Methodology descriptions (same across runs)")
        print("    - SHARED_*.RESULT_INTERPRETATION.md : Result interpretation guides")
        print()
        
        if args.generate_figures:
            print("  Publication-Quality Figures (in figures/ subfolder):")
            print("    Generated when --generate-figures is enabled")
            print()
            print("    1. term_signatures_distributions (.svg, .data.csv)")
            print("       Distribution of gene coverage across precomputed library terms")
            print("       Shows quality metrics: genes found per term and coverage fraction")
            print()
            print("    2. dimension_filtering (.svg, .data.csv)")
            print("       Standard deviation distribution across embedding dimensions")
            print("       Visualizes dimension selection threshold for filtering")
            print()
            print("    3. dimension_filtering_bar (.svg, .data.csv)")
            print("       Bar chart showing selected (orange) vs filtered (blue) dimensions")
            print("       Highlights which dimensions contribute to similarity matching")
            print()
            print("    4. top_matches_example (.svg, .data.csv)")
            print("       Horizontal bar chart of top matching terms from reference library")
            print("       Shows closest functional annotations with variance-weighted distances")
            print()
            print("    5. similarity_distributions (.svg, .data.csv)")
            print("       Dual histogram: distance and cosine similarity distributions")
            print("       Statistical overview of all term similarity scores")
            print()
            print("    6. distance_cosine_hexbin (.svg, .data.csv)")
            print("       2D hexagonal density plot of distance vs cosine similarity")
            print("       Shows correlation between the two similarity metrics")
            print()
            print("    7. pca_term_space (.svg, .data.csv)")
            print("       PCA projection of all term signatures in 2D space")
            print("       Top matches highlighted to show relative positioning")
            print()
            print("    8. signature_heatmap (.svg, .data.json)")
            print("       Side-by-side heatmaps of signature mean and std vectors")
            print("       Visual representation of embedding signature across dimensions")
            print()
            print("    Additional figures for opposite/different term matching:")
            print("    (Generated when --find-opposite-terms is enabled or set to 'auto' for applicable library groups)")
            print()
            print("    9. top_opposite_matches_example (.svg, .data.csv)")
            print("       Horizontal bar chart of top opposite/different terms from reference library")
            print("       Shows most dissimilar functional annotations based on negative cosine similarity")
            print()
            print("    10. distance_cosine_hexbin_opposite (.svg, .data.csv)")
            print("        2D hexagonal density plot with opposite matches highlighted")
            print("        Same base plot as #6 but highlights the most opposite terms")
            print()
            print("    11. pca_term_space_opposite (.svg, .data.csv)")
            print("        PCA projection with opposite matches highlighted")
            print("        Shows spatial relationship of opposite terms to input genes")
            print()
            print("    Each figure includes a companion data file (.csv or .json) for reproduction")
            print()
        
        print("=" * 80)
        
        # Create parent folder for this entity-id if output_folder is specified
        parent_output_folder = args.output_folder
        if args.output_folder:
            parent_output_folder = os.path.join(args.output_folder, args.entity_id)
        else:
            parent_output_folder = args.entity_id
        
        # Run analysis for each library separately
        results = []
        for idx, library_filename in enumerate(libraries, 1):
            library_name = library_filename.replace(".tsv", "")
            jsonl_filename = f"{library_name}.signatures.jsonl"
            jsonl_path = precomputed_dir / jsonl_filename
            
            print(f"\n[{idx}/{len(libraries)}] Processing library: {library_filename}")
            
            if not jsonl_path.exists():
                print(f"  WARNING: Precomputed signatures not found, skipping.")
                continue
            
            # Create library-specific entity_id and output folder
            library_entity_id = f"{args.entity_id}_{library_name}"
            # Create a subfolder for each library under the parent folder
            library_output_folder = os.path.join(parent_output_folder, library_entity_id)
            
            try:
                result = geneset_embedding(
                    embeddings_path=args.embeddings_path,
                    entity_id=library_entity_id,
                    associated_item=args.associated_item,
                    output_folder=library_output_folder,  # Library-specific subfolder
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
                    find_opposite_terms=args.find_opposite_terms,
                    enable_clustering=args.enable_clustering,
                    cluster_min_genes=args.cluster_min_genes,
                    random_state=args.random_state,
                    verbose_output=False,  # Suppress detailed OUTPUT: lines for batch processing
                    generate_figures=args.generate_figures,
                    figure_formats=args.figure_formats,
                    show_top5_matches=args.show_top5_matches,
                    shared_doc_dir=parent_output_folder,  # Shared docs go in parent folder
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
        # Print output location
        if args.output_folder:
            base_path = Path(args.output_folder) / args.entity_id
        else:
            base_path = Path(args.entity_id)
        print(f"\nOutput structure:")
        print(f"  {base_path}/")
        print(f"  ├── shared_doc/  (shared documentation for all libraries)")
        print(f"  │   ├── BACKGROUND.md")
        print(f"  │   ├── FIGURE_LEGEND.md")
        print(f"  │   └── SHARED_*.md files")
        print(f"  ├── {args.entity_id}_<library_1>/  (library 1 results)")
        print(f"  │   ├── figures/")
        print(f"  │   ├── *.csv, *.json, *.md files")
        print(f"  │   └── ...")
        print(f"  ├── {args.entity_id}_<library_2>/  (library 2 results)")
        print(f"  │   ├── figures/")
        print(f"  │   └── ...")
        print(f"  └── ...")
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
            associated_item=args.associated_item,
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
            find_opposite_terms=args.find_opposite_terms,
            enable_clustering=args.enable_clustering,
            cluster_min_genes=args.cluster_min_genes,
            random_state=args.random_state,
            generate_figures=args.generate_figures,
            figure_formats=args.figure_formats,
            show_top5_matches=args.show_top5_matches,
        )


if __name__ == "__main__":
    main()



