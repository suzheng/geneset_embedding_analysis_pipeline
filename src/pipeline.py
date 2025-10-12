from __future__ import annotations

import csv
import json
import os
import re
import sys
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from .cluster import maybe_cluster_vectors
from .embeddings_io import get_vectors_for_genes, load_embeddings
from .match import TermSignature, rank_terms, rank_terms_farthest, rank_terms_most_opposite
from .report import ensure_dir, write_json
from .signature import (
    Signature,
    compute_signature,
    compute_differential_signature,
    select_dims_by_cv,
    select_dims_by_quantile,
    select_dims_by_zscore,
    select_top_k_low_variance,
)


# ---------------------- Existing low-level API ----------------------

def compute_input_signatures(
    embeddings_path: str,
    input_genes: Iterable[str],
    enable_clustering: bool = True,
    cluster_min_genes: int = 20,
    reference_genes: Optional[Iterable[str]] = None,
) -> Tuple[List[Signature], List[str], List[str], int]:
    """
    Compute signatures for input genes, optionally as differential vs reference genes.
    
    Args:
        embeddings_path: Path to embeddings file
        input_genes: Input gene symbols
        enable_clustering: Whether to cluster (ignored if reference_genes provided)
        cluster_min_genes: Minimum genes for clustering
        reference_genes: Optional reference gene symbols for differential mode
        
    Returns:
        Tuple of (signatures, found_genes, missing_genes, embedding_dimension)
    """
    symbol_to_vector, _n, d = load_embeddings(embeddings_path)
    vectors, found, missing = get_vectors_for_genes(symbol_to_vector, input_genes)
    if vectors.shape[0] == 0:
        raise ValueError("No input genes found in embeddings.")

    # Differential mode: compute (input - reference)
    if reference_genes is not None:
        ref_vectors, ref_found, ref_missing = get_vectors_for_genes(symbol_to_vector, reference_genes)
        if ref_vectors.shape[0] == 0:
            raise ValueError("No reference genes found in embeddings.")
        
        # Compute differential signature (skip clustering in this mode)
        diff_sig = compute_differential_signature(vectors, ref_vectors, n_requested=len(input_genes))
        signatures: List[Signature] = [diff_sig]
        print(f"INFO: Differential mode - input genes: {len(found)}, reference genes: {len(ref_found)}, missing: {len(missing)} input + {len(ref_missing)} reference")
        return signatures, found, missing, d
    
    # Standard mode: cluster and compute signatures
    groups = maybe_cluster_vectors(
        vectors, enable_clustering, min_genes_for_clustering=cluster_min_genes
    )
    signatures = [compute_signature(g, n_requested=len(input_genes)) for g in groups]
    return signatures, found, missing, d


def load_precomputed_terms_jsonl(path: str) -> List[TermSignature]:
    results: List[TermSignature] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            results.append(
                TermSignature(
                    term_id=row["term_id"],
                    library=row["library"],
                    term=row["term"],
                    mean=np.array(row["mean"], dtype=np.float32),
                    std=np.array(row["std"], dtype=np.float32),
                    n_found=int(row["n_found"]),
                    coverage=float(row["coverage"]),
                )
            )
    return results


def make_dims_mask(
    method: Optional[str],
    mean: np.ndarray,
    std: np.ndarray,
    *,
    q_low: float = 0.3,
    cv_threshold: float = 1.0,
    z_threshold: float = 0.5,
    top_k: int = 8,
) -> Optional[np.ndarray]:
    if method is None or method.lower() == "none":
        return None
    name = method.lower()
    if name == "quantile":
        return select_dims_by_quantile(std, q_low=q_low)
    if name == "cv":
        return select_dims_by_cv(mean, std, cv_threshold=cv_threshold)
    if name == "zscore":
        return select_dims_by_zscore(std, z_threshold=z_threshold)
    if name == "topk":
        return select_top_k_low_variance(std, k=top_k)
    raise ValueError(f"Unknown dims filter method: {method}")


def match_input_to_terms(
    input_signature: Signature,
    term_signatures: List[TermSignature],
    dims_filter: Optional[str] = "quantile",
    **filter_kwargs,
):
    mask = make_dims_mask(dims_filter, input_signature.mean, input_signature.std, **filter_kwargs)
    return rank_terms(input_signature.mean, input_signature.std, term_signatures, dims_mask=mask)


def match_input_to_terms_farthest(
    input_signature: Signature,
    term_signatures: List[TermSignature],
    dims_filter: Optional[str] = "quantile",
    **filter_kwargs,
):
    mask = make_dims_mask(dims_filter, input_signature.mean, input_signature.std, **filter_kwargs)
    return rank_terms_farthest(
        input_signature.mean, input_signature.std, term_signatures, dims_mask=mask
    )


def match_input_to_terms_most_opposite(
    input_signature: Signature,
    term_signatures: List[TermSignature],
    dims_filter: Optional[str] = "quantile",
    **filter_kwargs,
):
    mask = make_dims_mask(dims_filter, input_signature.mean, input_signature.std, **filter_kwargs)
    return rank_terms_most_opposite(
        input_signature.mean, input_signature.std, term_signatures, dims_mask=mask
    )


# ---------------------- High-level pipeline API ----------------------

def sanitize_simple_id(entity_id: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_-]", "", str(entity_id))


def should_find_opposite_terms(entity_id: str, find_opposite_terms: str) -> bool:
    """
    Determine if opposite terms should be found based on the find_opposite_terms setting
    and the library group detected from entity_id.
    
    Args:
        entity_id: The entity ID which may contain library name
        find_opposite_terms: "true"/"false"/"auto" (case-insensitive)
        
    Returns:
        True if opposite terms should be found, False otherwise
    """
    setting = str(find_opposite_terms).lower().strip()
    
    if setting in ["true", "t", "1", "yes"]:
        return True
    elif setting in ["false", "f", "0", "no"]:
        return False
    elif setting == "auto":
        # Groups that should automatically enable opposite term finding
        auto_enable_groups = [
            "Drug_and_Chemical_Perturbations",
            "Gene_Perturbations", 
            "Protein_Domains_and_Complexes",
            "Gene_Regulation",
            "Miscellaneous"
        ]
        
        # Also check library filenames that belong to these groups
        auto_enable_libraries = [
            # Drug & Chemical Perturbations
            "LINCS_L1000_Chem_Pert_Consensus_Sigs",
            "RummaGEO_DrugPerturbations_2025",
            "DrugMatrix",
            "DSigDB",
            "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
            # Gene Perturbations
            "LINCS_L1000_CRISPR_KO_Consensus_Sigs",
            "RummaGEO_GenePerturbations_2025",
            "PerturbAtlas",
            # Protein Domains & Complexes
            "InterPro_Domains_2019",
            "Pfam_Domains_2019",
            "CORUM",
            # Gene Regulation
            "Rummagene_transcription_factors",
            "TargetScan_microRNA",
            "lncHUB_lncRNA_Co-Expression",
            "TF_Perturbations_Followed_by_Expression",
            # Miscellaneous
            "Epigenomics_Roadmap_HM_ChIP-seq",
            "Virus-Host_PPI_P-HIPSTer_2020",
            "HMDB_Metabolites",
            "ProteomicsDB_2020",
            "GlyGen_Glycosylated_Proteins_2022",
            "Rummagene_kinases",
            "Chromosome_Location",
        ]
        
        entity_lower = entity_id.lower()
        
        # Check if entity_id contains any of the auto-enable groups or libraries
        for group in auto_enable_groups:
            if group.lower() in entity_lower:
                return True
        
        for lib in auto_enable_libraries:
            if lib.lower() in entity_lower:
                return True
        
        return False
    else:
        # Default to False if unrecognized setting
        return False


def generate_output_folder(entity_id: str, pipeline_short_name: str) -> str:
    return f"{entity_id}_{pipeline_short_name}"


def _get_workspace() -> str:
    workspace = os.getenv("PROJECT_WORKSPACE_VOLUME_PATH")
    if not workspace:
        # Fallback to current working directory with a warning for reproducibility.
        print(
            # "WARN: PROJECT_WORKSPACE_VOLUME_PATH is not set; using current working directory",
            file=sys.stderr,
        )
        workspace = os.getcwd()
    return os.path.abspath(workspace)


def get_output_path(output_folder: str) -> str:
    workspace = _get_workspace()
    return os.path.join(workspace, output_folder)


def _gather_versions() -> Dict[str, str]:
    versions: Dict[str, str] = {}
    versions["python"] = sys.version.split()[0]
    try:
        import numpy as _np  # noqa: F401

        versions["numpy"] = _np.__version__
    except Exception:
        pass
    try:
        import sklearn as _sk  # noqa: F401

        versions["scikit_learn"] = _sk.__version__
    except Exception:
        pass
    try:
        import hdbscan as _hb  # noqa: F401

        versions["hdbscan"] = getattr(_hb, "__version__", "unknown")
    except Exception:
        pass
    return versions


def _write_text(path: str, text: str) -> None:
    ensure_dir(path)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def _write_method_doc(output_dir: str, entity_id: str, pipeline_short: str, step: str, text: str, verbose: bool = True) -> str:
    path = os.path.join(output_dir, f"{entity_id}_{pipeline_short}_{step}.METHOD.md")
    _write_text(path, text)
    _print_output(path, "md", f"{step} methodology", verbose=verbose)
    return path


def _write_shared_method_doc(entity_id: str, output_dir: str, step: str, text: str, verbose: bool = True) -> str:
    """Write shared methodology that is the same across all runs."""
    shared_dir = os.path.join(output_dir, "shared_doc")
    os.makedirs(shared_dir, exist_ok=True)
    path = os.path.join(shared_dir, f"SHARED_{step}.METHOD.md")
    _write_text(path, text)
    _print_output(path, "md", f"Shared {step} methodology", verbose=verbose)
    return path


def _write_run_params_doc(output_dir: str, entity_id: str, pipeline_short: str, step: str, text: str, verbose: bool = True) -> str:
    """Write run-specific parameters for a methodology step."""
    path = os.path.join(output_dir, f"{entity_id}_{pipeline_short}_{step}.PARAMS.md")
    _write_text(path, text)
    _print_output(path, "md", f"{step} run parameters", verbose=verbose)
    return path


def _write_result_doc(output_dir: str, entity_id: str, pipeline_short: str, step: str, text: str, verbose: bool = True) -> str:
    path = os.path.join(output_dir, f"{entity_id}_{pipeline_short}_{step}.RESULT.md")
    _write_text(path, text)
    _print_output(path, "md", f"{step} results", verbose=verbose)
    return path


def _write_shared_result_doc(entity_id: str, output_dir: str, step: str, text: str, verbose: bool = True) -> str:
    """Write shared result interpretation that is the same across all runs."""
    shared_dir = os.path.join(output_dir, "shared_doc")
    os.makedirs(shared_dir, exist_ok=True)
    path = os.path.join(shared_dir, f"SHARED_{step}.RESULT_INTERPRETATION.md")
    _write_text(path, text)
    _print_output(path, "md", f"Shared {step} result interpretation", verbose=verbose)
    return path


def _copy_background_doc(output_dir: str, verbose: bool = True) -> str:
    """Copy BACKGROUND.md from src directory to shared_doc folder."""
    shared_dir = os.path.join(output_dir, "shared_doc")
    os.makedirs(shared_dir, exist_ok=True)
    
    # Determine the source path for BACKGROUND.md in the src directory
    # Get the directory where this pipeline.py file is located
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    source_path = os.path.join(pipeline_dir, "BACKGROUND.md")
    dest_path = os.path.join(shared_dir, "BACKGROUND.md")
    
    # Copy the file if it exists
    if os.path.exists(source_path):
        import shutil
        shutil.copy2(source_path, dest_path)
        _print_output(dest_path, "md", "Pipeline background documentation", verbose=verbose)
        return dest_path
    else:
        if verbose:
            print(f"WARNING: BACKGROUND.md not found at {source_path}", file=sys.stderr)
        return ""


def _write_figure_legend_doc(entity_id: str, output_dir: str, library_name: str, find_opposite_terms: str, verbose: bool = True) -> str:
    """Write comprehensive figure legends for all generated figures."""
    shared_dir = os.path.join(output_dir, "shared_doc")
    os.makedirs(shared_dir, exist_ok=True)
    path = os.path.join(shared_dir, "FIGURE_LEGEND.md")
    
    # Determine if opposite figures are included
    has_opposite = should_find_opposite_terms(entity_id, find_opposite_terms)
    
    # Generate comprehensive figure legends following academic publication standards
    legend_text = f"""# Figure Legends

This document provides comprehensive figure titles and legends for all generated figures, following academic publication standards.

## Figure 1: Gene Set Coverage Distribution in Reference Library

**Title:** Distribution of Gene Coverage Across Precomputed Term Signatures in {library_name or 'Reference Library'}

**Legend:** 
(A) Histogram showing the distribution of the number of genes found per term signature in the reference library. Each bar represents the frequency of terms with a given number of genes successfully mapped to the embedding space. This distribution reflects the size variability of gene sets in the reference library and their representation in the pre-trained embedding model.
(B) Histogram depicting the coverage distribution (fraction of genes with available embeddings) across all term signatures. Coverage values range from 0 to 1, where higher values indicate better representation of gene sets in the embedding space. High coverage (>0.8) suggests comprehensive gene annotation, while lower coverage may indicate specialized or poorly characterized gene sets.

**File:** `term_signatures_distributions.svg`

---

## Figure 2: Embedding Dimension Standard Deviation Profile

**Title:** Standard Deviation Profile Across Embedding Dimensions for Query Gene Set

**Legend:** 
Line plot showing the standard deviation (σ) of the query gene set signature across all embedding dimensions, sorted in ascending order. The horizontal dashed line indicates the quantile threshold (q = 0.3) used for dimension selection. Dimensions below this threshold exhibit low within-set variability and are retained for downstream matching, as they represent stable, coherent semantic features. Dimensions above the threshold show high variability and are filtered out to focus the analysis on the most informative features. This dimension selection strategy enhances the signal-to-noise ratio in similarity calculations.

**File:** `dimension_filtering.svg`

---

## Figure 3: Dimension Selection Summary

**Title:** Selected and Filtered Embedding Dimensions Based on Variability Threshold

**Legend:** 
Bar chart visualizing the standard deviation of each embedding dimension for the query gene set. Dimensions are colored by their selection status: orange bars represent selected dimensions (low variability, below the 0.3 quantile threshold) that are retained for matching, while blue bars indicate filtered dimensions (high variability) that are excluded from similarity calculations. This visualization provides a comprehensive view of which semantic features contribute to the functional interpretation of the gene set.

**File:** `dimension_filtering_bar.svg`

---

## Figure 4: Top Functionally Similar Terms

**Title:** Top Matching Gene Set Terms Ranked by Variance-Weighted Distance

**Legend:** 
Horizontal bar chart displaying the top 10 most similar terms from the reference library ({library_name or 'reference library'}), ranked by variance-weighted distance (lower values indicate stronger similarity). Each bar represents the distance metric, which combines the Euclidean distance between signature mean vectors with inverse-variance weighting to prioritize stable dimensions. Term names are displayed on the y-axis, ordered from most similar (top) to least similar (bottom). This ranking identifies the closest functional relatives of the query gene set in the reference library, providing biological context and potential functional annotations.

**File:** `top_matches_example.svg`

---

## Figure 5: Similarity Metric Distributions

**Title:** Distribution of Distance and Cosine Similarity Scores Across All Reference Terms

**Legend:** 
(A) Distribution of variance-weighted distance scores computed between the query gene set signature and all reference term signatures. Lower distance values indicate greater similarity. The distribution typically shows a right-skewed pattern, with most terms at moderate distances and a tail of highly similar matches. Kernel density estimation (KDE) overlay highlights the probability density.
(B) Distribution of cosine similarity scores between the query signature and all reference terms. Cosine similarity ranges from -1 (opposite direction) to +1 (same direction), measuring directional alignment in embedding space. Values near +1 indicate strong functional similarity, while values near 0 suggest orthogonal (unrelated) gene sets. These distributions provide a global view of how the query gene set relates to the entire reference library.

**File:** `similarity_distributions.svg`

---

## Figure 6: Joint Distance-Similarity Landscape

**Title:** Two-Dimensional Density Plot of Distance vs. Cosine Similarity for Reference Terms

**Legend:** 
Hexagonal binning plot showing the joint distribution of variance-weighted distance (x-axis) and cosine similarity (y-axis) across all reference term signatures. The color intensity (blue gradient) represents the density of terms in each hexagonal bin, revealing the relationship between these two complementary similarity metrics. Overlaid magenta points highlight the top matching terms, which typically cluster in the low-distance, high-cosine-similarity region (bottom-right quadrant). This visualization helps identify whether top matches are driven primarily by distance minimization, directional alignment, or both.

**File:** `distance_cosine_hexbin.svg`

---

## Figure 7: PCA Projection of Reference Term Space

**Title:** Principal Component Analysis of Reference Term Signatures with Top Matches Highlighted

**Legend:** 
Two-dimensional PCA projection of all reference term signatures (mean embedding vectors), with PC1 and PC2 axes labeled by their explained variance ratios. Background points (light blue) represent all terms in the reference library, illustrating the global structure of functional gene set space. Green points indicate the top matching terms, showing their spatial proximity to the query gene set. The query gene set signature is shown as an orange star, positioned within the reference term landscape. Terms that cluster near the query are functionally related, while distant terms represent distinct biological processes. This projection provides an intuitive geometric view of functional similarity relationships.

**File:** `pca_term_space.svg`

---

## Figure 8: Gene Set Signature Heatmap

**Title:** Heatmap Visualization of Query Gene Set Signature Mean and Standard Deviation Vectors

**Legend:** 
(A) Heatmap displaying the mean embedding vector (μ) of the query gene set across all dimensions. Color intensity (coolwarm colormap) represents the magnitude and sign of each dimension's mean value, with blue indicating negative values and red indicating positive values. This vector captures the centroid or "average position" of the gene set in embedding space, encoding its dominant functional characteristics.
(B) Heatmap showing the standard deviation vector (σ) across all dimensions, visualized using a magma colormap (dark purple to bright yellow). Higher values (yellow) indicate dimensions with high within-set gene variability, suggesting functional heterogeneity, while lower values (purple) represent coherent dimensions where genes are tightly clustered. Together, these two vectors define the complete signature used for functional term matching.

**File:** `signature_heatmap.svg`

---
"""

    # Add opposite term figures if enabled
    if has_opposite:
        legend_text += f"""## Figure 9: Top Functionally Opposite Terms

**Title:** Top Opposite Gene Set Terms Ranked by Negative Cosine Similarity

**Legend:** 
Horizontal bar chart displaying the top 10 most functionally opposite terms from the reference library ({library_name or 'reference library'}), ranked by cosine similarity (most negative values indicate strongest opposition). Each bar represents the cosine similarity score, where values approaching -1 indicate gene sets with opposite directionality in embedding space. This analysis identifies biological processes or gene sets that are functionally antagonistic or contrasting to the query gene set, which is particularly valuable for perturbation studies (e.g., drug vs. disease, knockdown vs. overexpression) where identifying opposing biological states is informative.

**File:** `top_opposite_matches_example.svg`

---

## Figure 10: Joint Distance-Similarity Landscape with Opposite Matches

**Title:** Two-Dimensional Density Plot Highlighting Top Opposite Gene Set Terms

**Legend:** 
Hexagonal binning plot showing the joint distribution of variance-weighted distance (x-axis) and cosine similarity (y-axis) across all reference term signatures. The color intensity (blue gradient) represents term density. Overlaid orange points highlight the top opposite terms, which typically cluster in regions with negative cosine similarity values (upper portion of plot), indicating directional opposition to the query gene set. This visualization reveals whether opposite matches are characterized by large distances, negative directional alignment, or both, providing insight into the nature of functional opposition.

**File:** `distance_cosine_hexbin_opposite.svg`

---

## Figure 11: PCA Projection with Opposite Matches Highlighted

**Title:** Principal Component Analysis of Reference Terms with Functionally Opposite Matches

**Legend:** 
Two-dimensional PCA projection of all reference term signatures, with PC1 and PC2 axes labeled by their explained variance ratios. Background points (light blue) represent all terms in the reference library. Orange points indicate the top opposite matching terms, showing their spatial positioning relative to the query gene set (orange star). Terms positioned on the opposite side of the PCA space from the query typically represent contrasting biological processes or phenotypes. This geometric visualization helps interpret the nature of functional opposition in terms of the major axes of variation in the reference library.

**File:** `pca_term_space_opposite.svg`

---
"""

    # Add footer with general notes
    legend_text += """## General Notes

**Figure Format:** All figures are provided in multiple formats (SVG for publications, PNG for presentations) to ensure broad compatibility.

**Data Files:** Each figure is accompanied by a corresponding data file (`.data.csv` or `.data.json`) containing the underlying numerical data used to generate the visualization, ensuring full reproducibility and enabling custom reanalysis.

**Color Schemes:** Figures use a consistent, colorblind-friendly palette optimized for scientific publication. Background grid lines are removed per publication standards, and font sizes are optimized for readability in both print and digital formats.

**Statistical Methods:** All distance calculations use variance-weighted metrics that prioritize stable dimensions (low within-set variability) and down-weight noisy dimensions (high variability). Cosine similarity provides a complementary directional measure. Dimension selection (quantile filtering at q=0.3) focuses the analysis on the most informative semantic features.

**Reproducibility:** All figures are generated using standardized matplotlib configuration with embedded fonts (Type 42 PDF fonts, SVG none), ensuring consistent appearance across platforms and publications. Random seeds are fixed for all analyses involving stochastic components (e.g., PCA initialization).

---

**Generated by:** Gene Set Embedding Analysis Pipeline  
**Document Version:** 1.0  
**Last Updated:** {datetime.utcnow().strftime('%Y-%m-%d')}
"""

    _write_text(path, legend_text)
    _print_output(path, "md", "Figure legends and descriptions", verbose=verbose)
    return path


def _write_run_results_doc(output_dir: str, entity_id: str, pipeline_short: str, step: str, text: str, verbose: bool = True) -> str:
    """Write run-specific numerical results."""
    path = os.path.join(output_dir, f"{entity_id}_{pipeline_short}_{step}.RESULTS.md")
    _write_text(path, text)
    _print_output(path, "md", f"{step} numerical results", verbose=verbose)
    return path


def _read_gene_list_from_file(path: str) -> List[str]:
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Support comma or whitespace separated tokens per line
            parts = re.split(r"[\s,]+", line)
            for p in parts:
                token = (p or "").strip()
                if token:
                    genes.append(token)
    return genes


def _write_csv(path: str, header: List[str], rows: List[List[Any]]) -> None:
    ensure_dir(path)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for r in rows:
            writer.writerow(r)


def _print_output(path: str, ftype: str, desc: str, verbose: bool = True) -> None:
    """Print output file information.
    
    Args:
        path: File path
        ftype: File type (csv, json, md)
        desc: Description
        verbose: If False, suppress output (default True for backward compatibility)
    """
    if verbose:
        abs_path = os.path.abspath(path)
        print(f"OUTPUT: file={abs_path} type={ftype} desc=\"{desc}\"")


def _generate_figures(
    figures_dir: str,
    sig: Signature,
    term_signatures: Optional[List[TermSignature]],
    formats: List[str],
    entity_id: str,
    library_name: str = "",
    verbose: bool = True,
    find_opposite_terms: str = "auto"
) -> List[str]:
    """Generate publication-quality figures for the analysis.
    
    Returns list of paths to generated figure files.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        from scripts.generate_publication_figures import (
            _configure_matplotlib,
            _figure_term_signature_distributions,
            _figure_dimension_filtering,
            _figure_dimension_filtering_bar,
            _figure_top_matches,
            _figure_similarity_distributions,
            _figure_distance_cosine_hexbin,
            _figure_pca_term_space,
            _figure_signature_heatmap,
            _figure_top_matches_opposite,
            _figure_distance_cosine_hexbin_opposite,
            _figure_pca_term_space_opposite,
        )
        
        _configure_matplotlib()
        
        figure_paths: List[str] = []
        
        # Generate figures without progress messages (they're described in initial stdout)
        if term_signatures and len(term_signatures) > 0:
            fig_path, data_path = _figure_term_signature_distributions(figures_dir, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
        
        fig_path, data_path = _figure_dimension_filtering(figures_dir, sig, formats, library_name)
        figure_paths.extend([fig_path, data_path])
        
        fig_path, data_path = _figure_dimension_filtering_bar(figures_dir, sig, formats, library_name)
        figure_paths.extend([fig_path, data_path])
        
        if term_signatures and len(term_signatures) > 0:
            fig_path, data_path = _figure_top_matches(figures_dir, sig, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
            
            fig_path, data_path = _figure_similarity_distributions(figures_dir, sig, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
            
            fig_path, data_path = _figure_distance_cosine_hexbin(figures_dir, sig, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
            
            fig_path, data_path = _figure_pca_term_space(figures_dir, sig, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
        
        fig_path, data_path = _figure_signature_heatmap(figures_dir, sig, formats, library_name)
        figure_paths.extend([fig_path, data_path])
        
        # Generate opposite term figures if enabled
        if should_find_opposite_terms(entity_id, find_opposite_terms) and term_signatures and len(term_signatures) > 0:
            fig_path, data_path = _figure_top_matches_opposite(figures_dir, sig, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
            
            fig_path, data_path = _figure_distance_cosine_hexbin_opposite(figures_dir, sig, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
            
            fig_path, data_path = _figure_pca_term_space_opposite(figures_dir, sig, term_signatures, formats, library_name)
            figure_paths.extend([fig_path, data_path])
        
        # Create index file
        index_path = os.path.join(figures_dir, f"{entity_id}_figure_index.json")
        import json
        
        descriptions = {
            "term_signatures_distributions": "Distribution of gene coverage across precomputed term signatures",
            "dimension_filtering": "Standard deviation distribution across embedding dimensions",
            "dimension_filtering_bar": "Bar chart showing selected vs filtered dimensions",
            "top_matches_example": "Top matching terms from the reference library",
            "similarity_distributions": "Distribution of distance and cosine similarity scores",
            "distance_cosine_hexbin": "2D density plot of distance vs cosine similarity",
            "pca_term_space": "PCA projection of term signatures with top matches highlighted",
            "signature_heatmap": "Heatmap of signature mean and std vectors",
        }
        
        if should_find_opposite_terms(entity_id, find_opposite_terms):
            descriptions.update({
                "top_opposite_matches_example": "Top opposite/different terms from the reference library",
                "distance_cosine_hexbin_opposite": "2D density plot with opposite matches highlighted",
                "pca_term_space_opposite": "PCA projection with opposite matches highlighted",
            })
        
        with open(index_path, "w", encoding="utf-8") as f:
            json.dump({
                "entity_id": entity_id,
                "figures_directory": os.path.abspath(figures_dir),
                "figure_files": [os.path.abspath(p) for p in figure_paths],
                "formats": formats,
                "description": descriptions
            }, f, indent=2)
        figure_paths.append(index_path)
        
        # Brief summary without detailed descriptions (already shown in initial output)
        if verbose:
            print(f"  Generated {len(figure_paths)} figure and data files")
        
        return figure_paths
        
    except ImportError as e:
        print(f"WARNING: Could not generate figures - missing dependency: {e}")
        return []
    except Exception as e:
        print(f"WARNING: Figure generation failed: {e}")
        return []


def geneset_embedding(
    embeddings_path: str,
    entity_id: str,
    output_folder: Optional[str] = None,
    input_genes: Optional[Iterable[str]] = None,
    input_genes_file: Optional[str] = None,
    reference_genes: Optional[Iterable[str]] = None,
    reference_genes_file: Optional[str] = None,
    precomputed_terms_jsonl: Optional[str] = None,
    steps_to_run: Optional[str] = None,
    # Matching and signature options
    dims_filter: Optional[str] = "quantile",
    q_low: float = 0.3,
    cv_threshold: float = 1.0,
    z_threshold: float = 0.5,
    top_k: int = 25,
    find_opposite_terms: str = "auto",
    enable_clustering: bool = True,
    cluster_min_genes: int = 20,
    random_state: int = 17,
    verbose_output: bool = True,
    generate_figures: bool = False,
    figure_formats: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Production-oriented pipeline for gene set embedding matching with flexible inputs.

    Returns a compact dict including success flag, run_info, and artifact paths.
    """
    np.random.seed(random_state)

    clean_id = sanitize_simple_id(entity_id)
    if not clean_id:
        raise ValueError("ERROR: Invalid entity_id; field=entity_id; expected=alphanumeric/_/-")

    pipeline_short = "geneset_embed"
    # If output_folder not specified, use current_working_directory/{entity_id}
    folder_name = output_folder or clean_id
    output_dir = get_output_path(folder_name)
    os.makedirs(output_dir, exist_ok=True)

    if verbose_output:
        print(f"INFO: embeddings_path={os.path.abspath(embeddings_path)} entity_id={clean_id}")

    # Determine steps
    all_steps = [
        "preprocess",
        "qc",
        "signature",
        "term_matching",
        "cohort_integration",
        "group_comparison",
    ]
    steps_set = set([s.strip().lower() for s in (steps_to_run.split(",") if steps_to_run else all_steps) if s.strip()])

    # Load embeddings once
    symbol_to_vector, _n, d = load_embeddings(embeddings_path)

    # Prepare inputs: single sample mode
    samples: List[Dict[str, Any]] = []
    artifacts: List[str] = []
    versions = _gather_versions()

    # Single-sample mode
    if input_genes is None and not input_genes_file:
        raise ValueError(
            "ERROR: no input genes provided; field=input_genes|input_genes_file expected non-empty"
        )
    if input_genes is None:
        input_genes = _read_gene_list_from_file(str(input_genes_file))
    
    # Read reference genes if provided
    ref_genes_list: Optional[List[str]] = None
    if reference_genes is not None:
        ref_genes_list = list(reference_genes)
    elif reference_genes_file is not None:
        ref_genes_list = _read_gene_list_from_file(str(reference_genes_file))
    
    samples.append({"sample_id": clean_id, "genes": list(input_genes), "reference_genes": ref_genes_list})

    # Copy background documentation to shared_doc folder
    background_path = _copy_background_doc(output_dir, verbose=verbose_output)
    if background_path:
        artifacts.append(os.path.abspath(background_path))

    # Documentation: preprocess step
    if "preprocess" in steps_set:
        is_differential = ref_genes_list is not None and len(ref_genes_list) > 0
        
        # Shared methodology (same across all runs)
        shared_method_text = (
            "## Data Preprocessing\n\n"
            "### Gene Set Input and Normalization\n\n"
            "Input gene lists are provided as collections of human gene symbols. "
            "All gene symbols are normalized through a standardized procedure: "
            "whitespace trimming, uppercase conversion (e.g., 'brca1' → 'BRCA1'), "
            "and removal of empty entries. Malformed entries are filtered during embedding lookup.\n\n"
            "### Differential Embedding Mode (Optional)\n\n"
            "When a reference gene set is provided, the analysis computes differential embeddings. "
            "The differential signature is computed as the difference between the mean input gene set embedding "
            "and the mean reference gene set embedding: μ_diff = μ_input - μ_reference. "
            "This approach highlights the distinctive features of the input gene set relative to the reference. "
            "Clustering is automatically disabled in differential mode to maintain a single comparative signature.\n\n"
        )
        _write_shared_method_doc(clean_id, output_dir, "preprocess", shared_method_text, verbose=verbose_output)
        
        # Run-specific parameters
        hdbscan_text = ""
        if versions.get('hdbscan'):
            hdbscan_ver = versions.get('hdbscan', 'N/A')
            if hdbscan_ver and hdbscan_ver != 'unknown':
                hdbscan_text = f"HDBSCAN {hdbscan_ver} was used for density-based clustering. "
            else:
                hdbscan_text = "HDBSCAN was used for density-based clustering. "
        
        params_text = (
            f"## Run-Specific Parameters\n\n"
            f"### Software Environment\n\n"
            f"- **Python version**: {versions.get('python', 'N/A')}\n"
            f"- **NumPy version**: {versions.get('numpy', 'N/A')}\n"
            f"- **scikit-learn version**: {versions.get('scikit_learn', 'N/A')}\n"
        )
        if hdbscan_text:
            params_text += f"- **HDBSCAN version**: {versions.get('hdbscan', 'N/A')}\n"
        params_text += (
            f"\n### Analysis Parameters\n\n"
            f"- **Random seed**: {random_state}\n"
            f"- **Differential mode**: {'Enabled' if is_differential else 'Disabled'}\n"
        )
        if is_differential:
            params_text += f"- **Reference genes provided**: {len(ref_genes_list)}\n"
        params_text += "\n"
        
        _write_run_params_doc(output_dir, clean_id, pipeline_short, "preprocess", params_text, verbose=verbose_output)

        n_genes = len(samples[0].get('genes', [])) if samples else 0
        
        # Run-specific preprocessing results (minimal, no shared interpretation needed)
        result_text = (
            f"## Data Preprocessing Results\n\n"
            f"- **Input genes**: {n_genes}\n\n"
        )
        _write_run_results_doc(output_dir, clean_id, pipeline_short, "preprocess", result_text, verbose=verbose_output)

    # Per-sample processing
    per_sample_qc_rows: List[List[Any]] = []
    cohort_signature_matrix: List[List[Any]] = []
    cohort_header: List[str] = ["sample_id"] + [f"dim_{i}" for i in range(d)]
    term_signatures: Optional[List[TermSignature]] = None
    if precomputed_terms_jsonl and os.path.exists(precomputed_terms_jsonl):
        term_signatures = load_precomputed_terms_jsonl(precomputed_terms_jsonl)

    signature_files: List[str] = []
    term_matching_files: List[str] = []
    sample_match_summaries: List[str] = []
    total_files_written = 0  # Track all files written for accurate count

    for s in samples:
        sid = s["sample_id"]
        genes: List[str] = [g for g in s.get("genes", []) if g]
        if "preprocess" in steps_set:
            # nothing to do beyond parsing for this simple pipeline
            pass

        if not genes:
            print(f"WARN: sample {sid} has no genes after preprocessing; skipping", file=sys.stderr)
            continue

        # Compute signatures and QC
        if "signature" in steps_set or "qc" in steps_set:
            vectors, found, missing = get_vectors_for_genes(symbol_to_vector, genes)
            if vectors.shape[0] == 0:
                print(
                    f"WARN: sample {sid} contains no genes present in embeddings; skipping",
                    file=sys.stderr,
                )
                continue
            
            # Check if differential mode is enabled
            ref_genes = s.get("reference_genes")
            use_differential_mode = ref_genes is not None and len(ref_genes) > 0
            
            if use_differential_mode:
                # Differential mode: compute (input - reference)
                ref_vectors, ref_found, ref_missing = get_vectors_for_genes(symbol_to_vector, ref_genes)
                if ref_vectors.shape[0] == 0:
                    print(
                        f"WARN: sample {sid} contains no reference genes present in embeddings; skipping",
                        file=sys.stderr,
                    )
                    continue
                # Skip clustering in differential mode
                sample_sig = compute_differential_signature(vectors, ref_vectors, n_requested=len(genes))
                groups = [vectors]  # Single group for consistency
                if verbose_output:
                    print(f"INFO: Differential mode for {sid} - input: {len(found)}, reference: {len(ref_found)}")
            else:
                # Standard mode: cluster and compute signatures
                groups = maybe_cluster_vectors(
                    vectors, enable_clustering, min_genes_for_clustering=cluster_min_genes, random_state=random_state
                )
                # For reporting, use the aggregate signature of all found genes as the sample signature
                sample_sig = compute_signature(vectors, n_requested=len(genes))

            # QC row
            if "qc" in steps_set:
                per_sample_qc_rows.append([
                    sid,
                    len(genes),
                    len(found),
                    len(missing),
                    float(sample_sig.coverage),
                    int(len(groups)),
                ])
                if verbose_output:
                    print(
                        f"STEP: qc done sample={sid} n_input={len(genes)} n_found={len(found)} n_missing={len(missing)} dims={d} clusters={len(groups)}"
                    )

            # Save signature JSON for the sample
            if "signature" in steps_set:
                sig_path = os.path.join(
                    output_dir, f"{clean_id}_{pipeline_short}_{sid}_signature.json"
                )
                sig_payload = {
                    "sample_id": sid,
                    "n_requested": len(genes),
                    "n_found": int(sample_sig.n_genes),
                    "coverage": float(sample_sig.coverage),
                    "d": int(d),
                    "mean": sample_sig.mean.astype(float).tolist(),
                    "std": sample_sig.std.astype(float).tolist(),
                }
                write_json(sig_path, sig_payload)
                _print_output(sig_path, "json", f"Sample signature for {sid}", verbose=verbose_output)
                artifacts.append(os.path.abspath(sig_path))
                if verbose_output:
                    print(
                        f"STEP: signature done sample={sid} d={d}"
                    )
                signature_files.append(os.path.abspath(sig_path))

            # Cohort matrix row
            cohort_signature_matrix.append([sid] + [float(x) for x in sample_sig.mean.tolist()])
            s["signature"] = sample_sig
            s["found"] = found
            s["missing"] = missing

        # Term matching per sample
        if term_signatures is not None and "term_matching" in steps_set and s.get("signature") is not None:
            sig: Signature = s["signature"]
            dims_mask = make_dims_mask(
                dims_filter, sig.mean, sig.std, q_low=q_low, cv_threshold=cv_threshold, z_threshold=z_threshold
            )
            
            # Find closest terms
            ranked = rank_terms(
                sig.mean,
                sig.std,
                term_signatures,
                dims_mask=dims_mask,
                top_k=top_k,
            )
            match_csv = os.path.join(
                output_dir, f"{clean_id}_{pipeline_short}_{sid}_term_matches.csv"
            )
            rows = [["term_id", "library", "term", "distance", "cosine", "n_found", "coverage"]]
            for t, dist, cos in ranked:
                rows.append([t.term_id, t.library, t.term, f"{dist:.6g}", f"{cos:.6g}", t.n_found, f"{t.coverage:.3g}"])
            _write_csv(match_csv, rows[0], rows[1:])
            _print_output(match_csv, "csv", f"Top term matches for {sid}", verbose=verbose_output)
            artifacts.append(os.path.abspath(match_csv))
            if verbose_output:
                print(
                    f"STEP: term_matching done sample={sid} top_k={min(top_k, len(ranked))}"
                )
            term_matching_files.append(os.path.abspath(match_csv))
            # Summarize top-3 terms for results doc
            preview = ", ".join([f"{t.term[:40]} (d={dist:.3g}, cos={cos:.3g})" for t, dist, cos in ranked[:3]])
            sample_match_summaries.append(f"{sid}: {preview}")
            
            # Find opposite terms if enabled
            if should_find_opposite_terms(entity_id, find_opposite_terms):
                opposite_ranked = rank_terms_most_opposite(
                    sig.mean,
                    sig.std,
                    term_signatures,
                    dims_mask=dims_mask,
                    top_k=top_k,
                )
                opposite_match_csv = os.path.join(
                    output_dir, f"{clean_id}_{pipeline_short}_{sid}_opposite_term_matches.csv"
                )
                opposite_rows = [["term_id", "library", "term", "distance", "cosine", "n_found", "coverage"]]
                for t, dist, cos in opposite_ranked:
                    opposite_rows.append([t.term_id, t.library, t.term, f"{dist:.6g}", f"{cos:.6g}", t.n_found, f"{t.coverage:.3g}"])
                _write_csv(opposite_match_csv, opposite_rows[0], opposite_rows[1:])
                _print_output(opposite_match_csv, "csv", f"Top opposite term matches for {sid}", verbose=verbose_output)
                artifacts.append(os.path.abspath(opposite_match_csv))
                if verbose_output:
                    print(
                        f"STEP: opposite_term_matching done sample={sid} top_k={min(top_k, len(opposite_ranked))}"
                    )
                # Summarize top-3 opposite terms for results doc
                opposite_preview = ", ".join([f"{t.term[:40]} (d={dist:.3g}, cos={cos:.3g})" for t, dist, cos in opposite_ranked[:3]])
                sample_match_summaries.append(f"{sid} (opposite): {opposite_preview}")

    # Write QC summary
    if per_sample_qc_rows and "qc" in steps_set:
        qc_csv = os.path.join(output_dir, f"{clean_id}_{pipeline_short}_qc.csv")
        _write_csv(
            qc_csv,
            ["sample_id", "n_input_genes", "n_found", "n_missing", "coverage", "n_clusters"],
            per_sample_qc_rows,
        )
        _print_output(qc_csv, "csv", "Per-sample QC summary", verbose=verbose_output)
        artifacts.append(os.path.abspath(qc_csv))
        
        # Shared QC methodology
        shared_method_qc = (
            "## Quality Control\n\n"
            "### Gene Embedding Coverage\n\n"
            "For each gene set, normalized gene symbols are matched against the pre-trained embedding dictionary. "
            "Genes with available embeddings are extracted as D-dimensional vectors, while genes absent from "
            "the embedding space are tracked separately. Coverage is computed as the fraction of input genes "
            "successfully mapped to embeddings (coverage = n_found / n_requested). Gene sets with zero coverage "
            "are excluded from downstream analysis.\n\n"
            "### Optional Clustering Analysis\n\n"
            "To detect functional subpopulations or modules within gene sets, optional clustering can be performed "
            "when enabled and when gene sets contain sufficient genes (≥cluster_min_genes threshold). "
            "The primary clustering method uses HDBSCAN, a density-based algorithm that automatically determines "
            "the number of clusters without requiring pre-specification. The minimum cluster size is set to "
            "max(5, 0.2 × n_genes), and noise points (label = -1) are excluded. "
            "If HDBSCAN is unavailable or produces no valid clusters, K-means clustering is used as a fallback. "
            "For K-means, the optimal K is selected from the range [2, min(6, n_genes)) by maximizing the "
            "silhouette score, which measures cluster separation quality (range: [-1, 1], higher values indicate "
            "better-defined clusters). Only solutions with silhouette score > 0.35 and minimum cluster size > 2 genes "
            "are accepted. When clustering conditions are not met or fail, all genes are treated as a single group.\n\n"
        )
        _write_shared_method_doc(clean_id, output_dir, "qc", shared_method_qc, verbose=verbose_output)
        
        # Run-specific QC parameters
        params_qc = (
            f"## Run-Specific Parameters\n\n"
            f"- **Embedding dimensions (D)**: {d}\n"
            f"- **Total genes in embedding dictionary**: {_n}\n"
            f"- **Clustering enabled**: {enable_clustering}\n"
            f"- **Minimum genes for clustering**: {cluster_min_genes}\n"
            f"- **Random state**: {random_state}\n\n"
        )
        _write_run_params_doc(output_dir, clean_id, pipeline_short, "qc", params_qc, verbose=verbose_output)
        # RESULT.md for QC
        n_s = len(per_sample_qc_rows)
        total_input = sum(int(r[1]) for r in per_sample_qc_rows)
        total_found = sum(int(r[2]) for r in per_sample_qc_rows)
        total_missing = sum(int(r[3]) for r in per_sample_qc_rows)
        avg_cov = (sum(float(r[4]) for r in per_sample_qc_rows) / max(1, n_s)) if n_s else 0.0
        min_cov = min((float(r[4]) for r in per_sample_qc_rows), default=0.0)
        max_cov = max((float(r[4]) for r in per_sample_qc_rows), default=0.0)
        total_clusters = sum(int(r[5]) for r in per_sample_qc_rows)
        
        # Shared QC result interpretation
        shared_result_qc = (
            "## Quality Control Result Interpretation\n\n"
            "### Embedding Coverage\n\n"
            "Embedding coverage indicates the fraction of input genes successfully mapped to the pre-trained "
            "embedding space. High coverage (>80%) suggests the input gene set is well-represented in the "
            "embedding model. Low coverage may indicate species mismatch, non-standard gene nomenclature, "
            "or highly specialized genes not included in the training corpus.\n\n"
            "### Clustering Analysis\n\n"
            "When clustering is enabled and identifies multiple subgroups within a gene set, this suggests "
            "functional heterogeneity - the input genes may represent distinct biological processes or pathways. "
            "Conversely, when no clustering structure is detected (all genes form a single group), this indicates "
            "functional homogeneity, suggesting the genes are coherently related to a unified biological function.\n\n"
        )
        _write_shared_result_doc(clean_id, output_dir, "qc", shared_result_qc, verbose=verbose_output)
        
        # Run-specific QC results
        result_text_qc = (
            f"## Quality Control Results\n\n"
            f"### Gene Coverage\n\n"
            f"- **Total input genes**: {total_input}\n"
            f"- **Genes found in embeddings**: {total_found} ({(total_found / max(1, total_input)) * 100:.1f}%)\n"
            f"- **Genes not found**: {total_missing}\n"
        )
        
        # Add coverage statistics for multi-sample scenarios
        if n_s > 1:
            result_text_qc += (
                f"- **Coverage range across {n_s} samples**: {min_cov * 100:.1f}% - {max_cov * 100:.1f}%\n"
                f"- **Mean coverage**: {avg_cov * 100:.1f}%\n"
            )
        
        # Add clustering results
        result_text_qc += f"\n### Clustering Analysis\n\n"
        if enable_clustering and total_clusters > n_s:
            result_text_qc += (
                f"- **Total clusters detected**: {total_clusters} across {n_s} sample(s)\n"
                f"- **Average clusters per sample**: {(total_clusters / max(1, n_s)):.1f}\n"
                f"- **Interpretation**: Functional heterogeneity detected\n\n"
            )
        elif enable_clustering:
            result_text_qc += (
                f"- **Clusters detected**: {total_clusters} (no multi-cluster structure)\n"
                f"- **Interpretation**: Functional homogeneity\n\n"
            )
        else:
            result_text_qc += f"- **Status**: Clustering disabled\n\n"
        
        _write_run_results_doc(output_dir, clean_id, pipeline_short, "qc", result_text_qc, verbose=verbose_output)

    # Cohort integration (matrix)
    if cohort_signature_matrix and "cohort_integration" in steps_set:
        cohort_csv = os.path.join(output_dir, f"{clean_id}_{pipeline_short}_cohort_signatures.csv")
        _write_csv(cohort_csv, cohort_header, cohort_signature_matrix)
        _print_output(cohort_csv, "csv", "Cohort signature means (rows=samples, cols=dims)", verbose=verbose_output)
        artifacts.append(os.path.abspath(cohort_csv))
        if verbose_output:
            print(
                f"STEP: cohort_integration done n_samples={len(cohort_signature_matrix)} dims={d}"
            )
        # Shared cohort integration methodology
        shared_method_ci = (
            "## Multi-Sample Integration\n\n"
            "### Signature Matrix Construction\n\n"
            "For cohort-level analyses involving multiple samples, individual gene set signatures are aggregated "
            "into a unified matrix M ∈ ℝ^(n×D), where n is the number of samples "
            "and D is the embedding dimensionality. Each row M[i, :] corresponds to the mean embedding vector "
            "(μ) of sample i, creating a compact representation that preserves the geometric relationships in "
            "embedding space.\n\n"
            "### Normalization Considerations\n\n"
            "No additional normalization is applied at this stage, as the gene embeddings are pre-trained and "
            "pre-normalized, all samples are computed in the same embedding space, and the variance-weighted "
            "distance metric already accounts for scale differences across dimensions. This approach preserves "
            "the interpretability of the original embedding space coordinates. For visualization purposes "
            "(e.g., PCA, t-SNE, UMAP), post-hoc standardization can be applied if needed.\n\n"
        )
        _write_shared_method_doc(clean_id, output_dir, "cohort_integration", shared_method_ci, verbose=verbose_output)
        
        # Run-specific parameters
        params_ci = (
            f"## Run-Specific Parameters\n\n"
            f"- **Number of samples (n)**: {len(cohort_signature_matrix)}\n"
            f"- **Embedding dimensions (D)**: {d}\n"
            f"- **Matrix shape**: ({len(cohort_signature_matrix)}, {d})\n\n"
        )
        _write_run_params_doc(output_dir, clean_id, pipeline_short, "cohort_integration", params_ci, verbose=verbose_output)
        # Shared cohort integration result interpretation
        shared_result_ci = (
            "## Multi-Sample Integration Result Interpretation\n\n"
            "The cohort signature matrix serves as a unified representation of all samples in embedding space, "
            "enabling several downstream analyses:\n\n"
            "- **Dimensionality Reduction**: Apply PCA, t-SNE, or UMAP to visualize sample relationships in 2D/3D\n"
            "- **Clustering**: Identify groups of samples with similar gene set profiles\n"
            "- **Outlier Detection**: Detect samples with unusual embedding patterns\n"
            "- **Group Comparisons**: Compare centroids between experimental groups (e.g., case vs. control)\n"
            "- **Batch Effect Assessment**: Evaluate technical variation across batches or cohorts\n\n"
            "Matrix statistics (mean, SD, range) provide a summary of the overall embedding distribution "
            "and can be used to assess data quality and detect potential issues.\n\n"
        )
        _write_shared_result_doc(clean_id, output_dir, "cohort_integration", shared_result_ci, verbose=verbose_output)
        
        # Run-specific cohort integration results
        result_text_ci = (
            f"## Multi-Sample Integration Results\n\n"
            f"- **Number of samples**: {len(cohort_signature_matrix)}\n"
            f"- **Embedding dimensions**: {d}\n"
        )
        
        if cohort_signature_matrix and len(cohort_signature_matrix) > 1:
            matrix_data = np.array([row[1:] for row in cohort_signature_matrix], dtype=float)
            result_text_ci += (
                f"\n### Matrix Statistics\n\n"
                f"- **Mean embedding value**: {np.mean(matrix_data):.3f}\n"
                f"- **Standard deviation**: {np.std(matrix_data):.3f}\n"
                f"- **Range**: [{np.min(matrix_data):.3f}, {np.max(matrix_data):.3f}]\n"
            )
        
        result_text_ci += "\n"
        _write_run_results_doc(output_dir, clean_id, pipeline_short, "cohort_integration", result_text_ci, verbose=verbose_output)

    # METHOD/RESULT for signature step (after collecting signature files)
    if 'signature_files' in locals() and signature_files and "signature" in steps_set:
        # Shared signature methodology
        shared_method_sig = (
            "## Gene Set Signature Computation\n\n"
            "### Statistical Summary in Embedding Space\n\n"
            "For each gene set, a signature is computed as a statistical summary of the gene embedding vectors. "
            "Given an array of gene vectors V ∈ ℝ^(n×D), where n is the number of genes and D is the "
            "embedding dimensionality, we compute:\n\n"
            "1. **Mean vector (μ)**: The per-dimension mean across all gene vectors, representing the centroid "
            "of the gene set in embedding space:\n"
            "   ```\n"
            "   μ_d = (1/n) Σ_(i=1)^n v_i,d    for d ∈ [0, D-1]\n"
            "   ```\n\n"
            "2. **Standard deviation vector (σ)**: The per-dimension sample standard deviation with Bessel's correction:\n"
            "   ```\n"
            "   σ_d = sqrt( (1/(n-1)) Σ_(i=1)^n (v_i,d - μ_d)² )\n"
            "   ```\n\n"
            "The standard deviation captures the variability or coherence of the gene set in each embedding dimension. "
            "High σ in a dimension indicates genes are diverse in that semantic feature, while low σ indicates "
            "genes are coherent. To ensure numerical stability in downstream distance calculations, σ_d values "
            "below 1×10⁻⁸ are clamped to this minimum threshold.\n\n"
        )
        _write_shared_method_doc(clean_id, output_dir, "signature", shared_method_sig, verbose=verbose_output)
        
        # Run-specific parameters
        params_sig = (
            f"## Run-Specific Parameters\n\n"
            f"- **Embedding dimensions (D)**: {d}\n"
            f"- **Clustering enabled**: {enable_clustering}\n"
            f"- **Minimum genes for clustering**: {cluster_min_genes}\n"
            f"- **Random state**: {random_state}\n"
            f"- **Numerical stability threshold (ε)**: 1×10⁻⁸\n\n"
        )
        _write_run_params_doc(output_dir, clean_id, pipeline_short, "signature", params_sig, verbose=verbose_output)
        
        # Shared signature result interpretation
        shared_result_sig = (
            "## Gene Set Signature Result Interpretation\n\n"
            "Each gene set signature comprises two key components:\n\n"
            "- **Mean vector (μ)**: Represents the centroid or 'average position' of the gene set in embedding space. "
            "This vector captures the dominant functional characteristics shared by genes in the set.\n\n"
            "- **Standard deviation vector (σ)**: Captures within-set variability along each embedding dimension. "
            "Low σ values indicate genes are tightly clustered (functionally coherent) in that dimension, "
            "while high σ values suggest diversity or heterogeneity.\n\n"
            "Coverage indicates the fraction of input genes with available embeddings. High coverage (>80%) suggests "
            "comprehensive representation, while lower coverage may warrant inspection of missing genes.\n\n"
        )
        _write_shared_result_doc(clean_id, output_dir, "signature", shared_result_sig, verbose=verbose_output)
        
        # Calculate signature statistics across all samples
        n_sigs = len(signature_files)
        result_text_sig = (
            f"## Gene Set Signature Results\n\n"
            f"- **Signatures computed**: {n_sigs}\n"
            f"- **Embedding dimensions per signature**: {d}\n"
        )
        
        # Add aggregate statistics if we have per-sample data
        if samples:
            all_coverages = [s.get('signature').coverage for s in samples if s.get('signature') is not None]
            all_n_found = [s.get('signature').n_genes for s in samples if s.get('signature') is not None]
            if all_coverages:
                mean_cov = np.mean(all_coverages)
        result_text_sig += (
                    f"- **Average genes per signature**: {np.mean(all_n_found):.0f}\n"
                    f"- **Mean coverage**: {mean_cov * 100:.1f}%\n"
                )
        
        result_text_sig += "\n"
        _write_run_results_doc(output_dir, clean_id, pipeline_short, "signature", result_text_sig, verbose=verbose_output)

    # METHOD/RESULT for term matching step
    if 'term_matching_files' in locals() and term_matching_files and "term_matching" in steps_set:
        n_terms = len(term_signatures) if term_signatures else 0
        
        # Shared term matching methodology
        shared_method_tm = (
            "## Functional Term Matching\n\n"
            "### Overview\n\n"
            "To enable functional interpretation of gene set signatures, we identify the most similar precomputed "
            "reference signatures from curated gene set libraries (e.g., pathways, disease gene sets). Similarity "
            "is assessed using variance-weighted distance and cosine similarity metrics.\n\n"
            "### Variance-Weighted Distance\n\n"
            "The primary similarity metric weights each embedding dimension by the inverse of the combined variance "
            "between the query signature (S) and reference term signature (T):\n\n"
            "```\n"
            "D(S, T) = Σ_d w_d × (μ_S,d - μ_T,d)²\n"
            "```\n\n"
            "where μ_S,d and μ_T,d are the mean values at dimension d, σ_S,d and σ_T,d are the standard deviations, "
            "and w_d = 1 / (σ²_S,d + σ²_T,d + ε) with ε = 1×10⁻⁸ for numerical stability. "
            "This weighting scheme prioritizes dimensions with low variability (high confidence) and down-weights "
            "uncertain dimensions with high variability. Smaller distances indicate more similar signatures.\n\n"
            "### Cosine Similarity\n\n"
            "As a complementary metric, cosine similarity is computed to assess directional alignment:\n\n"
            "```\n"
            "cos(S, T) = (μ_S · μ_T) / (||μ_S|| × ||μ_T||)\n"
            "```\n\n"
            "Cosine similarity ranges from -1 (opposite) to +1 (perfectly aligned), with 0 indicating orthogonality. "
            "This metric captures directional similarity while ignoring magnitude differences.\n\n"
            "### Dimension Selection Methods\n\n"
            "To focus on the most stable and informative dimensions, several optional filtering methods are available:\n\n"
            "1. **Quantile Filter**: Retains dimensions with σ ≤ quantile(σ, q_low), selecting dimensions with lowest variability.\n\n"
            "2. **Coefficient of Variation (CV) Filter**: Computes CV_d = σ_d / (|μ_d| + ε) and retains dimensions with CV_d ≤ threshold.\n\n"
            "3. **Z-Score Filter**: Computes z_d = (σ_d - mean(σ)) / std(σ) and retains dimensions in the lower tail (z_d ≤ -threshold).\n\n"
            "4. **Top-K Filter**: Selects K dimensions with the lowest standard deviation.\n\n"
            "5. **No Filtering**: Uses all dimensions equally.\n\n"
            "### Term Ranking\n\n"
            "For each query signature, variance-weighted distances and cosine similarities are computed against "
            "all precomputed reference term signatures. Terms are ranked by distance (ascending), "
            "with cosine similarity used as a secondary criterion for tie-breaking. The top K most similar "
            "terms are retained for interpretation.\n\n"
        )
        _write_shared_method_doc(clean_id, output_dir, "term_matching", shared_method_tm, verbose=verbose_output)
        
        # Run-specific parameters (this is where library-specific info goes)
        params_tm = (
            f"## Run-Specific Parameters\n\n"
            f"### Reference Library\n\n"
            f"- **Number of precomputed terms**: {n_terms}\n"
        )
        if precomputed_terms_jsonl:
            library_name = os.path.basename(precomputed_terms_jsonl).replace('.signatures.jsonl', '')
            params_tm += f"- **Library**: {library_name}\n"
        
        params_tm += (
            f"\n### Dimension Selection\n\n"
            f"- **Method**: {dims_filter or 'none'}\n"
        )
        
        if dims_filter == "quantile":
            params_tm += f"- **Quantile threshold (q_low)**: {q_low} ({int(q_low * 100)}th percentile)\n"
        elif dims_filter == "cv":
            params_tm += f"- **CV threshold**: {cv_threshold}\n"
        elif dims_filter == "zscore":
            params_tm += f"- **Z-score threshold**: {z_threshold}\n"
        elif dims_filter == "topk":
            params_tm += f"- **Top-K dimensions**: {top_k}\n"
        
        params_tm += (
            f"\n### Ranking Parameters\n\n"
            f"- **Top-K matches returned**: {top_k}\n"
            f"- **Ranking metric**: Variance-weighted distance (primary), cosine similarity (tie-breaker)\n\n"
        )
        
        _write_run_params_doc(output_dir, clean_id, pipeline_short, "term_matching", params_tm, verbose=verbose_output)
        
        # Shared term matching result interpretation
        shared_result_tm = (
            "## Functional Term Matching Result Interpretation\n\n"
            "Term matching identifies known gene sets (pathways, diseases, functional categories) with similar "
            "embedding profiles to the query gene set. This provides biological context and functional interpretation.\n\n"
            "### Interpreting Match Metrics\n\n"
            "- **Distance (d)**: Variance-weighted distance. Lower values (typically 0-10) indicate stronger similarity. "
            "Distance <3 suggests very strong functional similarity, 3-6 moderate similarity, >6 weaker similarity.\n\n"
            "- **Cosine similarity (cos)**: Directional alignment. Values range from -1 to +1. "
            "Values >0.7 indicate strong alignment, 0.5-0.7 moderate alignment, <0.5 weak alignment.\n\n"
            "### Match Interpretation\n\n"
            "Top matches represent the closest functional relatives in the reference library. Multiple related terms "
            "in top matches suggest robust functional annotation. Diverse top matches may indicate:\n"
            "- The query gene set spans multiple biological processes\n"
            "- The gene set represents a novel or poorly characterized function\n"
            "- Low coverage or quality issues in the query gene set\n\n"
            "Coverage of matched terms indicates how well-represented those terms are in the embedding space.\n\n"
        )
        _write_shared_result_doc(clean_id, output_dir, "term_matching", shared_result_tm, verbose=verbose_output)
        
        # Run-specific term matching results
        library_name = "unknown"
        if precomputed_terms_jsonl:
            library_name = os.path.basename(precomputed_terms_jsonl).replace('.signatures.jsonl', '')
        
        result_text_tm = (
            f"## Functional Term Matching Results\n\n"
            f"- **Reference library**: {library_name}\n"
            f"- **Terms queried**: {n_terms}\n"
            f"- **Top matches returned per query**: {top_k}\n"
        )
        
        # Add preview of top matches if available
        if 'sample_match_summaries' in locals() and sample_match_summaries:
            result_text_tm += f"\n### Top Matches by Sample\n\n"
            for summary in sample_match_summaries[:5]:  # Show top 5 samples
                result_text_tm += f"- {summary}\n"
            if len(sample_match_summaries) > 5:
                result_text_tm += f"\n*Showing {min(5, len(sample_match_summaries))} of {len(sample_match_summaries)} samples*\n"
        
        result_text_tm += "\n"
        _write_run_results_doc(output_dir, clean_id, pipeline_short, "term_matching", result_text_tm, verbose=verbose_output)

    # Group comparison step not supported in single-sample mode

    # Run info
    run_info = {
        "entity_id": clean_id,
        "pipeline": pipeline_short,
        "output_dir": os.path.abspath(output_dir),
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "versions": _gather_versions(),
        "params": {
            "dims_filter": dims_filter,
            "q_low": q_low,
            "cv_threshold": cv_threshold,
            "z_threshold": z_threshold,
            "top_k": top_k,
            "enable_clustering": enable_clustering,
            "cluster_min_genes": cluster_min_genes,
            "random_state": random_state,
            "steps_to_run": sorted(list(steps_set)),
            "differential_mode": ref_genes_list is not None and len(ref_genes_list) > 0,
            "n_reference_genes": len(ref_genes_list) if ref_genes_list else 0,
        },
    }
    run_info_path = os.path.join(output_dir, f"{clean_id}_{pipeline_short}_run_info.json")
    write_json(run_info_path, run_info)
    _print_output(run_info_path, "json", "Run information and parameters", verbose=verbose_output)
    artifacts.append(os.path.abspath(run_info_path))

    # Generate figures if requested
    figure_paths: List[str] = []
    if generate_figures and samples and samples[0].get("signature") is not None:
        figures_dir = os.path.join(output_dir, "figures")
        os.makedirs(figures_dir, exist_ok=True)
        
        if figure_formats is None:
            figure_formats = ["svg", "png"]
        
        sample_sig = samples[0]["signature"]
        
        # Extract and clean library name for figure titles
        lib_name = ""
        if precomputed_terms_jsonl:
            lib_name = os.path.basename(precomputed_terms_jsonl).replace('.signatures.jsonl', '')
            # Import clean_library_name from generate_publication_figures
            from scripts.generate_publication_figures import clean_library_name
            lib_name = clean_library_name(lib_name)
        
        # Generate figures (descriptions already shown in initial output)
        figure_paths = _generate_figures(
            figures_dir, 
            sample_sig, 
            term_signatures, 
            figure_formats,
            clean_id,
            lib_name,
            verbose_output,
            find_opposite_terms
        )
        artifacts.extend(figure_paths)
        
        # Generate figure legend documentation
        legend_path = _write_figure_legend_doc(
            clean_id, 
            output_dir, 
            lib_name, 
            find_opposite_terms,
            verbose_output
        )
        artifacts.append(legend_path)

    # If verbose is off, print a summary line
    if not verbose_output:
        library_name = "unknown"
        if precomputed_terms_jsonl:
            library_name = os.path.basename(precomputed_terms_jsonl).replace('.signatures.jsonl', '')
        
        # Count all files in the output directory
        total_files = 0
        for root, dirs, files in os.walk(output_dir):
            total_files += len(files)
        
        print(f"✓ {total_files} result files for library {library_name} have been output to {os.path.abspath(output_dir)}")
    
    print(
        f"DONE: success=true output_dir={os.path.abspath(output_dir)} artifacts={len(artifacts)}"
    )
    return {"success": True, "run_info": run_info, "artifacts": artifacts, "figures": figure_paths}

