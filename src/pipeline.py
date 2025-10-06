from __future__ import annotations

import csv
import itertools
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
) -> Tuple[List[Signature], List[str], List[str], int]:
    symbol_to_vector, _n, d = load_embeddings(embeddings_path)
    vectors, found, missing = get_vectors_for_genes(symbol_to_vector, input_genes)
    if vectors.shape[0] == 0:
        raise ValueError("No input genes found in embeddings.")

    groups = maybe_cluster_vectors(
        vectors, enable_clustering, min_genes_for_clustering=cluster_min_genes
    )
    signatures: List[Signature] = [compute_signature(g, n_requested=len(input_genes)) for g in groups]
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


def generate_output_folder(entity_id: str, pipeline_short_name: str) -> str:
    return f"{entity_id}_{pipeline_short_name}"


def _get_workspace() -> str:
    workspace = os.getenv("PROJECT_WORKSPACE_VOLUME_PATH")
    if not workspace:
        # Fallback to current working directory with a warning for reproducibility.
        print(
            "WARN: PROJECT_WORKSPACE_VOLUME_PATH is not set; using current working directory",
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


def _write_method_doc(output_dir: str, entity_id: str, pipeline_short: str, step: str, text: str) -> str:
    path = os.path.join(output_dir, f"{entity_id}_{pipeline_short}_{step}.METHOD.md")
    _write_text(path, text)
    _print_output(path, "md", f"{step} methodology")
    return path


def _write_result_doc(output_dir: str, entity_id: str, pipeline_short: str, step: str, text: str) -> str:
    path = os.path.join(output_dir, f"{entity_id}_{pipeline_short}_{step}.RESULT.md")
    _write_text(path, text)
    _print_output(path, "md", f"{step} results")
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


def _parse_metadata_csv(path: str) -> Tuple[List[Dict[str, str]], List[str]]:
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        # Sanitize column names
        fieldnames = [re.sub(r"[^a-zA-Z0-9_-]", "", n) for n in (reader.fieldnames or [])]
        rows_raw = list(reader)
    rows: List[Dict[str, str]] = []
    for r in rows_raw:
        clean: Dict[str, str] = {}
        for k, v in r.items():
            ck = re.sub(r"[^a-zA-Z0-9_-]", "", k)
            clean[ck] = (v or "").strip()
        rows.append(clean)
    return rows, fieldnames


def _write_csv(path: str, header: List[str], rows: List[List[Any]]) -> None:
    ensure_dir(path)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for r in rows:
            writer.writerow(r)


def _print_output(path: str, ftype: str, desc: str) -> None:
    abs_path = os.path.abspath(path)
    print(f"OUTPUT: file={abs_path} type={ftype} desc=\"{desc}\"")


def geneset_embedding(
    embeddings_path: str,
    entity_id: str,
    output_folder: Optional[str] = None,
    input_genes: Optional[Iterable[str]] = None,
    input_genes_file: Optional[str] = None,
    precomputed_terms_jsonl: Optional[str] = None,
    metadata_csv: Optional[str] = None,
    grouping_column: Optional[str] = None,
    subgroup_column: Optional[str] = None,
    comparison_prefix: Optional[str] = None,
    steps_to_run: Optional[str] = None,
    # Matching and signature options
    dims_filter: Optional[str] = "quantile",
    q_low: float = 0.3,
    cv_threshold: float = 1.0,
    z_threshold: float = 0.5,
    top_k: int = 25,
    enable_clustering: bool = True,
    cluster_min_genes: int = 20,
    random_state: int = 17,
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
    folder_name = output_folder or generate_output_folder(clean_id, pipeline_short)
    output_dir = get_output_path(folder_name)
    os.makedirs(output_dir, exist_ok=True)

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

    # Prepare inputs: either single sample or metadata-defined cohort
    samples: List[Dict[str, Any]] = []
    artifacts: List[str] = []
    versions = _gather_versions()

    if metadata_csv:
        rows, fields = _parse_metadata_csv(metadata_csv)
        required_cols = ["sample_id", "filepath"]
        for req in required_cols:
            if req not in fields:
                raise ValueError(
                    f"ERROR: missing required column; field=metadata_csv expected includes {required_cols}"
                )
        for r in rows:
            sid = sanitize_simple_id(r.get("sample_id", ""))
            if not sid:
                continue
            samples.append({
                "sample_id": sid,
                "filepath": r.get("filepath", ""),
                "metadata": r,
            })
    else:
        # Single-sample mode
        if input_genes is None and not input_genes_file:
            raise ValueError(
                "ERROR: no input genes provided; field=input_genes|input_genes_file expected non-empty"
            )
        if input_genes is None:
            input_genes = _read_gene_list_from_file(str(input_genes_file))
        samples.append({"sample_id": clean_id, "genes": list(input_genes)})

    # Documentation: preprocess step
    if "preprocess" in steps_set:
        method_text = (
            f"## Preprocessing Methodology\n\n"
            f"### Overview\n\n"
            f"The preprocessing step parses and validates input data, preparing sample records for downstream analysis. "
            f"This pipeline supports two input modes: (1) single-sample mode with a direct gene list, and (2) cohort mode with metadata CSV/TSV.\n\n"
            f"### Input Mode Details\n\n"
            f"**Single-Sample Mode:**\n"
            f"- Genes can be provided via `--input-genes` (comma-separated) or `--input-genes-file` (text file)\n"
            f"- File format: one gene per line, or comma/space-separated on each line\n"
            f"- The sample_id is set to the entity_id provided\n\n"
            f"**Cohort/Metadata Mode:**\n"
            f"- Metadata file provided via `--metadata-csv` parameter\n"
            f"- Required columns: `sample_id`, `filepath`\n"
            f"- Optional columns: `{grouping_column or 'N/A'}` (grouping), `{subgroup_column or 'N/A'}` (subgroup)\n"
            f"- Each row represents one sample; filepath points to a gene list file\n"
            f"- Samples with invalid or missing filepaths are skipped with warnings\n\n"
            f"### Gene Symbol Normalization\n\n"
            f"All gene symbols are normalized using the following procedure:\n"
            f"1. Whitespace trimming (leading/trailing spaces removed)\n"
            f"2. Uppercase conversion (e.g., 'brca1' → 'BRCA1')\n"
            f"3. Empty strings are discarded\n"
            f"4. Malformed entries (non-alphanumeric special cases) are filtered during embedding lookup\n\n"
            f"This normalization is implemented in `embeddings_io._normalize_symbol()`.\n\n"
            f"### Run Configuration\n\n"
            f"- **entity_id**: {clean_id}\n"
            f"- **pipeline**: {pipeline_short}\n"
            f"- **n_samples**: {len(samples)}\n"
            f"- **output_directory**: {os.path.abspath(output_dir)}\n"
            f"- **embeddings_path**: {os.path.abspath(embeddings_path)}\n"
            f"- **precomputed_terms_jsonl**: {os.path.abspath(precomputed_terms_jsonl) if precomputed_terms_jsonl else 'None'}\n"
            f"- **metadata_csv**: {os.path.abspath(metadata_csv) if metadata_csv else 'None'}\n"
            f"- **grouping_column**: {grouping_column or 'None'}\n"
            f"- **subgroup_column**: {subgroup_column or 'None'}\n"
            f"- **comparison_prefix**: {comparison_prefix or 'comparison'}\n"
            f"- **steps_to_run**: {', '.join(sorted(steps_set))}\n"
            f"- **random_state**: {random_state}\n\n"
            f"### Software Environment\n\n"
            f"```json\n{json.dumps(versions, indent=2)}\n```\n\n"
            f"### File Format Specifications\n\n"
            f"**Embeddings File Format:**\n"
            f"- First line: `<num_genes> <num_dimensions>`\n"
            f"- Subsequent lines: `<GENE_SYMBOL> <float_1> <float_2> ... <float_D>`\n"
            f"- Whitespace-delimited (spaces or tabs)\n"
            f"- Malformed lines are silently skipped for robustness\n\n"
            f"**Precomputed Terms JSONL Format:**\n"
            f"- One JSON object per line\n"
            f"- Required fields: term_id, library, term, d, mean (list), std (list), n_found, coverage\n"
            f"- Generated by `precompute_terms.precompute_from_file()` or `precompute_from_dir()`\n\n"
        )
        _write_method_doc(output_dir, clean_id, pipeline_short, "preprocess", method_text)

        sample_details = "\n".join([f"  - {s['sample_id']}: {len(s.get('genes', []))} genes (initial)" for s in samples[:20]])
        if len(samples) > 20:
            sample_details += f"\n  ... and {len(samples) - 20} more samples"
        
        result_text = (
            f"## Preprocessing Results\n\n"
            f"### Summary\n\n"
            f"Successfully prepared {len(samples)} sample(s) for analysis. "
            f"In metadata mode, gene lists are loaded from individual filepaths during the QC step.\n\n"
            f"### Sample Details\n\n"
            f"{sample_details}\n\n"
            f"### Input Validation\n\n"
            f"- All sample IDs were sanitized (alphanumeric, underscores, hyphens only)\n"
            f"- Entity ID validated and cleaned: {clean_id}\n"
            f"- Output directory created: {os.path.abspath(output_dir)}\n\n"
        )
        _write_result_doc(output_dir, clean_id, pipeline_short, "preprocess", result_text)

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

    for s in samples:
        sid = s["sample_id"]
        # Load genes if provided by file path in metadata
        if "genes" not in s:
            gfile = s.get("filepath", "")
            if not gfile or not os.path.exists(gfile):
                print(
                    f"WARN: sample {sid} missing or invalid filepath; skipping",
                    file=sys.stderr,
                )
                continue
            s["genes"] = _read_gene_list_from_file(gfile)

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
                _print_output(sig_path, "json", f"Sample signature for {sid}")
                artifacts.append(os.path.abspath(sig_path))
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
            ranked = rank_terms(
                sig.mean,
                sig.std,
                term_signatures,
                dims_mask=make_dims_mask(
                    dims_filter, sig.mean, sig.std, q_low=q_low, cv_threshold=cv_threshold, z_threshold=z_threshold
                ),
                top_k=top_k,
            )
            match_csv = os.path.join(
                output_dir, f"{clean_id}_{pipeline_short}_{sid}_term_matches.csv"
            )
            rows = [["term_id", "library", "term", "distance", "cosine", "n_found", "coverage"]]
            for t, dist, cos in ranked:
                rows.append([t.term_id, t.library, t.term, f"{dist:.6g}", f"{cos:.6g}", t.n_found, f"{t.coverage:.3g}"])
            _write_csv(match_csv, rows[0], rows[1:])
            _print_output(match_csv, "csv", f"Top term matches for {sid}")
            artifacts.append(os.path.abspath(match_csv))
            print(
                f"STEP: term_matching done sample={sid} top_k={min(top_k, len(ranked))}"
            )
            term_matching_files.append(os.path.abspath(match_csv))
            # Summarize top-3 terms for results doc
            preview = ", ".join([f"{t.term[:40]} (d={dist:.3g}, cos={cos:.3g})" for t, dist, cos in ranked[:3]])
            sample_match_summaries.append(f"{sid}: {preview}")

    # Write QC summary
    if per_sample_qc_rows and "qc" in steps_set:
        qc_csv = os.path.join(output_dir, f"{clean_id}_{pipeline_short}_qc.csv")
        _write_csv(
            qc_csv,
            ["sample_id", "n_input_genes", "n_found", "n_missing", "coverage", "n_clusters"],
            per_sample_qc_rows,
        )
        _print_output(qc_csv, "csv", "Per-sample QC summary")
        artifacts.append(os.path.abspath(qc_csv))
        # METHOD.md for QC
        method_text_qc = (
            f"## Quality Control Methodology\n\n"
            f"### Overview\n\n"
            f"Quality control assesses the availability and quality of gene embeddings for each sample. "
            f"This step validates that sufficient genes are present in the embedding space and optionally identifies "
            f"distinct gene clusters within each sample.\n\n"
            f"### Gene Embedding Lookup\n\n"
            f"For each sample, genes are matched against the embedding dictionary:\n"
            f"1. Gene symbols are normalized (uppercase, trimmed) via `embeddings_io._normalize_symbol()`\n"
            f"2. Lookup in embedding dictionary `symbol_to_vector`\n"
            f"3. Found genes → vectors extracted and stacked into numpy array (shape: [n_found, D])\n"
            f"4. Missing genes → tracked separately for QC reporting\n\n"
            f"Implementation: `embeddings_io.get_vectors_for_genes()`\n\n"
            f"### Coverage Metric\n\n"
            f"Coverage = n_found / n_requested, where:\n"
            f"- n_found: number of genes with embeddings\n"
            f"- n_requested: total input genes for the sample\n\n"
            f"Samples with zero coverage are skipped with warnings.\n\n"
            f"### Clustering (Optional)\n\n"
            f"**Purpose:** Detect subpopulations or functional modules within gene lists.\n\n"
            f"**Algorithm Selection:**\n"
            f"- If `enable_clustering=True` and `n_genes >= cluster_min_genes`:\n"
            f"  - Primary: HDBSCAN (density-based, no preset K required)\n"
            f"    - min_cluster_size = max(5, 0.2 × n_genes)\n"
            f"    - Noise points (label = -1) are excluded\n"
            f"  - Fallback: KMeans if HDBSCAN unavailable or no valid clusters\n"
            f"    - Tests K ∈ [2, min(6, n_genes))\n"
            f"    - Selects K maximizing silhouette score\n"
            f"    - Requires silhouette score > 0.35 and cluster size > 2 genes\n"
            f"- If conditions not met, all genes treated as one cluster\n\n"
            f"Implementation: `cluster.maybe_cluster_vectors()`\n\n"
            f"**Silhouette Score:** Measures cluster quality, range [-1, 1]. Higher = better separated clusters.\n\n"
            f"### QC Metrics\n\n"
            f"Per sample, the following metrics are computed:\n"
            f"- **n_input_genes**: original input gene count\n"
            f"- **n_found**: genes with embeddings\n"
            f"- **n_missing**: genes without embeddings\n"
            f"- **coverage**: n_found / n_input_genes\n"
            f"- **n_clusters**: number of detected clusters (1 if clustering disabled/failed)\n\n"
            f"### Parameters\n\n"
            f"- **enable_clustering**: {enable_clustering}\n"
            f"- **cluster_min_genes**: {cluster_min_genes}\n"
            f"- **random_state**: {random_state}\n"
            f"- **embedding_dimensions**: {d}\n"
            f"- **total_genes_in_embeddings**: {_n}\n\n"
        )
        _write_method_doc(output_dir, clean_id, pipeline_short, "qc", method_text_qc)
        # RESULT.md for QC
        n_s = len(per_sample_qc_rows)
        total_input = sum(int(r[1]) for r in per_sample_qc_rows)
        total_found = sum(int(r[2]) for r in per_sample_qc_rows)
        total_missing = sum(int(r[3]) for r in per_sample_qc_rows)
        avg_cov = (sum(float(r[4]) for r in per_sample_qc_rows) / max(1, n_s)) if n_s else 0.0
        min_cov = min((float(r[4]) for r in per_sample_qc_rows), default=0.0)
        max_cov = max((float(r[4]) for r in per_sample_qc_rows), default=0.0)
        total_clusters = sum(int(r[5]) for r in per_sample_qc_rows)
        
        result_text_qc = (
            f"## Quality Control Results\n\n"
            f"### Summary Statistics\n\n"
            f"- **Samples processed**: {n_s}\n"
            f"- **Total input genes**: {total_input}\n"
            f"- **Total found in embeddings**: {total_found}\n"
            f"- **Total missing from embeddings**: {total_missing}\n"
            f"- **Global coverage**: {(total_found / max(1, total_input)):.4f}\n\n"
            f"### Coverage Statistics\n\n"
            f"- **Average coverage**: {avg_cov:.4f}\n"
            f"- **Min coverage**: {min_cov:.4f}\n"
            f"- **Max coverage**: {max_cov:.4f}\n\n"
            f"### Clustering Results\n\n"
            f"- **Total clusters detected**: {total_clusters} across {n_s} samples\n"
            f"- **Average clusters per sample**: {(total_clusters / max(1, n_s)):.2f}\n\n"
            f"### Per-Sample Details\n\n"
            f"Detailed QC metrics written to: {os.path.abspath(qc_csv)}\n\n"
            f"Columns:\n"
            f"- sample_id: Sample identifier\n"
            f"- n_input_genes: Input gene count\n"
            f"- n_found: Genes found in embeddings\n"
            f"- n_missing: Genes not in embeddings\n"
            f"- coverage: Fraction of genes found\n"
            f"- n_clusters: Number of detected clusters\n\n"
        )
        _write_result_doc(output_dir, clean_id, pipeline_short, "qc", result_text_qc)

    # Cohort integration (matrix)
    if cohort_signature_matrix and "cohort_integration" in steps_set:
        cohort_csv = os.path.join(output_dir, f"{clean_id}_{pipeline_short}_cohort_signatures.csv")
        _write_csv(cohort_csv, cohort_header, cohort_signature_matrix)
        _print_output(cohort_csv, "csv", "Cohort signature means (rows=samples, cols=dims)")
        artifacts.append(os.path.abspath(cohort_csv))
        print(
            f"STEP: cohort_integration done n_samples={len(cohort_signature_matrix)} dims={d}"
        )
        # METHOD.md for cohort integration
        method_text_ci = (
            f"## Cohort Integration Methodology\n\n"
            f"### Overview\n\n"
            f"Cohort integration aggregates per-sample signatures into a unified matrix representation, "
            f"enabling cross-sample comparisons, visualization (PCA, UMAP), and group-level analyses.\n\n"
            f"### Matrix Construction\n\n"
            f"**Structure:**\n"
            f"```\n"
            f"Cohort Matrix M ∈ ℝ^(n_samples × D)\n"
            f"M[i, :] = mean vector of sample i\n"
            f"```\n\n"
            f"where D = {d} (embedding dimensions)\n\n"
            f"**Construction Steps:**\n"
            f"1. For each sample, extract the signature mean vector μ ∈ ℝ^D\n"
            f"2. Stack mean vectors row-wise\n"
            f"3. Add sample_id as the first column for identification\n\n"
            f"**CSV Format:**\n"
            f"- Row 0 (header): sample_id, dim_0, dim_1, ..., dim_{d-1}\n"
            f"- Row i (data): <sample_id>, <mean_0>, <mean_1>, ..., <mean_{d-1}>\n\n"
            f"### Normalization Policy\n\n"
            f"**No additional normalization applied** at this stage because:\n"
            f"1. Gene embeddings are pre-trained and already normalized\n"
            f"2. Sample signatures are computed using the same embedding space\n"
            f"3. Variance-weighted distance already accounts for scale differences\n"
            f"4. Preserves interpretability of embedding space coordinates\n\n"
            f"If batch correction or standardization is needed for visualization, it should be applied "
            f"post-hoc (e.g., StandardScaler before PCA).\n\n"
            f"### Use Cases\n\n"
            f"The cohort matrix enables:\n"
            f"- **Dimensionality Reduction**: PCA, t-SNE, UMAP for visualization\n"
            f"- **Clustering**: K-means, hierarchical clustering across samples\n"
            f"- **Outlier Detection**: Identifying samples with unusual embedding profiles\n"
            f"- **Group Comparisons**: Centroid-based comparisons between metadata groups\n\n"
            f"### Parameters\n\n"
            f"- **n_samples**: {len(cohort_signature_matrix)}\n"
            f"- **n_dimensions**: {d}\n"
            f"- **matrix_shape**: ({len(cohort_signature_matrix)}, {d})\n\n"
        )
        _write_method_doc(output_dir, clean_id, pipeline_short, "cohort_integration", method_text_ci)
        # RESULT.md for cohort integration
        result_text_ci = (
            f"## Cohort Integration Results\n\n"
            f"### Summary\n\n"
            f"Successfully created cohort signature matrix.\n\n"
            f"### Matrix Details\n\n"
            f"- **Output file**: {os.path.abspath(cohort_csv)}\n"
            f"- **Samples**: {len(cohort_signature_matrix)}\n"
            f"- **Dimensions**: {d}\n"
            f"- **File format**: CSV (comma-separated values)\n"
            f"- **File size**: {os.path.getsize(cohort_csv) / 1024:.2f} KB\n\n"
            f"### Matrix Statistics\n\n"
        )
        if cohort_signature_matrix:
            matrix_data = np.array([row[1:] for row in cohort_signature_matrix], dtype=float)
            result_text_ci += (
                f"- **Global mean**: {np.mean(matrix_data):.6f}\n"
                f"- **Global std**: {np.std(matrix_data):.6f}\n"
                f"- **Global min**: {np.min(matrix_data):.6f}\n"
                f"- **Global max**: {np.max(matrix_data):.6f}\n\n"
            )
        result_text_ci += (
            f"### Sample IDs\n\n"
        )
        for i, row in enumerate(cohort_signature_matrix[:20], 1):
            result_text_ci += f"{i}. {row[0]}\n"
        if len(cohort_signature_matrix) > 20:
            result_text_ci += f"... and {len(cohort_signature_matrix) - 20} more samples\n"
        result_text_ci += "\n"
        _write_result_doc(output_dir, clean_id, pipeline_short, "cohort_integration", result_text_ci)

    # METHOD/RESULT for signature step (after collecting signature files)
    if 'signature_files' in locals() and signature_files and "signature" in steps_set:
        method_text_sig = (
            f"## Signature Computation Methodology\n\n"
            f"### Overview\n\n"
            f"A gene set signature is a statistical summary of the gene vectors in embedding space. "
            f"For each sample, we compute per-dimension mean and standard deviation, optionally after clustering.\n\n"
            f"### Algorithm\n\n"
            f"**Input:** Array of gene vectors V ∈ ℝ^(n×D), where n = number of genes, D = embedding dimensions\n\n"
            f"**Steps:**\n"
            f"1. **Optional Clustering** (see QC methodology for details):\n"
            f"   - If enabled and n >= cluster_min_genes, cluster vectors into groups\n"
            f"   - Each cluster gets its own signature (currently aggregated)\n"
            f"   - If disabled, all vectors form one group\n\n"
            f"2. **Per-Dimension Mean:**\n"
            f"   ```\n"
            f"   μ_d = (1/n) Σ_(i=1)^n v_i,d    for d ∈ [0, D-1]\n"
            f"   ```\n"
            f"   Implementation: `numpy.mean(vectors, axis=0)`\n\n"
            f"3. **Per-Dimension Sample Standard Deviation:**\n"
            f"   ```\n"
            f"   σ_d = sqrt( (1/(n-1)) Σ_(i=1)^n (v_i,d - μ_d)² )\n"
            f"   ```\n"
            f"   Implementation: `numpy.std(vectors, axis=0, ddof=1)` (Bessel's correction)\n\n"
            f"4. **Numerical Stability Clamping:**\n"
            f"   ```\n"
            f"   σ_d = max(σ_d, ε)    where ε = {1e-8}\n"
            f"   ```\n"
            f"   This prevents division by zero in downstream distance calculations.\n\n"
            f"5. **Coverage Calculation:**\n"
            f"   ```\n"
            f"   coverage = n_found / n_requested\n"
            f"   ```\n"
            f"   where n_found = genes with embeddings, n_requested = total input genes\n\n"
            f"Implementation: `signature.compute_signature()`\n\n"
            f"### Output Format\n\n"
            f"Each signature is saved as JSON with fields:\n"
            f"- **sample_id**: Sample identifier\n"
            f"- **n_requested**: Total genes in input list\n"
            f"- **n_found**: Genes with embeddings\n"
            f"- **coverage**: Fraction of genes found (n_found / n_requested)\n"
            f"- **d**: Number of embedding dimensions\n"
            f"- **mean**: Array of D mean values (one per dimension)\n"
            f"- **std**: Array of D standard deviation values (one per dimension)\n\n"
            f"### Interpretation\n\n"
            f"- **Mean vector (μ)**: Represents the \"centroid\" of the gene set in embedding space\n"
            f"- **Std vector (σ)**: Captures the spread/variability in each dimension\n"
            f"- High σ in a dimension → genes diverse in that semantic feature\n"
            f"- Low σ in a dimension → genes coherent/similar in that feature\n\n"
            f"### Parameters\n\n"
            f"- **enable_clustering**: {enable_clustering}\n"
            f"- **cluster_min_genes**: {cluster_min_genes}\n"
            f"- **random_state**: {random_state}\n"
            f"- **epsilon (ε)**: 1e-8 (numerical stability threshold)\n"
            f"- **embedding_dimensions (D)**: {d}\n\n"
        )
        _write_method_doc(output_dir, clean_id, pipeline_short, "signature", method_text_sig)
        result_text_sig = (
            f"## Signature Computation Results\n\n"
            f"### Summary\n\n"
            f"Successfully computed and saved {len(signature_files)} sample signature(s).\n\n"
            f"### Output Files\n\n"
        )
        for i, sf in enumerate(signature_files[:10], 1):
            result_text_sig += f"{i}. `{sf}`\n"
        if len(signature_files) > 10:
            result_text_sig += f"... and {len(signature_files) - 10} more files\n"
        result_text_sig += (
            f"\n### File Format\n\n"
            f"Each JSON file contains:\n"
            f"```json\n"
            f"{{\n"
            f'  "sample_id": "...",\n'
            f'  "n_requested": <int>,\n'
            f'  "n_found": <int>,\n'
            f'  "coverage": <float>,\n'
            f'  "d": {d},\n'
            f'  "mean": [<float>, ...],\n'
            f'  "std": [<float>, ...]\n'
            f"}}\n"
            f"```\n\n"
        )
        _write_result_doc(output_dir, clean_id, pipeline_short, "signature", result_text_sig)

    # METHOD/RESULT for term matching step
    if 'term_matching_files' in locals() and term_matching_files and "term_matching" in steps_set:
        n_terms = len(term_signatures) if term_signatures else 0
        method_text_tm = (
            f"## Term Matching Methodology\n\n"
            f"### Overview\n\n"
            f"Term matching identifies the most similar precomputed term signatures to each sample signature. "
            f"This enables functional interpretation by finding known gene sets (pathways, diseases, etc.) with similar embedding profiles.\n\n"
            f"### Distance Metric: Variance-Weighted Distance\n\n"
            f"**Formula:**\n"
            f"```\n"
            f"D(S, T) = Σ_d w_d × (μ_S,d - μ_T,d)²\n"
            f"```\n\n"
            f"where:\n"
            f"- S = sample signature, T = term signature\n"
            f"- μ_S,d, μ_T,d = mean values at dimension d\n"
            f"- σ_S,d, σ_T,d = standard deviations at dimension d\n"
            f"- w_d = 1 / (σ²_S,d + σ²_T,d + ε), ε = {1e-8}\n\n"
            f"**Interpretation:**\n"
            f"- Weights inversely proportional to combined variance\n"
            f"- Dimensions with low variability (high confidence) get more weight\n"
            f"- Dimensions with high variability (uncertain) get less weight\n"
            f"- Smaller distance → more similar signatures\n\n"
            f"Implementation: `match.variance_weighted_distance()`\n\n"
            f"### Cosine Similarity\n\n"
            f"**Formula:**\n"
            f"```\n"
            f"cos(S, T) = (μ_S · μ_T) / (||μ_S|| × ||μ_T||)\n"
            f"```\n\n"
            f"- Range: [-1, 1]\n"
            f"- 1 = perfectly aligned, 0 = orthogonal, -1 = opposite\n"
            f"- Measures directional similarity, ignoring magnitude\n\n"
            f"Implementation: `match.cosine_similarity()`\n\n"
            f"### Dimension Filtering (Optional)\n\n"
            f"To focus on the most informative dimensions, we optionally filter based on standard deviation.\n\n"
            f"**Method: {dims_filter or 'none'}**\n\n"
        )
        
        if dims_filter == "quantile":
            method_text_tm += (
                f"**Quantile Filter:**\n"
                f"- Keep dimensions where σ ≤ quantile(σ, q_low)\n"
                f"- q_low = {q_low} ({q_low * 100}th percentile)\n"
                f"- Selects dimensions with lowest variability (most stable)\n"
                f"- Implementation: `signature.select_dims_by_quantile()`\n\n"
            )
        elif dims_filter == "cv":
            method_text_tm += (
                f"**Coefficient of Variation (CV) Filter:**\n"
                f"- CV_d = σ_d / (|μ_d| + ε)\n"
                f"- Keep dimensions where CV_d ≤ {cv_threshold}\n"
                f"- Filters out dimensions with high relative variability\n"
                f"- Implementation: `signature.select_dims_by_cv()`\n\n"
            )
        elif dims_filter == "zscore":
            method_text_tm += (
                f"**Z-Score Filter:**\n"
                f"- Compute z_d = (σ_d - mean(σ)) / std(σ)\n"
                f"- Keep dimensions where z_d ≤ -{abs(z_threshold)}\n"
                f"- z_threshold = {z_threshold}\n"
                f"- Selects dimensions in the lower tail of variance distribution\n"
                f"- Implementation: `signature.select_dims_by_zscore()`\n\n"
            )
        elif dims_filter == "topk":
            method_text_tm += (
                f"**Top-K Filter:**\n"
                f"- Select K = {top_k} dimensions with lowest σ\n"
                f"- Fixed number of dimensions regardless of distribution\n"
                f"- Implementation: `signature.select_top_k_low_variance()`\n\n"
            )
        else:
            method_text_tm += (
                f"**No Filtering:**\n"
                f"- All {d} dimensions used equally\n\n"
            )
        
        method_text_tm += (
            f"### Ranking and Selection\n\n"
            f"1. Compute distance and cosine for all {n_terms} term signatures\n"
            f"2. Sort by distance (ascending), then by cosine (descending) as tie-breaker\n"
            f"3. Return top {top_k} matches\n\n"
            f"Implementation: `match.rank_terms()`\n\n"
            f"### Parameters Summary\n\n"
            f"- **dims_filter**: {dims_filter or 'none'}\n"
            f"- **q_low**: {q_low} (for quantile filter)\n"
            f"- **cv_threshold**: {cv_threshold} (for CV filter)\n"
            f"- **z_threshold**: {z_threshold} (for z-score filter)\n"
            f"- **top_k**: {top_k}\n"
            f"- **n_precomputed_terms**: {n_terms}\n"
            f"- **precomputed_terms_source**: {os.path.abspath(precomputed_terms_jsonl) if precomputed_terms_jsonl else 'N/A'}\n\n"
        )
        _write_method_doc(output_dir, clean_id, pipeline_short, "term_matching", method_text_tm)
        
        preview = "\n".join(sample_match_summaries[:10]) if 'sample_match_summaries' in locals() else ""
        result_text_tm = (
            f"## Term Matching Results\n\n"
            f"### Summary\n\n"
            f"Generated {len(term_matching_files)} per-sample match table(s), each containing top {top_k} term matches.\n\n"
            f"### Output Files\n\n"
        )
        for i, tmf in enumerate(term_matching_files[:10], 1):
            result_text_tm += f"{i}. `{tmf}`\n"
        if len(term_matching_files) > 10:
            result_text_tm += f"... and {len(term_matching_files) - 10} more files\n"
        
        result_text_tm += (
            f"\n### CSV Format\n\n"
            f"Each match file contains columns:\n"
            f"- **term_id**: Unique MD5 hash identifier for the term\n"
            f"- **library**: Source database/library (e.g., ClinVar_2019.tsv)\n"
            f"- **term**: Human-readable term name (pathway, disease, etc.)\n"
            f"- **distance**: Variance-weighted distance (lower = more similar)\n"
            f"- **cosine**: Cosine similarity (higher = more similar direction)\n"
            f"- **n_found**: Number of genes in term with embeddings\n"
            f"- **coverage**: Fraction of term genes with embeddings\n\n"
            f"### Top Matches Preview\n\n"
            f"First 10 samples:\n\n"
            f"{preview}\n\n"
        )
        _write_result_doc(output_dir, clean_id, pipeline_short, "term_matching", result_text_tm)

    # Group comparison (pairwise) using metadata groupings
    group_comp_files: List[str] = []
    if metadata_csv and grouping_column and "group_comparison" in steps_set and term_signatures is not None:
        rows, _fields = _parse_metadata_csv(metadata_csv)
        # Map sample -> group
        sid_to_group: Dict[str, str] = {}
        for r in rows:
            sid = sanitize_simple_id(r.get("sample_id", ""))
            grp = r.get(grouping_column, "")
            if sid and grp:
                sid_to_group[sid] = grp
        unique_groups = sorted(set(sid_to_group.values()))
        # Build group centroid signatures
        group_to_sigs: Dict[str, List[Signature]] = {}
        for s in samples:
            sig: Optional[Signature] = s.get("signature")  # type: ignore[assignment]
            if sig is None:
                continue
            sid = s["sample_id"]
            grp = sid_to_group.get(sid)
            if not grp:
                continue
            group_to_sigs.setdefault(grp, []).append(sig)

        def _centroid(sigs: List[Signature]) -> Signature:
            means = np.stack([x.mean for x in sigs], axis=0)
            stds = np.stack([x.std for x in sigs], axis=0)
            mean = means.mean(axis=0)
            # Aggregate std conservatively as mean std
            std = np.maximum(stds.mean(axis=0), 1e-8)
            return Signature(mean=mean.astype(np.float32), std=std.astype(np.float32), n_genes=0, coverage=1.0)

        centroids: Dict[str, Signature] = {g: _centroid(v) for g, v in group_to_sigs.items() if v}

        # Pairwise comparisons
        pair_count = 0
        for g1, g2 in itertools.combinations(sorted(centroids.keys()), 2):
            sig1 = centroids[g1]
            sig2 = centroids[g2]
            # Score each term by distance to centroid; compare delta = dist_g2 - dist_g1
            mask1 = make_dims_mask(dims_filter, sig1.mean, sig1.std, q_low=q_low, cv_threshold=cv_threshold, z_threshold=z_threshold)
            mask2 = make_dims_mask(dims_filter, sig2.mean, sig2.std, q_low=q_low, cv_threshold=cv_threshold, z_threshold=z_threshold)
            # Use mask from g1 for consistency if both are present and same size; otherwise no mask
            dims_mask = mask1 if (mask1 is not None and mask2 is not None and mask1.shape == mask2.shape) else None
            rows_out: List[List[Any]] = []
            header = ["term_id", "library", "term", f"dist_{g1}", f"dist_{g2}", "delta_g2_minus_g1"]
            for t in term_signatures:
                d1 = float(np.sum(((sig1.mean[dims_mask] if dims_mask is not None else sig1.mean) - (t.mean[dims_mask] if dims_mask is not None else t.mean)) ** 2))
                d2 = float(np.sum(((sig2.mean[dims_mask] if dims_mask is not None else sig2.mean) - (t.mean[dims_mask] if dims_mask is not None else t.mean)) ** 2))
                rows_out.append([t.term_id, t.library, t.term, f"{d1:.6g}", f"{d2:.6g}", f"{(d2 - d1):.6g}"])
            # Sort by absolute delta descending
            rows_out.sort(key=lambda r: -abs(float(r[-1])))
            comp_prefix = comparison_prefix or "comparison"
            comp_file = os.path.join(
                output_dir,
                f"{clean_id}_{pipeline_short}_{comp_prefix}_{g1}_vs_{g2}_differential.tsv",
            )
            ensure_dir(comp_file)
            with open(comp_file, "w", encoding="utf-8", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(header)
                for r in rows_out:
                    w.writerow(r)
            _print_output(comp_file, "tsv", f"Group comparison {g1} vs {g2}")
            artifacts.append(os.path.abspath(comp_file))
            pair_count += 1
            group_comp_files.append(os.path.abspath(comp_file))
        print(
            f"STEP: group_comparison done pairs={pair_count} groups={len(centroids)} grouping_column={grouping_column}"
        )
        # METHOD.md for group comparison
        method_text_gc = (
            f"## Group Comparison Methodology\n\n"
            f"### Overview\n\n"
            f"Group comparison identifies terms (pathways, diseases) differentially associated with metadata-defined "
            f"sample groups. This enables functional characterization of group differences.\n\n"
            f"### Centroid Computation\n\n"
            f"For each group, we compute a centroid signature representing the average embedding profile:\n\n"
            f"**Algorithm:**\n"
            f"1. Partition samples by grouping column: `{grouping_column}`\n"
            f"2. For each group G, collect sample signatures S_1, S_2, ..., S_n\n"
            f"3. Compute centroid mean: μ_G = (1/n) Σ_i μ_i\n"
            f"4. Compute centroid std: σ_G = (1/n) Σ_i σ_i (conservative aggregation)\n"
            f"5. Apply max(σ_G, 1e-8) for numerical stability\n\n"
            f"**Groups detected:**\n"
        )
        for grp in sorted(centroids.keys()):
            n_samples_in_grp = len(group_to_sigs.get(grp, []))
            method_text_gc += f"- **{grp}**: {n_samples_in_grp} samples\n"
        
        method_text_gc += (
            f"\n### Pairwise Differential Analysis\n\n"
            f"For each pair of groups (G1, G2), we identify terms showing preferential association:\n\n"
            f"**Distance Calculation:**\n"
            f"```\n"
            f"For each term T:\n"
            f"  d_G1 = Σ_d (μ_G1,d - μ_T,d)²    (squared Euclidean distance, optionally filtered)\n"
            f"  d_G2 = Σ_d (μ_G2,d - μ_T,d)²\n"
            f"  Δ = d_G2 - d_G1\n"
            f"```\n\n"
            f"**Interpretation of Δ:**\n"
            f"- Δ > 0: Term closer to G1 (G1-specific)\n"
            f"- Δ < 0: Term closer to G2 (G2-specific)\n"
            f"- |Δ| = magnitude of differential association\n\n"
            f"**Ranking:** Terms sorted by |Δ| (descending) to prioritize strong differential signals.\n\n"
            f"### Dimension Filtering\n\n"
            f"- **Method**: {dims_filter or 'none'}\n"
            f"- Applied consistently across both groups if filters are compatible\n"
            f"- Focuses comparison on stable, informative dimensions\n\n"
            f"### Parameters\n\n"
            f"- **grouping_column**: {grouping_column}\n"
            f"- **n_groups**: {len(centroids)}\n"
            f"- **n_comparisons**: {pair_count}\n"
            f"- **dims_filter**: {dims_filter or 'none'}\n"
            f"- **comparison_prefix**: {comparison_prefix or 'comparison'}\n\n"
        )
        _write_method_doc(output_dir, clean_id, pipeline_short, "group_comparison", method_text_gc)
        # RESULT.md for group comparison
        example_gc = group_comp_files[0] if group_comp_files else ""
        result_text_gc = (
            f"## Group Comparison Results\n\n"
            f"### Summary\n\n"
            f"Generated {len(group_comp_files)} pairwise differential tables for {len(centroids)} groups.\n\n"
            f"### Output Files\n\n"
        )
        for i, gcf in enumerate(group_comp_files, 1):
            result_text_gc += f"{i}. `{gcf}`\n"
        
        result_text_gc += (
            f"\n### TSV Format\n\n"
            f"Each file contains columns:\n"
            f"- **term_id**: Unique term identifier\n"
            f"- **library**: Source database\n"
            f"- **term**: Term name\n"
            f"- **dist_<group1>**: Distance from group1 centroid to term\n"
            f"- **dist_<group2>**: Distance from group2 centroid to term\n"
            f"- **delta_g2_minus_g1**: Δ = dist_group2 - dist_group1\n\n"
            f"### Interpretation Guide\n\n"
            f"**Top rows (largest |Δ|)** indicate terms most differentially associated with groups:\n"
            f"- **Positive Δ**: Terms closer to first group (first group-specific)\n"
            f"- **Negative Δ**: Terms closer to second group (second group-specific)\n\n"
            f"**Example Use Case:** If comparing diseased vs. healthy:\n"
            f"- Large positive Δ → pathway enriched in diseased\n"
            f"- Large negative Δ → pathway enriched in healthy\n\n"
        )
        _write_result_doc(output_dir, clean_id, pipeline_short, "group_comparison", result_text_gc)

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
        },
    }
    run_info_path = os.path.join(output_dir, f"{clean_id}_{pipeline_short}_run_info.json")
    write_json(run_info_path, run_info)
    _print_output(run_info_path, "json", "Run information and parameters")
    artifacts.append(os.path.abspath(run_info_path))

    print(
        f"DONE: success=true output_dir={os.path.abspath(output_dir)} artifacts={len(artifacts)}"
    )
    return {"success": True, "run_info": run_info, "artifacts": artifacts}

