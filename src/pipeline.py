from __future__ import annotations

import json
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from .cluster import maybe_cluster_vectors
from .embeddings_io import get_vectors_for_genes, load_embeddings
from .match import TermSignature, rank_terms, rank_terms_farthest, rank_terms_most_opposite
from .signature import (
    Signature,
    compute_signature,
    select_dims_by_cv,
    select_dims_by_quantile,
    select_dims_by_zscore,
    select_top_k_low_variance,
)


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

    groups = maybe_cluster_vectors(vectors, enable_clustering, min_genes_for_clustering=cluster_min_genes)
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
    return rank_terms_farthest(input_signature.mean, input_signature.std, term_signatures, dims_mask=mask)


def match_input_to_terms_most_opposite(
    input_signature: Signature,
    term_signatures: List[TermSignature],
    dims_filter: Optional[str] = "quantile",
    **filter_kwargs,
):
    mask = make_dims_mask(dims_filter, input_signature.mean, input_signature.std, **filter_kwargs)
    return rank_terms_most_opposite(input_signature.mean, input_signature.std, term_signatures, dims_mask=mask)


