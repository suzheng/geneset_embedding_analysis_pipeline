"""Embedding signature pipeline package."""

__all__ = [
    # Intentionally minimal to avoid heavy imports at package import time.
]


# Re-export main tool function for external use
from .pipeline import (
    compute_input_signatures,
    load_precomputed_terms_jsonl,
    make_dims_mask,
    match_input_to_terms,
    match_input_to_terms_farthest,
    match_input_to_terms_most_opposite,
    geneset_embedding,
)

__all__ += [
    "compute_input_signatures",
    "load_precomputed_terms_jsonl",
    "make_dims_mask",
    "match_input_to_terms",
    "match_input_to_terms_farthest",
    "match_input_to_terms_most_opposite",
    "geneset_embedding",
]


