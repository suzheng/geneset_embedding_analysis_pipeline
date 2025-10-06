from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


DEFAULT_EPSILON: float = 1e-8


@dataclass
class TermSignature:
    term_id: str
    library: str
    term: str
    mean: np.ndarray
    std: np.ndarray
    n_found: int
    coverage: float


def variance_weighted_distance(
    s_mean: np.ndarray,
    s_std: np.ndarray,
    t_mean: np.ndarray,
    t_std: np.ndarray,
    dims_mask: Optional[np.ndarray] = None,
    epsilon: float = DEFAULT_EPSILON,
) -> float:
    if dims_mask is not None:
        s_mean = s_mean[dims_mask]
        s_std = s_std[dims_mask]
        t_mean = t_mean[dims_mask]
        t_std = t_std[dims_mask]
    w = 1.0 / (np.square(s_std) + np.square(t_std) + epsilon)
    diff = s_mean - t_mean
    return float(np.sum(w * np.square(diff)))


def cosine_similarity(a: np.ndarray, b: np.ndarray, dims_mask: Optional[np.ndarray] = None) -> float:
    if dims_mask is not None:
        a = a[dims_mask]
        b = b[dims_mask]
    denom = (np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0.0:
        return 0.0
    return float(np.dot(a, b) / denom)


def rank_terms(
    s_mean: np.ndarray,
    s_std: np.ndarray,
    candidates: Sequence[TermSignature],
    dims_mask: Optional[np.ndarray] = None,
    top_k: int = 25,
) -> List[Tuple[TermSignature, float, float]]:
    """
    Return top_k candidates ranked by ascending variance-weighted distance.
    Also return cosine similarity as a secondary score.
    """
    scores: List[Tuple[TermSignature, float, float]] = []
    for t in candidates:
        d = variance_weighted_distance(s_mean, s_std, t.mean, t.std, dims_mask)
        cos = cosine_similarity(s_mean, t.mean, dims_mask)
        scores.append((t, d, cos))
    scores.sort(key=lambda x: (x[1], -x[2]))
    return scores[:top_k]


def rank_terms_farthest(
    s_mean: np.ndarray,
    s_std: np.ndarray,
    candidates: Sequence[TermSignature],
    dims_mask: Optional[np.ndarray] = None,
    top_k: int = 25,
) -> List[Tuple[TermSignature, float, float]]:
    """
    Return top_k candidates ranked by descending variance-weighted distance.
    Tie-breaker prefers more negative cosine (i.e., more opposite orientation).
    """
    scores: List[Tuple[TermSignature, float, float]] = []
    for t in candidates:
        d = variance_weighted_distance(s_mean, s_std, t.mean, t.std, dims_mask)
        cos = cosine_similarity(s_mean, t.mean, dims_mask)
        scores.append((t, d, cos))
    scores.sort(key=lambda x: (-x[1], x[2]))
    return scores[:top_k]


def rank_terms_most_opposite(
    s_mean: np.ndarray,
    s_std: np.ndarray,
    candidates: Sequence[TermSignature],
    dims_mask: Optional[np.ndarray] = None,
    top_k: int = 25,
) -> List[Tuple[TermSignature, float, float]]:
    """
    Return top_k candidates with the most negative cosine similarity.
    Tie-breaker prefers smaller distance (closer in magnitude but opposite direction).
    """
    scores: List[Tuple[TermSignature, float, float]] = []
    for t in candidates:
        d = variance_weighted_distance(s_mean, s_std, t.mean, t.std, dims_mask)
        cos = cosine_similarity(s_mean, t.mean, dims_mask)
        scores.append((t, d, cos))
    scores.sort(key=lambda x: (x[2], x[1]))  # cosine ascending (most negative first), then distance ascending
    return scores[:top_k]


