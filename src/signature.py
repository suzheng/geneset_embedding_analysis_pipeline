from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import numpy as np


DEFAULT_EPSILON: float = 1e-8


@dataclass
class Signature:
    mean: np.ndarray
    std: np.ndarray
    n_genes: int
    coverage: float


def compute_signature(
    vectors: np.ndarray, n_requested: int, epsilon: float = DEFAULT_EPSILON
) -> Signature:
    """
    Compute per-dimension mean and std for the provided gene vectors.
    """
    if vectors.size == 0:
        raise ValueError("No vectors provided to compute signature.")
    mean = vectors.mean(axis=0)
    std = vectors.std(axis=0, ddof=1)  # sample std
    std = np.where(std < epsilon, epsilon, std)
    n_genes = vectors.shape[0]
    coverage = float(n_genes) / float(max(1, n_requested))
    return Signature(mean=mean.astype(np.float32), std=std.astype(np.float32), n_genes=n_genes, coverage=coverage)


# ---- Optional dimension filtering utilities (std-based) ----

def select_dims_by_quantile(std: np.ndarray, q_low: float) -> np.ndarray:
    threshold = np.quantile(std, q_low)
    return std <= threshold


def select_dims_by_cv(mean: np.ndarray, std: np.ndarray, cv_threshold: float, epsilon: float = DEFAULT_EPSILON) -> np.ndarray:
    cv = std / (np.abs(mean) + epsilon)
    return cv <= cv_threshold


def select_dims_by_zscore(std: np.ndarray, z_threshold: float) -> np.ndarray:
    mu = float(np.mean(std))
    sigma = float(np.std(std))
    if sigma == 0.0:
        return np.ones_like(std, dtype=bool)
    z = (std - mu) / sigma
    return z <= -abs(z_threshold)


def select_top_k_low_variance(std: np.ndarray, k: int) -> np.ndarray:
    if k <= 0:
        return np.zeros_like(std, dtype=bool)
    if k >= std.size:
        return np.ones_like(std, dtype=bool)
    idx = np.argpartition(std, k)[:k]
    mask = np.zeros_like(std, dtype=bool)
    mask[idx] = True
    return mask


