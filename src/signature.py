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


def compute_differential_signature(
    input_vectors: np.ndarray, 
    reference_vectors: np.ndarray, 
    n_requested: int, 
    epsilon: float = DEFAULT_EPSILON
) -> Signature:
    """
    Compute a differential signature: (input_mean - reference_mean).
    
    The standard deviation is computed as the pooled std from both groups.
    This represents the signature for the differential embedding space.
    
    Args:
        input_vectors: Gene vectors for input genes (shape: [n_input, D])
        reference_vectors: Gene vectors for reference genes (shape: [n_ref, D])
        n_requested: Total number of input genes requested
        epsilon: Numerical stability threshold for std
    
    Returns:
        Signature with differential mean and pooled std
    """
    if input_vectors.size == 0:
        raise ValueError("No input vectors provided to compute differential signature.")
    if reference_vectors.size == 0:
        raise ValueError("No reference vectors provided to compute differential signature.")
    
    # Compute means for both groups
    input_mean = input_vectors.mean(axis=0)
    reference_mean = reference_vectors.mean(axis=0)
    
    # Differential mean
    diff_mean = input_mean - reference_mean
    
    # Compute pooled standard deviation
    # For two groups, pooled std = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
    n_input = input_vectors.shape[0]
    n_ref = reference_vectors.shape[0]
    
    input_std = input_vectors.std(axis=0, ddof=1)
    reference_std = reference_vectors.std(axis=0, ddof=1)
    
    # Pooled variance
    pooled_var = ((n_input - 1) * input_std**2 + (n_ref - 1) * reference_std**2) / (n_input + n_ref - 2)
    pooled_std = np.sqrt(pooled_var)
    
    # Clamp to epsilon for numerical stability
    pooled_std = np.where(pooled_std < epsilon, epsilon, pooled_std)
    
    # Coverage is based on input genes
    coverage = float(n_input) / float(max(1, n_requested))
    
    return Signature(
        mean=diff_mean.astype(np.float32), 
        std=pooled_std.astype(np.float32), 
        n_genes=n_input, 
        coverage=coverage
    )


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



