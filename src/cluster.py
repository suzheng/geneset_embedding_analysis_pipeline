from __future__ import annotations

from typing import List, Optional

import numpy as np

try:
    import hdbscan  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    hdbscan = None

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def maybe_cluster_vectors(
    vectors: np.ndarray,
    enable_clustering: bool,
    min_genes_for_clustering: int = 20,
    random_state: int = 17,
) -> List[np.ndarray]:
    """
    Optionally cluster vectors into groups; return a list of arrays, one per cluster.
    If clustering is disabled or not enough vectors, returns [vectors].
    """
    n = vectors.shape[0]
    if not enable_clustering or n < min_genes_for_clustering:
        return [vectors]

    if hdbscan is not None and n >= max(10, min_genes_for_clustering):
        clusterer = hdbscan.HDBSCAN(min_cluster_size=max(5, int(0.2 * n)))
        labels = clusterer.fit_predict(vectors)
        unique = [lab for lab in sorted(set(labels)) if lab != -1]
        if not unique:
            return [vectors]
        clusters = [vectors[labels == lab] for lab in unique]
        return clusters

    # Fallback to KMeans with small grid of K
    best_k = None
    best_score = -1.0
    best_labels: Optional[np.ndarray] = None
    for k in range(2, min(6, n)):
        try:
            km = KMeans(n_clusters=k, n_init=10, random_state=random_state)
            labels = km.fit_predict(vectors)
            score = silhouette_score(vectors, labels)
        except Exception:
            continue
        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    if best_labels is None or best_score < 0.35:
        return [vectors]

    groups: List[np.ndarray] = [vectors[best_labels == i] for i in range(int(best_k))]
    # Filter tiny clusters (<= 2 items) back into a single group if needed
    groups = [g for g in groups if g.shape[0] > 2]
    if not groups:
        return [vectors]
    return groups


