import numpy as np

from src.match import (
    TermSignature,
    cosine_similarity,
    rank_terms,
    rank_terms_farthest,
    rank_terms_most_opposite,
    variance_weighted_distance,
)


def test_variance_weighted_distance_and_cosine():
    s_mean = np.array([0.0, 0.0], dtype=np.float32)
    s_std = np.array([1.0, 1.0], dtype=np.float32)
    t_mean = np.array([1.0, 0.0], dtype=np.float32)
    t_std = np.array([1.0, 1.0], dtype=np.float32)

    d = variance_weighted_distance(s_mean, s_std, t_mean, t_std)
    assert d > 0

    cos = cosine_similarity(s_mean, t_mean)
    assert abs(cos - 0.0) < 1e-6


def test_rank_terms():
    s_mean = np.array([0.0, 0.0], dtype=np.float32)
    s_std = np.array([1.0, 1.0], dtype=np.float32)
    candidates = [
        TermSignature(term_id="a", library="L", term="T1", mean=np.array([0.0, 0.0], dtype=np.float32), std=np.array([1.0, 1.0], dtype=np.float32), n_found=10, coverage=1.0),
        TermSignature(term_id="b", library="L", term="T2", mean=np.array([2.0, 0.0], dtype=np.float32), std=np.array([1.0, 1.0], dtype=np.float32), n_found=10, coverage=1.0),
    ]
    ranked = rank_terms(s_mean, s_std, candidates, top_k=2)
    assert len(ranked) == 2
    # First candidate should be the zero vector (closer)
    assert ranked[0][0].term_id == "a"


def test_rank_terms_farthest_and_most_opposite():
    s_mean = np.array([1.0, 0.0], dtype=np.float32)
    s_std = np.array([1.0, 1.0], dtype=np.float32)
    candidates = [
        TermSignature(term_id="near", library="L", term="near", mean=np.array([1.2, 0.0], dtype=np.float32), std=np.array([1.0, 1.0], dtype=np.float32), n_found=5, coverage=1.0),
        TermSignature(term_id="far", library="L", term="far", mean=np.array([10.0, 0.0], dtype=np.float32), std=np.array([1.0, 1.0], dtype=np.float32), n_found=5, coverage=1.0),
        TermSignature(term_id="oppo", library="L", term="oppo", mean=np.array([-1.0, 0.0], dtype=np.float32), std=np.array([1.0, 1.0], dtype=np.float32), n_found=5, coverage=1.0),
    ]

    # Farthest by distance should prefer 'far'
    farthest = rank_terms_farthest(s_mean, s_std, candidates, top_k=2)
    assert farthest[0][0].term_id == "far"

    # Most opposite by cosine should prefer 'oppo'
    most_opposite = rank_terms_most_opposite(s_mean, s_std, candidates, top_k=2)
    assert most_opposite[0][0].term_id == "oppo"


