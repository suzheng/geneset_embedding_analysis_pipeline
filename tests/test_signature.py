import numpy as np

from src.signature import (
    compute_signature,
    select_dims_by_quantile,
    select_dims_by_cv,
    select_dims_by_zscore,
    select_top_k_low_variance,
)


def test_compute_signature_and_coverage():
    vectors = np.array([[0.0, 0.0], [2.0, 2.0], [4.0, 4.0]], dtype=np.float32)
    sig = compute_signature(vectors, n_requested=4)
    assert sig.n_genes == 3
    assert abs(sig.coverage - 0.75) < 1e-6
    assert np.allclose(sig.mean, np.array([2.0, 2.0], dtype=np.float32))
    assert np.allclose(sig.std, np.array([2.0, 2.0], dtype=np.float32))


def test_std_based_dimension_filters():
    std = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    mean = np.array([2.0, 2.0, 2.0], dtype=np.float32)

    mask_q = select_dims_by_quantile(std, q_low=0.5)
    assert mask_q.tolist() == [True, True, False]

    mask_cv = select_dims_by_cv(mean, std, cv_threshold=0.6)
    assert mask_cv.tolist() == [True, False, False]

    mask_z = select_dims_by_zscore(std, z_threshold=0.5)
    assert mask_z.tolist() == [True, False, False]

    mask_topk = select_top_k_low_variance(std, k=2)
    assert mask_topk.sum() == 2
    assert mask_topk[0] and mask_topk[1]


