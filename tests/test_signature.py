import numpy as np

from src.signature import (
    compute_signature,
    compute_differential_signature,
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


def test_compute_differential_signature_basic():
    """Test basic differential signature computation."""
    # Input genes: mean [3, 3], reference genes: mean [1, 1]
    # Expected differential: [2, 2]
    input_vectors = np.array([[2.0, 2.0], [4.0, 4.0]], dtype=np.float32)
    reference_vectors = np.array([[0.0, 0.0], [2.0, 2.0]], dtype=np.float32)
    
    diff_sig = compute_differential_signature(input_vectors, reference_vectors, n_requested=2)
    
    # Check differential mean
    expected_diff_mean = np.array([2.0, 2.0], dtype=np.float32)
    assert np.allclose(diff_sig.mean, expected_diff_mean, atol=1e-5)
    
    # Check that pooled std is computed correctly
    # Input std: sqrt(((4-3)^2 + (2-3)^2)/1) = sqrt(2)
    # Reference std: sqrt(((2-1)^2 + (0-1)^2)/1) = sqrt(2)
    # Pooled std: sqrt((1*2 + 1*2)/(2+2-2)) = sqrt(4/2) = sqrt(2) ≈ 1.414
    expected_pooled_std = np.sqrt(2.0)
    assert np.allclose(diff_sig.std, expected_pooled_std, atol=1e-5)
    
    # Check metadata
    assert diff_sig.n_genes == 2
    assert diff_sig.coverage == 1.0


def test_compute_differential_signature_coverage():
    """Test that coverage is computed correctly in differential mode."""
    input_vectors = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]], dtype=np.float32)
    reference_vectors = np.array([[0.5, 0.5], [1.5, 1.5]], dtype=np.float32)
    
    # n_requested=5 but only 3 found
    diff_sig = compute_differential_signature(input_vectors, reference_vectors, n_requested=5)
    
    assert diff_sig.n_genes == 3
    assert abs(diff_sig.coverage - 0.6) < 1e-6  # 3/5 = 0.6


def test_compute_differential_signature_opposite_direction():
    """Test differential signature when reference has higher values."""
    # Input genes: mean [1, 1], reference genes: mean [3, 3]
    # Expected differential: [-2, -2]
    input_vectors = np.array([[0.0, 0.0], [2.0, 2.0]], dtype=np.float32)
    reference_vectors = np.array([[2.0, 2.0], [4.0, 4.0]], dtype=np.float32)
    
    diff_sig = compute_differential_signature(input_vectors, reference_vectors, n_requested=2)
    
    expected_diff_mean = np.array([-2.0, -2.0], dtype=np.float32)
    assert np.allclose(diff_sig.mean, expected_diff_mean, atol=1e-5)


def test_compute_differential_signature_multidimensional():
    """Test differential signature with more dimensions."""
    input_vectors = np.array([
        [1.0, 2.0, 3.0, 4.0],
        [2.0, 3.0, 4.0, 5.0],
        [3.0, 4.0, 5.0, 6.0]
    ], dtype=np.float32)
    
    reference_vectors = np.array([
        [0.5, 1.0, 1.5, 2.0],
        [1.5, 2.0, 2.5, 3.0]
    ], dtype=np.float32)
    
    diff_sig = compute_differential_signature(input_vectors, reference_vectors, n_requested=3)
    
    # Input mean: [2.0, 3.0, 4.0, 5.0]
    # Reference mean: [1.0, 1.5, 2.0, 2.5]
    # Differential: [1.0, 1.5, 2.0, 2.5]
    expected_diff_mean = np.array([1.0, 1.5, 2.0, 2.5], dtype=np.float32)
    assert np.allclose(diff_sig.mean, expected_diff_mean, atol=1e-5)
    
    assert diff_sig.mean.shape == (4,)
    assert diff_sig.std.shape == (4,)


def test_compute_differential_signature_empty_input():
    """Test that empty input vectors raise error."""
    input_vectors = np.array([], dtype=np.float32).reshape(0, 2)
    reference_vectors = np.array([[1.0, 1.0]], dtype=np.float32)
    
    try:
        compute_differential_signature(input_vectors, reference_vectors, n_requested=1)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "No input vectors" in str(e)


def test_compute_differential_signature_empty_reference():
    """Test that empty reference vectors raise error."""
    input_vectors = np.array([[1.0, 1.0]], dtype=np.float32)
    reference_vectors = np.array([], dtype=np.float32).reshape(0, 2)
    
    try:
        compute_differential_signature(input_vectors, reference_vectors, n_requested=1)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "No reference vectors" in str(e)


def test_compute_differential_signature_numerical_stability():
    """Test that std is clamped to epsilon for numerical stability."""
    # Create vectors with zero variance in some dimensions
    input_vectors = np.array([[1.0, 2.0], [1.0, 3.0]], dtype=np.float32)
    reference_vectors = np.array([[1.0, 1.0], [1.0, 2.0]], dtype=np.float32)
    
    diff_sig = compute_differential_signature(input_vectors, reference_vectors, n_requested=2, epsilon=1e-8)
    
    # First dimension has zero variance in both groups, should be clamped to epsilon
    assert diff_sig.std[0] >= 1e-8
    # Second dimension has non-zero variance
    assert diff_sig.std[1] > 1e-8



