"""Integration tests for differential signature mode in the pipeline."""

import tempfile
from pathlib import Path

import numpy as np

from src.pipeline import compute_input_signatures


def test_compute_input_signatures_standard_mode(tmp_path):
    """Test standard mode without reference genes."""
    # Create test embeddings file
    emb_file = tmp_path / "embeddings.txt"
    content = "\n".join([
        "5 3",
        "GENE1 1.0 2.0 3.0",
        "GENE2 2.0 3.0 4.0",
        "GENE3 3.0 4.0 5.0",
        "GENE4 4.0 5.0 6.0",
        "GENE5 5.0 6.0 7.0",
        "",
    ])
    emb_file.write_text(content, encoding="utf-8")
    
    input_genes = ["GENE1", "GENE2", "GENE3"]
    
    signatures, found, missing, d = compute_input_signatures(
        embeddings_path=str(emb_file),
        input_genes=input_genes,
        enable_clustering=False,
        cluster_min_genes=20,
        reference_genes=None,
    )
    
    assert len(signatures) == 1
    assert len(found) == 3
    assert len(missing) == 0
    assert d == 3
    
    # Check that the mean is the average of the three genes
    # (1+2+3)/3 = 2, (2+3+4)/3 = 3, (3+4+5)/3 = 4
    expected_mean = np.array([2.0, 3.0, 4.0], dtype=np.float32)
    assert np.allclose(signatures[0].mean, expected_mean, atol=1e-5)


def test_compute_input_signatures_differential_mode(tmp_path):
    """Test differential mode with reference genes."""
    # Create test embeddings file
    emb_file = tmp_path / "embeddings.txt"
    content = "\n".join([
        "5 3",
        "GENE1 1.0 2.0 3.0",
        "GENE2 2.0 3.0 4.0",
        "GENE3 3.0 4.0 5.0",
        "REF1 0.0 1.0 2.0",
        "REF2 1.0 2.0 3.0",
        "",
    ])
    emb_file.write_text(content, encoding="utf-8")
    
    input_genes = ["GENE1", "GENE2", "GENE3"]
    reference_genes = ["REF1", "REF2"]
    
    signatures, found, missing, d = compute_input_signatures(
        embeddings_path=str(emb_file),
        input_genes=input_genes,
        enable_clustering=True,  # Should be ignored in differential mode
        cluster_min_genes=2,
        reference_genes=reference_genes,
    )
    
    # In differential mode, should return only one signature (no clustering)
    assert len(signatures) == 1
    assert len(found) == 3
    assert len(missing) == 0
    assert d == 3
    
    # Input mean: (1+2+3)/3 = 2, (2+3+4)/3 = 3, (3+4+5)/3 = 4
    # Reference mean: (0+1)/2 = 0.5, (1+2)/2 = 1.5, (2+3)/2 = 2.5
    # Differential: 2-0.5 = 1.5, 3-1.5 = 1.5, 4-2.5 = 1.5
    expected_diff_mean = np.array([1.5, 1.5, 1.5], dtype=np.float32)
    assert np.allclose(signatures[0].mean, expected_diff_mean, atol=1e-5)


def test_compute_input_signatures_differential_mode_missing_genes(tmp_path):
    """Test differential mode with some missing genes."""
    # Create test embeddings file
    emb_file = tmp_path / "embeddings.txt"
    content = "\n".join([
        "3 2",
        "GENE1 1.0 2.0",
        "GENE2 2.0 3.0",
        "REF1 0.5 1.0",
        "",
    ])
    emb_file.write_text(content, encoding="utf-8")
    
    input_genes = ["GENE1", "GENE2", "MISSING_GENE"]
    reference_genes = ["REF1"]
    
    signatures, found, missing, d = compute_input_signatures(
        embeddings_path=str(emb_file),
        input_genes=input_genes,
        enable_clustering=False,
        cluster_min_genes=20,
        reference_genes=reference_genes,
    )
    
    assert len(signatures) == 1
    assert len(found) == 2
    assert len(missing) == 1
    assert "MISSING_GENE" in missing
    assert d == 2
    
    # Input mean: (1+2)/2 = 1.5, (2+3)/2 = 2.5
    # Reference mean: 0.5, 1.0
    # Differential: 1.5-0.5 = 1.0, 2.5-1.0 = 1.5
    expected_diff_mean = np.array([1.0, 1.5], dtype=np.float32)
    assert np.allclose(signatures[0].mean, expected_diff_mean, atol=1e-5)
    
    # Check coverage is based on input genes
    assert signatures[0].n_genes == 2
    assert abs(signatures[0].coverage - 2.0/3.0) < 1e-6


def test_compute_input_signatures_differential_mode_no_clustering():
    """Verify that clustering is skipped in differential mode even when enabled."""
    # This test verifies the conceptual behavior by checking that only one signature is returned
    with tempfile.TemporaryDirectory() as tmpdir:
        emb_file = Path(tmpdir) / "embeddings.txt"
        
        # Create many genes to potentially trigger clustering
        lines = ["10 3"]
        for i in range(10):
            lines.append(f"GENE{i} {i}.0 {i+1}.0 {i+2}.0")
        for i in range(5):
            lines.append(f"REF{i} {i*0.5} {i*0.5+0.5} {i*0.5+1.0}")
        lines.append("")
        
        emb_file.write_text("\n".join(lines), encoding="utf-8")
        
        input_genes = [f"GENE{i}" for i in range(10)]
        reference_genes = [f"REF{i}" for i in range(5)]
        
        signatures, found, missing, d = compute_input_signatures(
            embeddings_path=str(emb_file),
            input_genes=input_genes,
            enable_clustering=True,  # Enabled but should be ignored
            cluster_min_genes=5,  # Low enough to allow clustering
            reference_genes=reference_genes,
        )
        
        # Should return only 1 signature (no clustering in differential mode)
        assert len(signatures) == 1
        assert len(found) == 10
        assert signatures[0].mean.shape == (3,)


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])

