import io

import numpy as np

from src.embeddings_io import get_vectors_for_genes, load_embeddings


def test_load_embeddings_and_lookup(tmp_path):
    p = tmp_path / "emb.txt"
    content = "\n".join(
        [
            "3 3",
            "GENE1 1 2 3",
            "GENE2 2 3 4",
            "GENE3 3 4 5",
            "",
        ]
    )
    p.write_text(content, encoding="utf-8")

    mapping, n, d = load_embeddings(str(p))
    assert d == 3
    assert len(mapping) == 3
    assert np.allclose(mapping["GENE1"], np.array([1, 2, 3], dtype=np.float32))

    vectors, found, missing = get_vectors_for_genes(mapping, ["gene1", "geneX", "GENE3"]) 
    assert vectors.shape == (2, 3)
    assert found == ["GENE1", "GENE3"]
    assert missing == ["GENEX"]
    assert np.allclose(vectors[0], np.array([1, 2, 3], dtype=np.float32))
    assert np.allclose(vectors[1], np.array([3, 4, 5], dtype=np.float32))


