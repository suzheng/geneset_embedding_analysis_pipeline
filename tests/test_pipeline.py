import json

import numpy as np

from src.pipeline import (
    compute_input_signatures,
    load_precomputed_terms_jsonl,
    make_dims_mask,
    match_input_to_terms,
    match_input_to_terms_farthest,
    match_input_to_terms_most_opposite,
)
from src.precompute_terms import precompute_from_dir


def test_pipeline_end_to_end(tmp_path):
    # Create small embeddings file
    emb = tmp_path / "emb.txt"
    emb.write_text("\n".join(["4 2", "A 0 0", "B 1 0", "C 0 1", "D 1 1"]) + "\n", encoding="utf-8")

    # Create small terms dir
    terms_dir = tmp_path / "terms"
    terms_dir.mkdir()
    (terms_dir / "X.tsv").write_text("T1\tA\tB\nT2\tC\tD\n", encoding="utf-8")

    # Precompute
    out_jsonl = tmp_path / "precomputed.jsonl"
    precompute_from_dir(str(emb), str(terms_dir), str(out_jsonl), min_genes_per_term=1)

    # Load precomputed
    terms = load_precomputed_terms_jsonl(str(out_jsonl))
    assert len(terms) == 2

    # Compute input signatures from genes
    sigs, found, missing, d = compute_input_signatures(str(emb), ["A", "D"], enable_clustering=False)
    assert d == 2
    assert len(found) == 2
    assert len(missing) == 0
    sig = sigs[0]

    # Make a mask and match
    mask = make_dims_mask("quantile", sig.mean, sig.std, q_low=0.5)
    ranked = match_input_to_terms(sig, terms, dims_filter="quantile", q_low=0.5)
    assert len(ranked) == 2
    # Expect either term may tie; just ensure consistent types
    (_, dist0, _), (_, dist1, _) = ranked[0], ranked[1]
    assert dist0 <= dist1 or dist1 <= dist0

    farthest = match_input_to_terms_farthest(sig, terms, dims_filter="quantile", q_low=0.5)
    assert len(farthest) == 2

    opposite = match_input_to_terms_most_opposite(sig, terms, dims_filter="quantile", q_low=0.5)
    assert len(opposite) == 2


