from __future__ import annotations

import json
import os
from typing import List

from src.pipeline import compute_input_signatures, load_precomputed_terms_jsonl, match_input_to_terms
from src.precompute_terms import precompute_from_file


def main() -> None:
    # Paths
    project_root = os.path.dirname(os.path.dirname(__file__))
    embeddings_path = os.path.join(project_root, "data", "GeneRAIN-vec.32d.txt")
    clinvar_path = \
        "/Users/suzheng/noBackupData/BioBot/backend/data/tools/tool_data/database/gene_enrichment_analysis/ClinVar_2019.tsv"
    out_jsonl = os.path.join(project_root, "precomputed", "ClinVar_2019.signatures.jsonl")
    os.makedirs(os.path.dirname(out_jsonl), exist_ok=True)

    # Precompute ClinVar
    print("Precomputing ClinVar_2019 signatures...")
    precompute_from_file(embeddings_path, clinvar_path, out_jsonl, min_genes_per_term=3)
    print(f"Wrote: {out_jsonl}")

    # Example input gene list (pick a few common ones likely present)
    example_genes: List[str] = [
        "TP53", "EGFR", "BRCA1", "BRCA2", "PTEN", "KRAS", "APOE", "MYC", "ALK", "MTHFR", "CFTR", "TNF"
    ]

    # Compute input signature(s)
    sigs, found, missing, d = compute_input_signatures(
        embeddings_path,
        example_genes,
        enable_clustering=False,
    )
    print(f"Found {len(found)} genes in embeddings; missing: {len(missing)}")
    if missing:
        print("Missing symbols:", ", ".join(sorted(missing)))

    # Load precomputed ClinVar term signatures
    term_signatures = load_precomputed_terms_jsonl(out_jsonl)
    print(f"Loaded {len(term_signatures)} term signatures from ClinVar_2019")

    # Match and show top results for the first signature
    sig = sigs[0]
    ranked = match_input_to_terms(sig, term_signatures, dims_filter="quantile", q_low=0.3)
    print("Top 10 matches:")
    for i, (t, dist, cos) in enumerate(ranked[:10], start=1):
        print(f"{i:2d}. {t.library} | {t.term[:80]} | dist={dist:.4f} | cos={cos:.4f} | n={t.n_found} cov={t.coverage:.2f}")


if __name__ == "__main__":
    main()


