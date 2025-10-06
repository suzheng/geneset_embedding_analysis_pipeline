from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from .embeddings_io import get_vectors_for_genes, load_embeddings
from .signature import compute_signature
from .terms_io import iter_terms_from_tsv, list_tsv_files


@dataclass
class StoredTermSignature:
    term_id: str
    library: str
    term: str
    d: int
    mean: List[float]
    std: List[float]
    n_found: int
    coverage: float


def _term_id(library: str, term: str) -> str:
    key = f"{library}::{term}".encode("utf-8")
    return hashlib.md5(key).hexdigest()


def precompute_from_dir(
    embeddings_path: str,
    terms_dir: str,
    output_jsonl_path: str,
    min_genes_per_term: int = 5,
) -> None:
    symbol_to_vector, _n, d = load_embeddings(embeddings_path)

    out_dir = os.path.dirname(output_jsonl_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    with open(output_jsonl_path, "w", encoding="utf-8") as out:
        for tsv_path in list_tsv_files(terms_dir):
            library = os.path.basename(tsv_path)
            for term, genes in iter_terms_from_tsv(tsv_path):
                vectors, found, missing = get_vectors_for_genes(symbol_to_vector, genes)
                if vectors.shape[0] < min_genes_per_term:
                    continue
                sig = compute_signature(vectors, n_requested=len(genes))
                record = StoredTermSignature(
                    term_id=_term_id(library, term),
                    library=library,
                    term=term,
                    d=d,
                    mean=sig.mean.astype(float).tolist(),
                    std=sig.std.astype(float).tolist(),
                    n_found=sig.n_genes,
                    coverage=float(sig.coverage),
                )
                out.write(json.dumps(asdict(record)) + "\n")


def precompute_from_file(
    embeddings_path: str,
    tsv_file_path: str,
    output_jsonl_path: str,
    min_genes_per_term: int = 5,
) -> None:
    symbol_to_vector, _n, d = load_embeddings(embeddings_path)

    out_dir = os.path.dirname(output_jsonl_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    library = os.path.basename(tsv_file_path)
    with open(output_jsonl_path, "w", encoding="utf-8") as out:
        for term, genes in iter_terms_from_tsv(tsv_file_path):
            vectors, found, missing = get_vectors_for_genes(symbol_to_vector, genes)
            if vectors.shape[0] < min_genes_per_term:
                continue
            sig = compute_signature(vectors, n_requested=len(genes))
            record = StoredTermSignature(
                term_id=_term_id(library, term),
                library=library,
                term=term,
                d=d,
                mean=sig.mean.astype(float).tolist(),
                std=sig.std.astype(float).tolist(),
                n_found=sig.n_genes,
                coverage=float(sig.coverage),
            )
            out.write(json.dumps(asdict(record)) + "\n")


