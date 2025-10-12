from __future__ import annotations

import io
from typing import Dict, Iterable, List, Tuple

import numpy as np


def _normalize_symbol(symbol: str) -> str:
    return symbol.strip().upper()


def load_embeddings(embeddings_path: str) -> Tuple[Dict[str, np.ndarray], int, int]:
    """
    Load gene embeddings from a whitespace-delimited text file.

    The first line should be "<num_records> <num_dimensions>". Subsequent lines are
    "<GENE_SYMBOL> <d1> <d2> ... <dD>".

    Returns a tuple (symbol_to_vector, num_records, num_dimensions).
    """
    symbol_to_vector: Dict[str, np.ndarray] = {}

    with io.open(embeddings_path, "r", encoding="utf-8") as f:
        header = f.readline().strip().split()
        if len(header) < 2:
            raise ValueError(
                f"Embeddings header must contain two integers: count and dimension. Got: {header}"
            )
        try:
            _n = int(header[0])
            d = int(header[1])
        except ValueError as exc:
            raise ValueError(
                f"Embeddings header must contain integers. Got: {header}"
            ) from exc

        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            symbol = _normalize_symbol(parts[0])
            values = parts[1:]
            if len(values) != d:
                # Skip malformed lines silently to be robust.
                continue
            try:
                vec = np.array([float(x) for x in values], dtype=np.float32)
            except ValueError:
                continue
            symbol_to_vector[symbol] = vec

    return symbol_to_vector, len(symbol_to_vector), d


def get_vectors_for_genes(
    symbol_to_vector: Dict[str, np.ndarray], genes: Iterable[str]
) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Given a symbol->vector dictionary and an iterable of gene symbols, return:
      - stacked vectors (shape: [n_found, D])
      - found symbols (normalized)
      - missing symbols (normalized)
    """
    found_vectors: List[np.ndarray] = []
    found_symbols: List[str] = []
    missing_symbols: List[str] = []

    for gene in genes:
        symbol = _normalize_symbol(gene)
        vec = symbol_to_vector.get(symbol)
        if vec is None:
            missing_symbols.append(symbol)
        else:
            found_symbols.append(symbol)
            found_vectors.append(vec)

    if found_vectors:
        stacked = np.stack(found_vectors, axis=0)
    else:
        stacked = np.empty((0, 0), dtype=np.float32)

    return stacked, found_symbols, missing_symbols



