from __future__ import annotations

import csv
import io
import os
from typing import Generator, Iterable, List, Tuple


def _normalize_symbol(symbol: str) -> str:
    return symbol.strip().upper()


def iter_terms_from_tsv(file_path: str) -> Generator[Tuple[str, List[str]], None, None]:
    """
    Yield (term, genes) from an Enrichr-style TSV where the first column is the term
    and subsequent columns are gene symbols. Empty cells are ignored.
    """
    def _parse_gene_cell(cell: str) -> str:
        # Accept forms like "GENE" or "GENE,VALUE"; ignore the numeric part if present
        token = (cell or "").strip()
        if "," in token:
            token = token.split(",", 1)[0]
        return _normalize_symbol(token)

    with io.open(file_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            term = row[0].strip()
            raw_genes = [c for c in row[1:] if c and c.strip()]
            genes = []
            for c in raw_genes:
                g = _parse_gene_cell(c)
                if g:
                    genes.append(g)
            if term:
                yield term, genes


def list_tsv_files(terms_dir: str) -> List[str]:
    paths: List[str] = []
    for name in sorted(os.listdir(terms_dir)):
        if name.lower().endswith(".tsv"):
            paths.append(os.path.join(terms_dir, name))
    return paths


