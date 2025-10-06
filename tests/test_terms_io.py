import os

from src.terms_io import iter_terms_from_tsv, list_tsv_files


def test_iter_terms_from_tsv(tmp_path):
    p = tmp_path / "terms.tsv"
    p.write_text(
        "TERM1\tGENE1\tGENE2\nTERM2\tGENE3,-1.23\t\t\n\n",
        encoding="utf-8",
    )
    rows = list(iter_terms_from_tsv(str(p)))
    assert rows == [("TERM1", ["GENE1", "GENE2"]), ("TERM2", ["GENE3"])]


def test_list_tsv_files(tmp_path):
    a = tmp_path / "A.tsv"
    b = tmp_path / "B.txt"
    c = tmp_path / "C.TSV"
    a.write_text("", encoding="utf-8")
    b.write_text("", encoding="utf-8")
    c.write_text("", encoding="utf-8")
    listed = list_tsv_files(str(tmp_path))
    # Sorted order, includes .tsv and .TSV
    assert os.path.basename(listed[0]) in {"A.tsv", "C.TSV"}
    assert os.path.basename(listed[-1]) in {"A.tsv", "C.TSV"}
    assert len(listed) == 2


