from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import List, Optional, Tuple


def _read_text_if_exists(file_path: Path) -> str:
    if not file_path.exists():
        return ""
    try:
        return file_path.read_text(encoding="utf-8")
    except Exception:
        return ""


def _collect_shared_docs(entity_dir: Path) -> Tuple[str, str]:
    """
    Collect background and shared method docs, and figure legend text.

    Returns a tuple of (background_methods_text, figure_legend_text).
    """
    shared_dir = entity_dir / "shared_doc"
    parts: List[str] = []

    # Background
    parts.append("===== BACKGROUND =====\n")
    parts.append(_read_text_if_exists(shared_dir / "BACKGROUND.md"))

    # Methods (shared)
    parts.append("\n===== METHODS =====\n")
    for fname in [
        "SHARED_qc.METHOD.md",
        "SHARED_cohort_integration.METHOD.md",
        "SHARED_signature.METHOD.md",
        "SHARED_term_matching.METHOD.md",
    ]:
        text = _read_text_if_exists(shared_dir / fname)
        if text:
            parts.append(f"\n---\n# {fname}\n")
            parts.append(text)

    methods_text = "\n".join([p for p in parts if p])
    figure_legend_text = _read_text_if_exists(shared_dir / "FIGURE_LEGEND.md")
    return methods_text, figure_legend_text


def _find_first_term_matches_csv(lib_dir: Path) -> Optional[Path]:
    candidates = sorted(lib_dir.glob("*_term_matches.csv"))
    return candidates[0] if candidates else None


def _read_top_matches_term_distance(csv_path: Path, max_rows: int = 5) -> List[Tuple[str, str]]:
    """
    Read the CSV produced by term matching and extract top rows for columns 'term' and 'distance'.
    Emulates: csv2t | cut -f3-4 | head -6 (header + top 5 rows).
    """
    top_rows: List[Tuple[str, str]] = []
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, [])
        if not header:
            return top_rows
        # Map columns
        try:
            term_idx = header.index("term")
        except ValueError:
            term_idx = None
        try:
            distance_idx = header.index("distance")
        except ValueError:
            distance_idx = None
        if term_idx is None or distance_idx is None:
            return top_rows
        for i, row in enumerate(reader):
            if i >= max_rows:
                break
            if len(row) <= max(term_idx, distance_idx):
                continue
            term = row[term_idx]
            distance = row[distance_idx]
            top_rows.append((term, distance))
    return top_rows


def _read_llm_interpretation(lib_dir: Path) -> str:
    candidates = sorted(lib_dir.glob("*_llm_interpretation.txt"))
    if not candidates:
        return ""
    try:
        return candidates[0].read_text(encoding="utf-8")
    except Exception:
        return ""


def _list_svg_figures(lib_dir: Path) -> List[str]:
    fig_dir = lib_dir / "figures"
    if not fig_dir.exists():
        return []
    return sorted([p.name for p in fig_dir.glob("*.svg")])


def generate_manuscript_prompt(entity_dir: str, output_filename: Optional[str] = None) -> str:
    """
    Generate a consolidated manuscript-writing prompt from shared docs and per-library outputs.

    Writes a text file under the entity directory and returns its path.
    """
    entity_path = Path(entity_dir).resolve()
    if not entity_path.exists():
        raise FileNotFoundError(f"Entity directory not found: {entity_dir}")

    entity_id = entity_path.name

    # Collect shared background/methods and figure legend
    background_methods_text, figure_legend_text = _collect_shared_docs(entity_path)

    # Compose prompt header
    lines: List[str] = []
    lines.append(
        "Please read the background, method, results and the interpretation below, and write a manuscript, "
        "keep references as [URL of reference], below is the related information:\n\n"
    )

    # Add Background + Methods
    if background_methods_text:
        lines.append(background_methods_text)
        lines.append("\n\n")

    # Results & Interpretation per library
    lines.append("===== RESULTS AND INTERPRETATION =====\n")

    # Enumerate library subfolders (exclude shared_doc)
    subdirs = [p for p in entity_path.iterdir() if p.is_dir() and p.name != "shared_doc"]
    subdirs.sort(key=lambda p: p.name)

    for lib_dir in subdirs:
        # Skip directories that do not look like per-library outputs
        has_tm = any(lib_dir.glob("*_term_matches.csv"))
        has_fig = (lib_dir / "figures").exists()
        has_any = has_tm or has_fig or any(lib_dir.glob("*_llm_interpretation.txt"))
        if not has_any:
            continue

        library_name = lib_dir.name.replace(f"{entity_id}_", "")
        lines.append(f"\n---\n### Library: {library_name}\n")

        # Top matches (term, distance)
        tm_csv = _find_first_term_matches_csv(lib_dir)
        if tm_csv:
            lines.append("Top Matches (term\tdistance):\n")
            lines.append("term\tdistance")
            for term, distance in _read_top_matches_term_distance(tm_csv, max_rows=5):
                lines.append(f"{term}\t{distance}")
            lines.append("")

        # LLM interpretation
        interp = _read_llm_interpretation(lib_dir)
        if interp:
            lines.append("Interpretation:\n")
            lines.append(interp)
            lines.append("")

        # Figures
        svg_files = _list_svg_figures(lib_dir)
        if svg_files:
            lines.append("Figures available (SVG):\n")
            for fn in svg_files:
                lines.append(f"- {fn}")
            lines.append("")

    # General figure guidance (canonical list) and legends
    lines.append(
        "Keep in mind, each library has the following canonical figures (if generated):\n"
        "figures/dimension_filtering.svg\n"
        "figures/dimension_filtering_bar.svg\n"
        "figures/distance_cosine_hexbin.svg\n"
        "figures/pca_term_space.svg\n"
        "figures/signature_heatmap.svg\n"
        "figures/similarity_distributions.svg\n"
        "figures/term_signatures_distributions.svg\n"
        "figures/top_matches_example.svg\n\n"
        "Please add appropriate figures for the significant results to appropriate positions in the manuscript.\n\n"
    )

    if figure_legend_text:
        lines.append("===== FIGURE LEGENDS =====\n")
        lines.append(figure_legend_text)
        lines.append("\n")

    # Write output file
    if output_filename:
        out_path = entity_path / output_filename
    else:
        out_path = entity_path / f"{entity_id}_manuscript_prompt.txt"

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return str(out_path)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate manuscript prompt from analysis outputs")
    p.add_argument("entity_dir", help="Path to the entity output directory (contains shared_doc and per-library folders)")
    p.add_argument("--output-file", default=None, help="Name of the output text file (default: <entity_id>_manuscript_prompt.txt)")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out_path = generate_manuscript_prompt(args.entity_dir, args.output_file)
    print(out_path)


if __name__ == "__main__":
    main()









