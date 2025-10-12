from __future__ import annotations

import json
import os
from typing import Any, Dict, Iterable, List, Tuple
import argparse
import textwrap

import numpy as np


def _ensure_dir(path: str) -> None:
    d = os.path.dirname(path)
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


def _replace_underscores(text: str) -> str:
    return (text or "").replace("_", " ")


def clean_library_name(library_name: str) -> str:
    """
    Clean library name for display in figures.
    
    - Replace underscores with spaces
    - Remove trailing _YYYY (year patterns like _2019, _2024)
    
    Args:
        library_name: Raw library name (e.g., "BioPlanet_2019" or "Reactome_Pathways_2024")
        
    Returns:
        Cleaned library name (e.g., "BioPlanet" or "Reactome Pathways")
    """
    import re
    # Remove trailing _YYYY pattern (year from 1900 to 2099)
    cleaned = re.sub(r'_\d{4}$', '', library_name)
    # Replace remaining underscores with spaces
    cleaned = cleaned.replace("_", " ")
    return cleaned


def wrap_title(title: str, max_chars_per_line: int = 60) -> str:
    """
    Wrap title text to prevent it from being too long.
    
    Args:
        title: The title text to wrap
        max_chars_per_line: Maximum characters per line
        
    Returns:
        Title with newlines inserted for wrapping
    """
    import textwrap
    lines = textwrap.wrap(title, width=max_chars_per_line)
    return "\n".join(lines)


def _configure_matplotlib() -> None:
    import matplotlib as mpl
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    mpl.rcParams["axes.titlesize"] = 8
    mpl.rcParams["axes.labelsize"] = 8
    mpl.rcParams["xtick.labelsize"] = 7.5
    mpl.rcParams["ytick.labelsize"] = 7.5
    mpl.rcParams["legend.fontsize"] = 7.5
    mpl.rcParams["figure.dpi"] = 300

    try:
        import seaborn as sns  # type: ignore
        palette = [
            "#0072b2",
            "#d55e00",
            "#009e73",
            "#cc79a7",
            "#56b4e9",
            "#e69f00",
            "#f0e442",
        ]
        # Remove gridlines globally per request
        sns.set_theme(style="white", context="paper")
        sns.set_palette(palette)
    except Exception:
        # Seaborn not available; matplotlib config above suffices
        pass


def _disable_grid(fig) -> None:
    try:
        for ax in fig.axes:
            ax.grid(False)
    except Exception:
        pass


def _wrap_or_truncate_labels(labels: List[str], max_chars_per_line: int = 26, max_lines: int = 2) -> List[str]:
    wrapped: List[str] = []
    for lab in labels:
        lab = _replace_underscores(lab)
        chunks = textwrap.wrap(lab, width=max_chars_per_line)
        if not chunks:
            wrapped.append("")
            continue
        if len(chunks) > max_lines:
            chunks = chunks[: max_lines]
            # add ellipsis to the last chunk if truncation happened
            if not lab.endswith(chunks[-1]):
                if len(chunks[-1]) > 3:
                    chunks[-1] = chunks[-1][:-3] + "..."
                else:
                    chunks[-1] = chunks[-1] + "..."
        wrapped.append("\n".join(chunks))
    return wrapped


def _project_root() -> str:
    return os.path.dirname(os.path.dirname(__file__))


def _load_example_signature(embeddings_path: str, example_genes: Iterable[str]):
    from src.embeddings_io import load_embeddings, get_vectors_for_genes
    from src.signature import compute_signature

    symbol_to_vector, _n, d = load_embeddings(embeddings_path)
    vectors, found, missing = get_vectors_for_genes(symbol_to_vector, example_genes)
    if vectors.shape[0] == 0:
        raise RuntimeError("No example genes found in embeddings.")
    sig = compute_signature(vectors, n_requested=len(list(example_genes)))
    return sig, found, missing, d


def _load_precomputed_terms(jsonl_path: str):
    from src.pipeline import load_precomputed_terms_jsonl

    return load_precomputed_terms_jsonl(jsonl_path)


def _save_figure(fig, figures_dir: str, base_filename: str, formats: List[str]) -> Dict[str, str]:
    paths: Dict[str, str] = {}
    for fmt in formats:
        out_path = os.path.join(figures_dir, f"{base_filename}.{fmt}")
        _ensure_dir(out_path)
        if fmt.lower() == "png":
            fig.savefig(out_path, format=fmt, dpi=300)
        else:
            fig.savefig(out_path, format=fmt)
        paths[fmt] = out_path
    return paths


def _figure_pipeline_overview(figures_dir: str, formats: List[str]) -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    width_in, height_in = 7.0, 3.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)
    ax.axis("off")

    steps = [
        "Preprocess",
        "QC",
        "Signature",
        "Term matching",
        "Cohort integration",
        "Group comparison",
    ]

    data_json_path = os.path.join(figures_dir, "pipeline_overview.data.json")
    _ensure_dir(data_json_path)
    with open(data_json_path, "w", encoding="utf-8") as f:
        json.dump({"nodes": steps, "edges": [[i, i + 1] for i in range(len(steps) - 1)]}, f, indent=2)

    n = len(steps)
    x0, x1 = 0.05, 0.95
    xs = np.linspace(x0, x1, n)
    y = 0.55
    box_w, box_h = 0.12, 0.28

    for i, (label, x) in enumerate(zip(steps, xs)):
        rect = patches.FancyBboxPatch(
            (x - box_w / 2.0, y - box_h / 2.0),
            box_w,
            box_h,
            boxstyle="round,pad=0.02,rounding_size=0.015",
            linewidth=1.0,
            facecolor="#56b4e9" if i % 2 else "#0072b2",
            edgecolor="#333333",
        )
        ax.add_patch(rect)
        ax.text(
            x,
            y,
            _replace_underscores(label),
            ha="center",
            va="center",
            color="white",
            fontsize=8,
            fontweight="bold",
        )
        if i < n - 1:
            x_next = xs[i + 1]
            ax.annotate(
                "",
                xy=(x_next - box_w / 2.5, y),
                xytext=(x + box_w / 2.0, y),
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#333333"),
            )

    saved = _save_figure(fig, figures_dir, "pipeline_overview", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_json_path


def _figure_term_signature_distributions(figures_dir: str, term_signatures, formats: List[str], library_name: str = "") -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    n_found = np.array([int(t.n_found) for t in term_signatures], dtype=np.int64)
    coverage = np.array([float(t.coverage) for t in term_signatures], dtype=np.float64)

    data_csv_path = os.path.join(figures_dir, "term_signatures_distributions.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("n_found,coverage\n")
        for nf, cov in zip(n_found.tolist(), coverage.tolist()):
            f.write(f"{nf},{cov}\n")

    width_in, height_in = 7.0, 3.4
    # Tweak right margin by providing wider figure and using tight layout
    fig, axes = plt.subplots(1, 2, figsize=(width_in + 0.6, height_in), constrained_layout=True)

    try:
        import seaborn as sns  # type: ignore
        sns.histplot(n_found, bins=30, ax=axes[0], color="#0072b2")
        sns.histplot(coverage, bins=30, ax=axes[1], color="#d55e00")
    except Exception:
        axes[0].hist(n_found, bins=30, color="#0072b2", edgecolor="#333333")
        axes[1].hist(coverage, bins=30, color="#d55e00", edgecolor="#333333")

    # Move "Term signature genes found" to y-axis, decrease title width
    if library_name:
        axes[0].set_title(wrap_title(_replace_underscores(library_name), max_chars_per_line=35))
    axes[0].set_xlabel(_replace_underscores("Genes found per term"))
    axes[0].set_ylabel(_replace_underscores("Term signature genes found"))

    if library_name:
        axes[1].set_title(wrap_title(_replace_underscores(library_name), max_chars_per_line=35))
    axes[1].set_xlabel(_replace_underscores("Coverage"))
    axes[1].set_ylabel(_replace_underscores("Term signature genes found"))
    axes[1].set_xlim(0.0, 1.0)

    # Panel labels
    axes[0].text(-0.12, 1.05, "A", transform=axes[0].transAxes, fontsize=9, fontweight="bold")
    axes[1].text(-0.12, 1.05, "B", transform=axes[1].transAxes, fontsize=9, fontweight="bold")

    for ax in axes:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(7.5)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "term_signatures_distributions", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_dimension_filtering(figures_dir: str, sig, formats: List[str], library_name: str = "") -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from src.pipeline import make_dims_mask

    std = np.asarray(sig.std, dtype=float)
    mean = np.asarray(sig.mean, dtype=float)
    d = std.size
    idx = np.arange(d)

    mask_quant = make_dims_mask("quantile", mean, std, q_low=0.3)
    threshold = float(np.quantile(std, 0.3))

    data_csv_path = os.path.join(figures_dir, "dimension_filtering.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("dim_index,std,mask_quantile\n")
        for i in range(d):
            f.write(f"{i},{float(std[i])},{int(bool(mask_quant[i]))}\n")

    width_in, height_in = 3.4, 3.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)

    order = np.argsort(std)
    std_sorted = std[order]
    ax.plot(np.arange(d), std_sorted, color="#0072b2", lw=1.5)
    ax.axhline(threshold, color="#d55e00", lw=1.0, ls="--", label=f"Quantile 0.3 thr = {threshold:.3g}")
    # Only show library name in title
    if library_name:
        ax.set_title(wrap_title(_replace_underscores(library_name), max_chars_per_line=40))
    ax.set_xlabel(_replace_underscores("Sorted dimension index"))
    ax.set_ylabel(_replace_underscores("Std"))
    ax.legend(frameon=False)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(7.5)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "dimension_filtering", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_dimension_filtering_bar(figures_dir: str, sig, formats: List[str], library_name: str = "") -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from src.pipeline import make_dims_mask

    std = np.asarray(sig.std, dtype=float)
    mean = np.asarray(sig.mean, dtype=float)
    d = std.size
    idx = np.arange(d)

    mask_quant = make_dims_mask("quantile", mean, std, q_low=0.3)

    data_csv_path = os.path.join(figures_dir, "dimension_filtering_bar.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("dim_index,std,mask_quantile\n")
        for i in range(d):
            f.write(f"{i},{float(std[i])},{int(bool(mask_quant[i]))}\n")

    # Two column width
    width_in, height_in = 7.0, 3.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)

    # Use third color (#009e73) for selected dimensions
    colors = np.array(["#56b4e9"] * d)
    colors[mask_quant] = "#d55e00"
    ax.bar(idx, std, color=colors, width=1.0)
    lib_suffix = f" ({library_name})" if library_name else ""
    ax.set_title(wrap_title(_replace_underscores(f"Dimension filtering (orange = selected){lib_suffix}")))
    ax.set_xlabel(_replace_underscores("Dimension index"))
    ax.set_ylabel(_replace_underscores("Std"))

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(7.5)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "dimension_filtering_bar", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_top_matches(figures_dir: str, sig, term_signatures, formats: List[str], library_name: str = "") -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from src.pipeline import match_input_to_terms

    ranked = match_input_to_terms(sig, term_signatures, dims_filter="quantile", q_low=0.3)
    top = ranked[:10]

    rows: List[Dict[str, Any]] = []
    for i, (t, dist, cos) in enumerate(top, start=1):
        rows.append(
            {
                "rank": i,
                "term_id": t.term_id,
                "library": t.library,
                "term": t.term,
                "distance": float(dist),
                "cosine": float(cos),
                "n_found": int(t.n_found),
                "coverage": float(t.coverage),
            }
        )

    data_csv_path = os.path.join(figures_dir, "top_matches_example.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("rank,term_id,library,term,distance,cosine,n_found,coverage\n")
        for r in rows:
            term_clean = r["term"].replace("\"", "'")
            f.write(
                f"{r['rank']},{r['term_id']},{r['library']},\"{term_clean}\",{r['distance']},{r['cosine']},{r['n_found']},{r['coverage']}\n"
            )

    width_in, height_in = 3.4, 4.5
    # Use tight_layout instead of constrained_layout to avoid warning with subplots_adjust
    fig, ax = plt.subplots(figsize=(width_in, height_in))

    # Wrap or truncate long labels to avoid overflow
    labels_raw = [r[0].term for r in top]
    labels = _wrap_or_truncate_labels(labels_raw, max_chars_per_line=26, max_lines=2)
    distances = [float(r[1]) for r in top]
    cosines = [float(r[2]) for r in top]
    y = np.arange(len(top))

    ax.barh(y, distances, color="#0072b2")
    ax.invert_yaxis()
    ax.set_xlabel(_replace_underscores("Variance weighted distance"))
    # Move "Top matches" to y-axis label
    ax.set_ylabel(_replace_underscores("Top matches"))
    ax.set_yticks(y)
    ax.set_yticklabels([_replace_underscores(l) for l in labels], fontsize=7.5)
    # Only keep library name in title, reduce title width
    if library_name:
        ax.set_title(wrap_title(_replace_underscores(library_name), max_chars_per_line=35))

    # Removed cosine annotations per request; values remain in data CSV

    # Dynamically expand left margin based on the longest label (approximate width)
    max_label_chars = max(len(l) for l in labels)
    # Heuristic: each character ~0.012 figure fraction at this font size
    extra_left = min(0.35, 0.012 * max_label_chars)
    plt.subplots_adjust(left=0.2 + extra_left)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "top_matches_example", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_similarity_distributions(figures_dir: str, sig, term_signatures, formats: List[str], library_name: str = "") -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns  # type: ignore
    except Exception:
        sns = None
    from src.match import variance_weighted_distance, cosine_similarity

    # Compute distance and cosine for all terms
    distances: List[float] = []
    cosines: List[float] = []
    for t in term_signatures:
        distances.append(variance_weighted_distance(sig.mean, sig.std, t.mean, t.std))
        cosines.append(cosine_similarity(sig.mean, t.mean))

    # Save input data
    data_csv_path = os.path.join(figures_dir, "similarity_distributions.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("distance,cosine\n")
        for dval, cval in zip(distances, cosines):
            f.write(f"{dval},{cval}\n")

    width_in, height_in = 7.0, 3.4
    fig, axes = plt.subplots(1, 2, figsize=(width_in, height_in), constrained_layout=True)
    
    if sns is not None:
        sns.histplot(distances, bins=50, ax=axes[0], color="#0072b2", kde=True, stat="density")
    else:
        axes[0].hist(distances, bins=50, color="#0072b2", edgecolor="#333333", density=True)
    # Remove "distribution" from title, decrease title width
    if library_name:
        axes[0].set_title(wrap_title(_replace_underscores(f"Distance ({library_name})"), max_chars_per_line=35))
    else:
        axes[0].set_title(_replace_underscores("Distance"))
    axes[0].set_xlabel(_replace_underscores("Variance weighted distance"))
    axes[0].set_ylabel(_replace_underscores("Density"))

    if sns is not None:
        sns.histplot(cosines, bins=50, ax=axes[1], color="#d55e00", kde=True, stat="density")
    else:
        axes[1].hist(cosines, bins=50, color="#d55e00", edgecolor="#333333", density=True)
    # Remove "distribution" from title, decrease title width
    if library_name:
        axes[1].set_title(wrap_title(_replace_underscores(f"Cosine similarity ({library_name})"), max_chars_per_line=35))
    else:
        axes[1].set_title(_replace_underscores("Cosine similarity"))
    axes[1].set_xlabel(_replace_underscores("Cosine similarity"))
    axes[1].set_ylabel(_replace_underscores("Density"))
    axes[1].set_xlim(-1.0, 1.0)

    # Panel labels
    axes[0].text(-0.12, 1.05, "A", transform=axes[0].transAxes, fontsize=9, fontweight="bold")
    axes[1].text(-0.12, 1.05, "B", transform=axes[1].transAxes, fontsize=9, fontweight="bold")

    for ax in axes:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(7.5)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "similarity_distributions", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_distance_cosine_hexbin(figures_dir: str, sig, term_signatures, formats: List[str], library_name: str = "", top_k: int = 15) -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from src.pipeline import match_input_to_terms
    from src.match import variance_weighted_distance, cosine_similarity

    # Compute scatter data
    distances: List[float] = []
    cosines: List[float] = []
    for t in term_signatures:
        distances.append(variance_weighted_distance(sig.mean, sig.std, t.mean, t.std))
        cosines.append(cosine_similarity(sig.mean, t.mean))

    # Save input data
    data_csv_path = os.path.join(figures_dir, "distance_cosine_hexbin.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("distance,cosine\n")
        for dval, cval in zip(distances, cosines):
            f.write(f"{dval},{cval}\n")

    # Top matches to highlight
    ranked = match_input_to_terms(sig, term_signatures, dims_filter="quantile", q_low=0.3)
    top = ranked[:top_k]

    width_in, height_in = 3.4, 3.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)
    hb = ax.hexbin(distances, cosines, gridsize=40, cmap="Blues", mincnt=1)
    cb = fig.colorbar(hb, ax=ax)
    cb.ax.set_ylabel(_replace_underscores("Density"))

    # Overlay top matches
    x_top = [float(r[1]) for r in top]
    y_top = [float(r[2]) for r in top]
    ax.scatter(
        x_top,
        y_top,
        c="#cc79a7",
        s=18,
        edgecolors="#333333",
        linewidths=0.5,
        zorder=3,
        label=_replace_underscores("Top matches"),
    )

    # Only show library name in title
    if library_name:
        ax.set_title(wrap_title(_replace_underscores(library_name), max_chars_per_line=40))
    ax.set_xlabel(_replace_underscores("Variance weighted distance"))
    ax.set_ylabel(_replace_underscores("Cosine similarity"))
    ax.set_ylim(-1.0, 1.0)
    ax.legend(loc="upper right", frameon=False)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "distance_cosine_hexbin", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_pca_term_space(figures_dir: str, sig, term_signatures, formats: List[str], library_name: str = "", top_k: int = 20) -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from src.pipeline import match_input_to_terms

    # Stack term means
    means = np.stack([np.asarray(t.mean, dtype=float) for t in term_signatures], axis=0)
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(means)
    
    # Project input signature into PCA space
    input_mean = np.asarray(sig.mean, dtype=float).reshape(1, -1)
    input_coord = pca.transform(input_mean)

    # Save input data
    data_csv_path = os.path.join(figures_dir, "pca_term_space.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("term_id,library,term,x,y,is_top,is_input\n")
        ranked = match_input_to_terms(sig, term_signatures, dims_filter="quantile", q_low=0.3)
        top_ids = set(id(r[0]) for r in ranked[:top_k])
        for (t, xy) in zip(term_signatures, coords):
            is_top = int(id(t) in top_ids)
            term_clean = t.term.replace("\"", "'")
            f.write(f"{t.term_id},{t.library},\"{term_clean}\",{xy[0]},{xy[1]},{is_top},0\n")
        # Add input signature
        f.write(f"INPUT,INPUT,\"Input genes\",{input_coord[0, 0]},{input_coord[0, 1]},0,1\n")

    width_in, height_in = 3.4, 3.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)

    # Plot background terms
    ax.scatter(coords[:, 0], coords[:, 1], s=6, c="#56b4e9", alpha=0.35, edgecolors="none")

    # Highlight top matches
    ranked = match_input_to_terms(sig, term_signatures, dims_filter="quantile", q_low=0.3)
    top = ranked[:top_k]
    top_indices = []
    # Map identity to index
    term_to_index = {id(t): i for i, t in enumerate(term_signatures)}
    for t, _d, _c in top:
        idx = term_to_index.get(id(t))
        if idx is not None:
            top_indices.append(idx)
    if top_indices:
        ax.scatter(
            coords[top_indices, 0],
            coords[top_indices, 1],
            s=16,
            c="#009e73",
            edgecolors="#333333",
            linewidths=0.5,
            label=_replace_underscores("Top matches"),
        )
    
    # Plot input genes embedding
    ax.scatter(
        input_coord[0, 0],
        input_coord[0, 1],
        s=80,
        c="#d55e00",
        marker="*",
        edgecolors="#333333",
        linewidths=0.8,
        zorder=5,
        label=_replace_underscores("Input genes"),
    )

    # Only show library name in title
    if library_name:
        ax.set_title(wrap_title(_replace_underscores(library_name), max_chars_per_line=40))
    
    # Add explained variance to axis labels
    var_pc1 = pca.explained_variance_ratio_[0] * 100
    var_pc2 = pca.explained_variance_ratio_[1] * 100
    ax.set_xlabel(f"PC1 ({var_pc1:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_pc2:.1f}%)")
    ax.legend(loc="best", frameon=False)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "pca_term_space", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_top_matches_opposite(figures_dir: str, sig, term_signatures, formats: List[str], library_name: str = "") -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from src.pipeline import match_input_to_terms_most_opposite

    ranked = match_input_to_terms_most_opposite(sig, term_signatures, dims_filter="quantile", q_low=0.3)
    top = ranked[:10]

    rows: List[Dict[str, Any]] = []
    for i, (t, dist, cos) in enumerate(top, start=1):
        rows.append(
            {
                "rank": i,
                "term_id": t.term_id,
                "library": t.library,
                "term": t.term,
                "distance": float(dist),
                "cosine": float(cos),
                "n_found": int(t.n_found),
                "coverage": float(t.coverage),
            }
        )

    data_csv_path = os.path.join(figures_dir, "top_opposite_matches_example.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("rank,term_id,library,term,distance,cosine,n_found,coverage\n")
        for r in rows:
            term_clean = r["term"].replace("\"", "'")
            f.write(
                f"{r['rank']},{r['term_id']},{r['library']},\"{term_clean}\",{r['distance']},{r['cosine']},{r['n_found']},{r['coverage']}\n"
            )

    width_in, height_in = 5.0, 4.5  # Increased width from 3.4 to 5.0
    fig, ax = plt.subplots(figsize=(width_in, height_in))

    labels_raw = [r[0].term for r in top]
    labels = _wrap_or_truncate_labels(labels_raw, max_chars_per_line=35, max_lines=2)  # Increased from 26 to 35
    # For opposite matches, we show cosine similarity (more negative = more opposite)
    cosines = [float(r[2]) for r in top]
    y = np.arange(len(top))

    ax.barh(y, cosines, color="#e69f00")  # Orange color for opposite
    ax.invert_yaxis()
    ax.set_xlabel(_replace_underscores("Cosine similarity (most negative)"))
    ax.set_ylabel(_replace_underscores("Top opposite matches"))
    ax.set_yticks(y)
    ax.set_yticklabels([_replace_underscores(l) for l in labels], fontsize=7.5)
    if library_name:
        ax.set_title(wrap_title(_replace_underscores(f"{library_name} (opposite)"), max_chars_per_line=45))

    max_label_chars = max(len(l) for l in labels)
    extra_left = min(0.35, 0.012 * max_label_chars)
    plt.subplots_adjust(left=0.2 + extra_left)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "top_opposite_matches_example", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_distance_cosine_hexbin_opposite(figures_dir: str, sig, term_signatures, formats: List[str], library_name: str = "", top_k: int = 15) -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from src.pipeline import match_input_to_terms_most_opposite
    from src.match import variance_weighted_distance, cosine_similarity

    # Compute scatter data
    distances: List[float] = []
    cosines: List[float] = []
    for t in term_signatures:
        distances.append(variance_weighted_distance(sig.mean, sig.std, t.mean, t.std))
        cosines.append(cosine_similarity(sig.mean, t.mean))

    # Save input data
    data_csv_path = os.path.join(figures_dir, "distance_cosine_hexbin_opposite.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("distance,cosine\n")
        for dval, cval in zip(distances, cosines):
            f.write(f"{dval},{cval}\n")

    # Top opposite matches to highlight
    ranked = match_input_to_terms_most_opposite(sig, term_signatures, dims_filter="quantile", q_low=0.3)
    top = ranked[:top_k]

    width_in, height_in = 3.4, 3.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)
    hb = ax.hexbin(distances, cosines, gridsize=40, cmap="Blues", mincnt=1)
    cb = fig.colorbar(hb, ax=ax)
    cb.ax.set_ylabel(_replace_underscores("Density"))

    # Overlay top opposite matches
    x_top = [float(r[1]) for r in top]
    y_top = [float(r[2]) for r in top]
    ax.scatter(
        x_top,
        y_top,
        c="#e69f00",  # Orange for opposite
        s=18,
        edgecolors="#333333",
        linewidths=0.5,
        zorder=3,
        label=_replace_underscores("Top opposite matches"),
    )

    if library_name:
        ax.set_title(wrap_title(_replace_underscores(f"{library_name} (opposite)"), max_chars_per_line=40))
    ax.set_xlabel(_replace_underscores("Variance weighted distance"))
    ax.set_ylabel(_replace_underscores("Cosine similarity"))
    ax.set_ylim(-1.0, 1.0)
    ax.legend(loc="upper right", frameon=False)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "distance_cosine_hexbin_opposite", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_pca_term_space_opposite(figures_dir: str, sig, term_signatures, formats: List[str], library_name: str = "", top_k: int = 20) -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from src.pipeline import match_input_to_terms_most_opposite

    # Stack term means
    means = np.stack([np.asarray(t.mean, dtype=float) for t in term_signatures], axis=0)
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(means)
    
    # Project input signature into PCA space
    input_mean = np.asarray(sig.mean, dtype=float).reshape(1, -1)
    input_coord = pca.transform(input_mean)

    # Save input data
    data_csv_path = os.path.join(figures_dir, "pca_term_space_opposite.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("term_id,library,term,x,y,is_top_opposite,is_input\n")
        ranked = match_input_to_terms_most_opposite(sig, term_signatures, dims_filter="quantile", q_low=0.3)
        top_ids = set(id(r[0]) for r in ranked[:top_k])
        for (t, xy) in zip(term_signatures, coords):
            is_top = int(id(t) in top_ids)
            term_clean = t.term.replace("\"", "'")
            f.write(f"{t.term_id},{t.library},\"{term_clean}\",{xy[0]},{xy[1]},{is_top},0\n")
        # Add input signature
        f.write(f"INPUT,INPUT,\"Input genes\",{input_coord[0, 0]},{input_coord[0, 1]},0,1\n")

    width_in, height_in = 3.4, 3.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)

    # Plot background terms
    ax.scatter(coords[:, 0], coords[:, 1], s=6, c="#56b4e9", alpha=0.35, edgecolors="none")

    # Highlight top opposite matches
    ranked = match_input_to_terms_most_opposite(sig, term_signatures, dims_filter="quantile", q_low=0.3)
    top = ranked[:top_k]
    top_indices = []
    term_to_index = {id(t): i for i, t in enumerate(term_signatures)}
    for t, _d, _c in top:
        idx = term_to_index.get(id(t))
        if idx is not None:
            top_indices.append(idx)
    if top_indices:
        ax.scatter(
            coords[top_indices, 0],
            coords[top_indices, 1],
            s=16,
            c="#e69f00",  # Orange for opposite
            edgecolors="#333333",
            linewidths=0.5,
            label=_replace_underscores("Top opposite matches"),
        )
    
    # Plot input genes embedding (same color as in normal PCA)
    ax.scatter(
        input_coord[0, 0],
        input_coord[0, 1],
        s=80,
        c="#d55e00",
        marker="*",
        edgecolors="#333333",
        linewidths=0.8,
        zorder=5,
        label=_replace_underscores("Input genes"),
    )

    if library_name:
        ax.set_title(wrap_title(_replace_underscores(f"{library_name} (opposite)"), max_chars_per_line=40))
    
    # Add explained variance to axis labels
    var_pc1 = pca.explained_variance_ratio_[0] * 100
    var_pc2 = pca.explained_variance_ratio_[1] * 100
    ax.set_xlabel(f"PC1 ({var_pc1:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_pc2:.1f}%)")
    ax.legend(loc="best", frameon=False)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "pca_term_space_opposite", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _write_figure_legend(figures_dir: str, library_name: str = "", has_opposite: bool = False) -> str:
    """Write comprehensive figure legends to FIGURE_LEGEND.md in the figures directory."""
    from datetime import datetime
    
    path = os.path.join(figures_dir, "FIGURE_LEGEND.md")
    
    # Generate comprehensive figure legends following academic publication standards
    legend_text = f"""# Figure Legends

This document provides comprehensive figure titles and legends for all generated figures, following academic publication standards.

## Figure 1: Gene Set Coverage Distribution in Reference Library

**Title:** Distribution of Gene Coverage Across Precomputed Term Signatures in {library_name or 'Reference Library'}

**Legend:** 
(A) Histogram showing the distribution of the number of genes found per term signature in the reference library. Each bar represents the frequency of terms with a given number of genes successfully mapped to the embedding space. This distribution reflects the size variability of gene sets in the reference library and their representation in the pre-trained embedding model.
(B) Histogram depicting the coverage distribution (fraction of genes with available embeddings) across all term signatures. Coverage values range from 0 to 1, where higher values indicate better representation of gene sets in the embedding space. High coverage (>0.8) suggests comprehensive gene annotation, while lower coverage may indicate specialized or poorly characterized gene sets.

**File:** `term_signatures_distributions.svg`

---

## Figure 2: Embedding Dimension Standard Deviation Profile

**Title:** Standard Deviation Profile Across Embedding Dimensions for Query Gene Set

**Legend:** 
Line plot showing the standard deviation (σ) of the query gene set signature across all embedding dimensions, sorted in ascending order. The horizontal dashed line indicates the quantile threshold (q = 0.3) used for dimension selection. Dimensions below this threshold exhibit low within-set variability and are retained for downstream matching, as they represent stable, coherent semantic features. Dimensions above the threshold show high variability and are filtered out to focus the analysis on the most informative features. This dimension selection strategy enhances the signal-to-noise ratio in similarity calculations.

**File:** `dimension_filtering.svg`

---

## Figure 3: Dimension Selection Summary

**Title:** Selected and Filtered Embedding Dimensions Based on Variability Threshold

**Legend:** 
Bar chart visualizing the standard deviation of each embedding dimension for the query gene set. Dimensions are colored by their selection status: orange bars represent selected dimensions (low variability, below the 0.3 quantile threshold) that are retained for matching, while blue bars indicate filtered dimensions (high variability) that are excluded from similarity calculations. This visualization provides a comprehensive view of which semantic features contribute to the functional interpretation of the gene set.

**File:** `dimension_filtering_bar.svg`

---

## Figure 4: Top Functionally Similar Terms

**Title:** Top Matching Gene Set Terms Ranked by Variance-Weighted Distance

**Legend:** 
Horizontal bar chart displaying the top 10 most similar terms from the reference library ({library_name or 'reference library'}), ranked by variance-weighted distance (lower values indicate stronger similarity). Each bar represents the distance metric, which combines the Euclidean distance between signature mean vectors with inverse-variance weighting to prioritize stable dimensions. Term names are displayed on the y-axis, ordered from most similar (top) to least similar (bottom). This ranking identifies the closest functional relatives of the query gene set in the reference library, providing biological context and potential functional annotations.

**File:** `top_matches_example.svg`

---

## Figure 5: Similarity Metric Distributions

**Title:** Distribution of Distance and Cosine Similarity Scores Across All Reference Terms

**Legend:** 
(A) Distribution of variance-weighted distance scores computed between the query gene set signature and all reference term signatures. Lower distance values indicate greater similarity. The distribution typically shows a right-skewed pattern, with most terms at moderate distances and a tail of highly similar matches. Kernel density estimation (KDE) overlay highlights the probability density.
(B) Distribution of cosine similarity scores between the query signature and all reference terms. Cosine similarity ranges from -1 (opposite direction) to +1 (same direction), measuring directional alignment in embedding space. Values near +1 indicate strong functional similarity, while values near 0 suggest orthogonal (unrelated) gene sets. These distributions provide a global view of how the query gene set relates to the entire reference library.

**File:** `similarity_distributions.svg`

---

## Figure 6: Joint Distance-Similarity Landscape

**Title:** Two-Dimensional Density Plot of Distance vs. Cosine Similarity for Reference Terms

**Legend:** 
Hexagonal binning plot showing the joint distribution of variance-weighted distance (x-axis) and cosine similarity (y-axis) across all reference term signatures. The color intensity (blue gradient) represents the density of terms in each hexagonal bin, revealing the relationship between these two complementary similarity metrics. Overlaid magenta points highlight the top matching terms, which typically cluster in the low-distance, high-cosine-similarity region (bottom-right quadrant). This visualization helps identify whether top matches are driven primarily by distance minimization, directional alignment, or both.

**File:** `distance_cosine_hexbin.svg`

---

## Figure 7: PCA Projection of Reference Term Space

**Title:** Principal Component Analysis of Reference Term Signatures with Top Matches Highlighted

**Legend:** 
Two-dimensional PCA projection of all reference term signatures (mean embedding vectors), with PC1 and PC2 axes labeled by their explained variance ratios. Background points (light blue) represent all terms in the reference library, illustrating the global structure of functional gene set space. Green points indicate the top matching terms, showing their spatial proximity to the query gene set. The query gene set signature is shown as an orange star, positioned within the reference term landscape. Terms that cluster near the query are functionally related, while distant terms represent distinct biological processes. This projection provides an intuitive geometric view of functional similarity relationships.

**File:** `pca_term_space.svg`

---

## Figure 8: Gene Set Signature Heatmap

**Title:** Heatmap Visualization of Query Gene Set Signature Mean and Standard Deviation Vectors

**Legend:** 
(A) Heatmap displaying the mean embedding vector (μ) of the query gene set across all dimensions. Color intensity (coolwarm colormap) represents the magnitude and sign of each dimension's mean value, with blue indicating negative values and red indicating positive values. This vector captures the centroid or "average position" of the gene set in embedding space, encoding its dominant functional characteristics.
(B) Heatmap showing the standard deviation vector (σ) across all dimensions, visualized using a magma colormap (dark purple to bright yellow). Higher values (yellow) indicate dimensions with high within-set gene variability, suggesting functional heterogeneity, while lower values (purple) represent coherent dimensions where genes are tightly clustered. Together, these two vectors define the complete signature used for functional term matching.

**File:** `signature_heatmap.svg`

---
"""

    # Add opposite term figures if enabled
    if has_opposite:
        legend_text += f"""## Figure 9: Top Functionally Opposite Terms

**Title:** Top Opposite Gene Set Terms Ranked by Negative Cosine Similarity

**Legend:** 
Horizontal bar chart displaying the top 10 most functionally opposite terms from the reference library ({library_name or 'reference library'}), ranked by cosine similarity (most negative values indicate strongest opposition). Each bar represents the cosine similarity score, where values approaching -1 indicate gene sets with opposite directionality in embedding space. This analysis identifies biological processes or gene sets that are functionally antagonistic or contrasting to the query gene set, which is particularly valuable for perturbation studies (e.g., drug vs. disease, knockdown vs. overexpression) where identifying opposing biological states is informative.

**File:** `top_opposite_matches_example.svg`

---

## Figure 10: Joint Distance-Similarity Landscape with Opposite Matches

**Title:** Two-Dimensional Density Plot Highlighting Top Opposite Gene Set Terms

**Legend:** 
Hexagonal binning plot showing the joint distribution of variance-weighted distance (x-axis) and cosine similarity (y-axis) across all reference term signatures. The color intensity (blue gradient) represents term density. Overlaid orange points highlight the top opposite terms, which typically cluster in regions with negative cosine similarity values (upper portion of plot), indicating directional opposition to the query gene set. This visualization reveals whether opposite matches are characterized by large distances, negative directional alignment, or both, providing insight into the nature of functional opposition.

**File:** `distance_cosine_hexbin_opposite.svg`

---

## Figure 11: PCA Projection with Opposite Matches Highlighted

**Title:** Principal Component Analysis of Reference Terms with Functionally Opposite Matches

**Legend:** 
Two-dimensional PCA projection of all reference term signatures, with PC1 and PC2 axes labeled by their explained variance ratios. Background points (light blue) represent all terms in the reference library. Orange points indicate the top opposite matching terms, showing their spatial positioning relative to the query gene set (orange star). Terms positioned on the opposite side of the PCA space from the query typically represent contrasting biological processes or phenotypes. This geometric visualization helps interpret the nature of functional opposition in terms of the major axes of variation in the reference library.

**File:** `pca_term_space_opposite.svg`

---
"""

    # Add footer with general notes
    legend_text += """## General Notes

**Figure Format:** All figures are provided in multiple formats (SVG for publications, PNG for presentations) to ensure broad compatibility.

**Data Files:** Each figure is accompanied by a corresponding data file (`.data.csv` or `.data.json`) containing the underlying numerical data used to generate the visualization, ensuring full reproducibility and enabling custom reanalysis.

**Color Schemes:** Figures use a consistent, colorblind-friendly palette optimized for scientific publication. Background grid lines are removed per publication standards, and font sizes are optimized for readability in both print and digital formats.

**Statistical Methods:** All distance calculations use variance-weighted metrics that prioritize stable dimensions (low within-set variability) and down-weight noisy dimensions (high variability). Cosine similarity provides a complementary directional measure. Dimension selection (quantile filtering at q=0.3) focuses the analysis on the most informative semantic features.

**Reproducibility:** All figures are generated using standardized matplotlib configuration with embedded fonts (Type 42 PDF fonts, SVG none), ensuring consistent appearance across platforms and publications. Random seeds are fixed for all analyses involving stochastic components (e.g., PCA initialization).

---

**Generated by:** Gene Set Embedding Analysis Pipeline  
**Document Version:** 1.0  
**Last Updated:** {datetime.utcnow().strftime('%Y-%m-%d')}
"""

    _ensure_dir(path)
    with open(path, "w", encoding="utf-8") as f:
        f.write(legend_text)
    
    return path


def _figure_signature_heatmap(figures_dir: str, sig, formats: List[str], library_name: str = "") -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from src.pipeline import make_dims_mask

    mean = np.asarray(sig.mean, dtype=float)
    std = np.asarray(sig.std, dtype=float)
    d = mean.size
    mask = make_dims_mask("quantile", mean, std, q_low=0.3)

    # Save input data
    data_json_path = os.path.join(figures_dir, "signature_heatmap.data.json")
    _ensure_dir(data_json_path)
    with open(data_json_path, "w", encoding="utf-8") as f:
        json.dump({"mean": mean.tolist(), "std": std.tolist(), "mask_quantile": [bool(x) for x in mask.tolist()]}, f, indent=2)

    width_in, height_in = 7.0, 2.8
    # Wider figure and explicit margins to prevent colorbar label cutoff
    fig, axes = plt.subplots(1, 2, figsize=(width_in + 1.0, height_in))
    
    # Heatmap of mean (as single-row image)
    im0 = axes[0].imshow(mean[np.newaxis, :], aspect="auto", cmap="coolwarm", extent=[1, d+1, 0, 1])
    # Remove "by dimension" from title, decrease title width
    if library_name:
        axes[0].set_title(wrap_title(_replace_underscores(f"Signature mean ({library_name})"), max_chars_per_line=35))
    else:
        axes[0].set_title(_replace_underscores("Signature mean"))
    axes[0].set_yticks([])
    axes[0].set_xlabel(_replace_underscores("Dimension index"))
    divider0 = make_axes_locatable(axes[0])
    cax0 = divider0.append_axes("right", size="2.5%", pad=0.10)
    fig.colorbar(im0, cax=cax0)

    # Heatmap of std
    im1 = axes[1].imshow(std[np.newaxis, :], aspect="auto", cmap="magma", extent=[1, d+1, 0, 1])
    # Remove "by dimension" from title, decrease title width
    if library_name:
        axes[1].set_title(wrap_title(_replace_underscores(f"Signature std ({library_name})"), max_chars_per_line=35))
    else:
        axes[1].set_title(_replace_underscores("Signature std"))
    axes[1].set_yticks([])
    axes[1].set_xlabel(_replace_underscores("Dimension index"))
    divider1 = make_axes_locatable(axes[1])
    cax1 = divider1.append_axes("right", size="2.5%", pad=0.10)
    fig.colorbar(im1, cax=cax1)

    # Panel labels
    axes[0].text(-0.08, 1.05, "A", transform=axes[0].transAxes, fontsize=9, fontweight="bold")
    axes[1].text(-0.08, 1.05, "B", transform=axes[1].transAxes, fontsize=9, fontweight="bold")

    # Adjust margins - increase top margin to 0.85 for more space, prevent text cutoff on right
    plt.subplots_adjust(left=0.08, right=0.90, top=0.85, bottom=0.18)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "signature_heatmap", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_json_path


def main(argv: List[str] | None = None) -> None:
    _configure_matplotlib()

    project_root = _project_root()
    figures_dir = os.path.join(project_root, "figures")
    embeddings_path = os.path.join(project_root, "data", "GeneRAIN-vec.200d.txt")
    precomputed_jsonl = os.path.join(project_root, "precomputed", "ClinVar_2019.signatures.jsonl")

    parser = argparse.ArgumentParser(description="Generate publication-quality figures for the pipeline")
    parser.add_argument("--also-png", action="store_true", help="Also save each figure as PNG (in addition to SVG)")
    args = parser.parse_args(argv)
    formats: List[str] = ["svg", "png"] if args.also_png else ["svg"]

    if not os.path.exists(embeddings_path):
        raise FileNotFoundError(f"Embeddings file not found: {embeddings_path}")
    if not os.path.exists(precomputed_jsonl):
        raise FileNotFoundError(f"Precomputed term signatures not found: {precomputed_jsonl}")

    example_genes: List[str] = [
        "TP53", "EGFR", "BRCA1", "BRCA2", "PTEN", "KRAS", "APOE", "MYC", "ALK", "MTHFR", "CFTR", "TNF"
    ]

    sig, found, missing, d = _load_example_signature(embeddings_path, example_genes)
    term_signatures = _load_precomputed_terms(precomputed_jsonl)

    metadata_path = os.path.join(figures_dir, "figure_generation.metadata.json")
    _ensure_dir(metadata_path)
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "embeddings_path": os.path.abspath(embeddings_path),
                "precomputed_jsonl": os.path.abspath(precomputed_jsonl),
                "example_genes": example_genes,
                "n_genes_found": len(found),
                "n_genes_missing": len(missing),
                "embedding_dim": int(d),
                "notes": "Figures generated for scientific publication; fonts and sizes configured per spec.",
            },
            f,
            indent=2,
        )

    # Extract library name from precomputed_jsonl
    library_name = ""
    if precomputed_jsonl:
        library_name = os.path.basename(precomputed_jsonl).replace('.signatures.jsonl', '')
        library_name = clean_library_name(library_name)
    
    outputs: List[Tuple[str, str]] = []
    # Pipeline overview removed per user request
    outputs.append(_figure_term_signature_distributions(figures_dir, term_signatures, formats, library_name))
    outputs.append(_figure_dimension_filtering(figures_dir, sig, formats, library_name))
    outputs.append(_figure_dimension_filtering_bar(figures_dir, sig, formats, library_name))
    outputs.append(_figure_top_matches(figures_dir, sig, term_signatures, formats, library_name))
    # Advanced/new figures
    outputs.append(_figure_similarity_distributions(figures_dir, sig, term_signatures, formats, library_name))
    outputs.append(_figure_distance_cosine_hexbin(figures_dir, sig, term_signatures, formats, library_name))
    outputs.append(_figure_pca_term_space(figures_dir, sig, term_signatures, formats, library_name))
    outputs.append(_figure_signature_heatmap(figures_dir, sig, formats, library_name))

    # Generate figure legend documentation
    legend_path = _write_figure_legend(figures_dir, library_name, has_opposite=False)
    print(f" - Figure legends: {os.path.abspath(legend_path)}")

    index_path = os.path.join(figures_dir, "index.json")
    _ensure_dir(index_path)
    with open(index_path, "w", encoding="utf-8") as f:
        entries: List[Dict[str, Any]] = []
        for fig, data in outputs:
            entry: Dict[str, Any] = {"figure_svg": os.path.abspath(fig), "data": os.path.abspath(data)}
            if "png" in formats:
                entry["figure_png"] = os.path.abspath(os.path.splitext(fig)[0] + ".png")
            entries.append(entry)
        json.dump({"figures": entries}, f, indent=2)

    print("Generated figures:")
    for fig, data in outputs:
        print(f" - {os.path.abspath(fig)}\n   data: {os.path.abspath(data)}")
        if "png" in formats:
            print(f"   png:  {os.path.abspath(os.path.splitext(fig)[0] + '.png')}")


if __name__ == "__main__":
    main()


