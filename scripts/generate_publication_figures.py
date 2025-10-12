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


def _figure_term_signature_distributions(figures_dir: str, term_signatures, formats: List[str]) -> Tuple[str, str]:
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

    axes[0].set_title(_replace_underscores("Term signature genes found"))
    axes[0].set_xlabel(_replace_underscores("Genes found per term"))
    axes[0].set_ylabel(_replace_underscores("Count"))

    axes[1].set_title(_replace_underscores("Term signature coverage"))
    axes[1].set_xlabel(_replace_underscores("Coverage"))
    axes[1].set_ylabel(_replace_underscores("Count"))
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


def _figure_dimension_filtering(figures_dir: str, sig, formats: List[str]) -> Tuple[str, str]:
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
    ax.set_title(_replace_underscores("Std distribution across dimensions"))
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


def _figure_dimension_filtering_bar(figures_dir: str, sig, formats: List[str]) -> Tuple[str, str]:
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
    ax.set_title(_replace_underscores("Dimension filtering (orange = selected)"))
    ax.set_xlabel(_replace_underscores("Dimension index"))
    ax.set_ylabel(_replace_underscores("Std"))

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(7.5)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "dimension_filtering_bar", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_top_matches(figures_dir: str, sig, term_signatures, formats: List[str]) -> Tuple[str, str]:
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
    fig, ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=True)

    # Wrap or truncate long labels to avoid overflow
    labels_raw = [r[0].term for r in top]
    labels = _wrap_or_truncate_labels(labels_raw, max_chars_per_line=26, max_lines=2)
    distances = [float(r[1]) for r in top]
    cosines = [float(r[2]) for r in top]
    y = np.arange(len(top))

    ax.barh(y, distances, color="#0072b2")
    ax.invert_yaxis()
    ax.set_xlabel(_replace_underscores("Variance weighted distance"))
    ax.set_ylabel("")
    ax.set_yticks(y)
    ax.set_yticklabels([_replace_underscores(l) for l in labels], fontsize=7.5)
    ax.set_title(_replace_underscores("Top matches for sample signature"))

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


def _figure_similarity_distributions(figures_dir: str, sig, term_signatures, formats: List[str]) -> Tuple[str, str]:
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
    axes[0].set_title(_replace_underscores("Distance distribution across all terms"))
    axes[0].set_xlabel(_replace_underscores("Variance weighted distance"))
    axes[0].set_ylabel(_replace_underscores("Density"))

    if sns is not None:
        sns.histplot(cosines, bins=50, ax=axes[1], color="#d55e00", kde=True, stat="density")
    else:
        axes[1].hist(cosines, bins=50, color="#d55e00", edgecolor="#333333", density=True)
    axes[1].set_title(_replace_underscores("Cosine similarity distribution across all terms"))
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


def _figure_distance_cosine_hexbin(figures_dir: str, sig, term_signatures, formats: List[str], top_k: int = 15) -> Tuple[str, str]:
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

    ax.set_title(_replace_underscores("Distance vs cosine hexbin"))
    ax.set_xlabel(_replace_underscores("Variance weighted distance"))
    ax.set_ylabel(_replace_underscores("Cosine similarity"))
    ax.set_ylim(-1.0, 1.0)
    ax.legend(loc="upper right", frameon=False)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "distance_cosine_hexbin", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_pca_term_space(figures_dir: str, sig, term_signatures, formats: List[str], top_k: int = 20) -> Tuple[str, str]:
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from src.pipeline import match_input_to_terms

    # Stack term means
    means = np.stack([np.asarray(t.mean, dtype=float) for t in term_signatures], axis=0)
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(means)

    # Save input data
    data_csv_path = os.path.join(figures_dir, "pca_term_space.data.csv")
    _ensure_dir(data_csv_path)
    with open(data_csv_path, "w", encoding="utf-8") as f:
        f.write("term_id,library,term,x,y,is_top\n")
        ranked = match_input_to_terms(sig, term_signatures, dims_filter="quantile", q_low=0.3)
        top_ids = set(id(r[0]) for r in ranked[:top_k])
        for (t, xy) in zip(term_signatures, coords):
            is_top = int(id(t) in top_ids)
            term_clean = t.term.replace("\"", "'")
            f.write(f"{t.term_id},{t.library},\"{term_clean}\",{xy[0]},{xy[1]},{is_top}\n")

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

    ax.set_title(_replace_underscores("PCA of term means (top matches highlighted)"))
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(loc="best", frameon=False)

    _disable_grid(fig)
    saved = _save_figure(fig, figures_dir, "pca_term_space", formats)
    plt.close(fig)
    primary = saved[formats[0]]
    return primary, data_csv_path


def _figure_signature_heatmap(figures_dir: str, sig, formats: List[str]) -> Tuple[str, str]:
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
    axes[0].set_title(_replace_underscores("Signature mean by dimension"))
    axes[0].set_yticks([])
    axes[0].set_xlabel(_replace_underscores("Dimension index"))
    divider0 = make_axes_locatable(axes[0])
    cax0 = divider0.append_axes("right", size="2.5%", pad=0.10)
    fig.colorbar(im0, cax=cax0)

    # Heatmap of std
    im1 = axes[1].imshow(std[np.newaxis, :], aspect="auto", cmap="magma", extent=[1, d+1, 0, 1])
    axes[1].set_title(_replace_underscores("Signature std by dimension"))
    axes[1].set_yticks([])
    axes[1].set_xlabel(_replace_underscores("Dimension index"))
    divider1 = make_axes_locatable(axes[1])
    cax1 = divider1.append_axes("right", size="2.5%", pad=0.10)
    fig.colorbar(im1, cax=cax1)

    # Panel labels
    axes[0].text(-0.08, 1.05, "A", transform=axes[0].transAxes, fontsize=9, fontweight="bold")
    axes[1].text(-0.08, 1.05, "B", transform=axes[1].transAxes, fontsize=9, fontweight="bold")

    # Adjust margins to prevent text cutoff on right
    plt.subplots_adjust(left=0.08, right=0.90, top=0.90, bottom=0.18)

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

    outputs: List[Tuple[str, str]] = []
    # Pipeline overview removed per user request
    outputs.append(_figure_term_signature_distributions(figures_dir, term_signatures, formats))
    outputs.append(_figure_dimension_filtering(figures_dir, sig, formats))
    outputs.append(_figure_dimension_filtering_bar(figures_dir, sig, formats))
    outputs.append(_figure_top_matches(figures_dir, sig, term_signatures, formats))
    # Advanced/new figures
    outputs.append(_figure_similarity_distributions(figures_dir, sig, term_signatures, formats))
    outputs.append(_figure_distance_cosine_hexbin(figures_dir, sig, term_signatures, formats))
    outputs.append(_figure_pca_term_space(figures_dir, sig, term_signatures, formats))
    outputs.append(_figure_signature_heatmap(figures_dir, sig, formats))

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


