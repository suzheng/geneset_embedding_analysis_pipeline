## Embedding Signature Pipeline — Design Plan

### 1) Goal and scope
- **Goal**: Given a user-provided gene list (e.g., 12 genes) presumed to share a phenotype, compute one or more embedding "signatures" that capture stable dimensions (low variability across the genes). Then, find the closest precomputed term signatures (from Enrichr-style libraries) to the input signature(s).
- **Scope**:
  - Read gene embeddings from a file like `data/GeneRAIN-vec.32d.txt` (and support other dimensions, e.g., 200D), where the first line declares counts and dimension.
  - Support Enrichr-style term files where column 1 is term name, columns 2..N are associated gene symbols: directory `
    /Users/suzheng/noBackupData/BioBot/backend/data/tools/tool_data/database/gene_enrichment_analysis/`.
  - Precompute per-term signature statistics (mean and std per dimension) to enable fast matching.
  - Provide similarity search from input signature(s) to term signatures.

### 2) Data sources and formats
- **Embeddings file** (example: `data/GeneRAIN-vec.32d.txt`):
  - First line: `<num_records> <num_dimensions>` (e.g., `31769 32`).
  - Subsequent lines: `<GENE_SYMBOL> <d1> <d2> ... <dD>` with whitespace-separated floats.
  - We will auto-detect dimension `D` from the first line to keep the pipeline dimension-agnostic.

- **Term libraries** (Enrichr-style TSVs):
  - Location: `
    /Users/suzheng/noBackupData/BioBot/backend/data/tools/tool_data/database/gene_enrichment_analysis/`.
  - Each file: first column is term (e.g., disease/pathway), columns 2..N are gene symbols associated with the term. Variable width rows are expected.
  - Examples: `KEGG_2021_Human.tsv`, `OMIM_Disease.tsv`, `GO_Biological_Process_2023.tsv`, etc. (actual filenames per provided directory listing).
  - We will treat empty trailing columns as absent genes and ignore them.

### 3) Definitions
- **Signature (of a gene set or cluster)**: Per-dimension statistics summarizing the set of gene embeddings:
  - `mean`: length-`D` mean vector over included genes.
  - `std`: length-`D` standard deviation vector over included genes.

- **Dimension selection (optional, std-based)**: Dimensions can be filtered dynamically at analysis time using rules on the `std` vector (e.g., lowest-variance quantile), rather than persisting any boolean mask.

- **Clusters**: If the input gene list is sufficiently large (exceeds a configurable threshold), genes may form multiple coherent subsets in embedding space. Each subset yields a separate signature.

### 4) Core computations
1) **Load embeddings index**
   - Read the first line to detect `(N, D)`.
   - Build a dictionary `symbol_to_vector: str -> np.ndarray[D]`.
   - Normalize gene symbols consistently (e.g., uppercase; strip whitespace). Optionally de-duplicate if the embedding file contains aliases.
   - Metrics: report `coverage = (#genes_found / #genes_requested)`.

2) **Collect input genes' vectors**
   - Input can be provided via a file (plain-text or CSV/TSV) or directly as a Python list to the API (no CLI dependency).
   - For each symbol, look up the embedding; drop symbols not present and record them as missing.
   - If fewer than `min_genes_for_signature` (default 5) are found, abort politely with a clear message.

3) **Optional clustering of input genes**
   - Trigger clustering if `len(genes_found) >= cluster_min_genes` (default 20).
   - Recommended default: robust clustering whose `k` does not need to be known, e.g., HDBSCAN (if available). Fallback: choose `k` via silhouette on KMeans for `k ∈ {2..5}` and accept if `silhouette >= min_silhouette`.
   - If clustering yields 1 cluster (or clustering disabled), compute a single signature. If multiple clusters, compute one signature per cluster.

4) **Compute signature statistics**
   - For each group (all genes or each cluster):
     - `mean`: arithmetic mean per dimension.
     - `std`: default to sample std per dimension with small `epsilon` floor for numerical stability.
   - Optionally support robust alternatives (configurable):
     - `robust_std = IQR / 1.349` or `median_abs_deviation * 1.4826`.

5) **Optional dimension filtering (std-based)**
   - Provide multiple interchangeable rules for selecting dimensions at analysis time (no persisted masks):
     - Quantile rule: choose dims where `std <= quantile(std, q_low)`.
     - CV rule: choose dims where `std / (|mean| + epsilon) <= cv_threshold`.
     - Z-score rule: choose dims whose std z-score ≤ `-z_threshold`.
     - Top-K rule: choose the `K` lowest-variance dims.
   - The chosen rule and parameters are applied dynamically when matching/signature analysis is performed.

### 5) Term signatures (precomputation)
- For each TSV file in the Enrichr directory:
  1. Parse rows: `term = row[0]`; `genes = row[1:]` (filter empties, normalize symbols).
  2. Lookup embedding vectors; require `min_genes_per_term` (default 5) after filtering to compute a signature.
  3. Compute per-dimension `mean` and `std` using the same statistics as for input sets (no persisted masks).
  4. Record additional metadata: `library_file`, `term`, `n_genes_total`, `n_genes_with_embeddings`, `coverage`.
  5. Store to a precomputed artifact for fast lookup:
     - Preferred: columnar `parquet` with columns `term_id` (hash of file+term), `library`, `term`, `D`, arrays for `mean`, `std`, `n_found`, `coverage`.
     - Secondary: `jsonl` for human inspectability.

### 6) Similarity and matching
- Compute similarity between an input signature `S = (μ_S, σ_S)` and a term signature `T = (μ_T, σ_T)` using a variance-aware distance. Recommended default:
  - **Variance-weighted squared distance**:
    - Optionally select dimensions using an std-based filter on `σ_S` (e.g., lowest-variance quantile). If no filter is used, all dims are considered.
    - For each dimension `i`, use variance weight `w_i = 1 / (σ_S,i^2 + σ_T,i^2 + ε)`.
    - Distance: `d(S,T) = Σ_i w_i (μ_S,i − μ_T,i)^2`.
  - Optionally compute a **cosine similarity** on `μ` restricted to the selected dims as a tie-breaker.
  - Return a ranked list of closest terms with distances/similarities and diagnostics (selected dims count, coverage, library source).

### 7) Outputs
- For input gene list:
  - One or more signatures with: `mean`, `std`, `n_genes`, `coverage`, clustering diagnostics.
  - Ranked term matches (top `K`, default 25) across selected libraries or all.
  - Optional CSV/JSON report and a compact human-readable Markdown summary.

- For precomputation:
  - Per-library parquet files under `precomputed/terms/` (e.g., `KEGG_2021_Human.signatures.parquet`).
  - An aggregated parquet index combining all libraries under `precomputed/index/signatures.parquet` for single-pass search.

### 8) Configuration (with reasonable defaults)
- Paths
  - `embeddings_path`: default `data/GeneRAIN-vec.32d.txt` (override-able).
  - `terms_dir`: default `
    /Users/suzheng/noBackupData/BioBot/backend/data/tools/tool_data/database/gene_enrichment_analysis/`.

- Thresholds & rules
  - `min_genes_for_signature = 5`
  - `cluster_min_genes = 20`
  - `dims_filter` (optional): e.g., `quantile` with `q_low = 0.3`, or `none`
  - `epsilon = 1e-8` for numerical stability
  - `min_genes_per_term = 5`
  - `top_k_matches = 25`

- Clustering (if enabled)
  - Method: `hdbscan` if available else `kmeans`
  - KMeans selection: `k ∈ {2..5}`, pick by silhouette ≥ 0.35; if none, use single cluster
  - Minimum cluster size: `max(5, 0.2 * len(genes))` for HDBSCAN (if used)

### 9) Edge cases and data hygiene
- Case sensitivity: normalize gene symbols to uppercase when matching.
- Missing genes: track and report; exclude from statistics.
- Duplicates: if duplicates occur, keep the first or average duplicates deterministically.
- Very small sets: below `min_genes_for_signature` → report insufficient data.
- Outliers: robust std options can help when one or two genes are extreme.
- Dimension-agnostic: always derive `D` from the embeddings header; validate each line has exactly `D` floats.

### 10) Performance considerations
- Build an in-memory dictionary for embeddings (fits for typical sizes like ~30K genes × 32/200 dims).
- Precomputation stored in parquet for vectorized I/O; lazy-load only needed columns.
- For matching across many terms, vectorize distance computation across all candidates using NumPy; if dims are filtered, restrict computations to those dims for speed.
- Optional: maintain a whitened mean vector per term (`μ_T / (σ_T + ε)`) to approximate the weighted metric; useful for ANN indexing in FAISS/Annoy if needed later.

### 11) Validation & QA
- Unit tests for: parsing embeddings header and rows, term TSV parsing, std-based dimension selection rules, distance computations.
- Sanity checks: the same gene list repeatedly yields the same signature; adding duplicated genes does not change results materially; removing a small number of genes does not radically change selected dims when filters are used.
- Spot-check: known pathways/diseases should rank high for curated gene lists.

### 12) Module design (no CLI dependency)
- `src/embeddings_io.py`: load embeddings; normalize symbols; expose `get_vectors(genes)`.
- `src/signature.py`: compute mean/std; provide std-based dimension selection utilities; robust std helpers (optional).
- `src/cluster.py`: clustering utilities (optional; graceful fallback if dependencies absent).
- `src/terms_io.py`: parse Enrichr TSVs; yield `(library, term, genes)`.
- `src/precompute_terms.py`: batch compute and persist term signatures (parquet/jsonl).
- `src/match.py`: variance-aware distance and ranking with optional std-based dimension filtering.
- `src/report.py`: render JSON/CSV/Markdown outputs.
- `src/pipeline.py`: high-level functions to compute input signatures, precompute term signatures, and perform matching.

### 13) Matching workflow (end-to-end)
1. Load embeddings (detect `D`).
2. Load input gene list, map to vectors, compute coverage.
3. If enough genes, cluster → produce 1..C signatures; else produce one signature.
4. For each signature, compute `mean` and `std` (no persisted masks).
5. Load precomputed term signatures (or compute on-the-fly if absent).
6. For each signature, compute stability-weighted distance to all terms; rank and report top `K`.
7. Output reports and diagnostics (missing genes, coverage, stable dims count, etc.).

### 14) Deliverables and artifacts
- `precomputed/terms/*.parquet` — per-library term signatures.
- `precomputed/index/signatures.parquet` — combined index (optional).
- `reports/<timestamp>/` — result bundles for input lists (JSON/CSV/Markdown).

### 15) Implementation notes
- Keep the code dimension-agnostic by deriving `D` from the embeddings header; never hardcode 32 or 200.
- Ensure numerical stability with `epsilon` wherever division by std or CV occurs.
- Persist the configuration used to construct each signature alongside outputs for reproducibility.

### 16) Near-term steps to implement
1. Implement embeddings loader and simple signature computation (mean/std + quantile-based stable dims).
2. Implement term TSV parser and per-term signature precomputation to parquet.
3. Implement stability-weighted distance and ranking across precomputed terms.
4. Add optional clustering; if absent, everything else should still work.
5. Wrap with CLI subcommands; produce JSON/CSV/Markdown outputs.


