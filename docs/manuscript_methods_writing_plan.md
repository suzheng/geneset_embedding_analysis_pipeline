# Methods: Step-by-Step Writing Plan (Clinician-Facing, Reproducible, No Results)

## Scope and guardrails for this section

This Methods plan is designed to produce a standard journal Methods section that is reproducible but not overly technical for clinical readers. The Methods text must not include any ranked terms, any “top matches,” any distances, any cosine similarities, or any outcome statements. Any numeric values or parameter settings that are reported must be copied from the run artifacts, not inferred.

The primary sources for Methods content in this repo are the shared method descriptions in `sepsis/shared_doc/SHARED_*.METHOD.md`, the per-library parameter files `sepsis/sepsis_<LIBRARY>/<prefix>_*.PARAMS.md`, and the run metadata `sepsis/sepsis_<LIBRARY>/<prefix>_run_info.json`. If the manuscript is limited to a subset of clinically focused libraries, the Methods must still describe the same pipeline while clearly stating which reference libraries were used for the primary analyses.

## Step-by-step drafting workflow

### Step 1: Define the analysis unit and the input gene set without implying any outcome

Start the Methods with a short paragraph defining what the input represents, using neutral language such as “a curated sepsis-associated gene set” or “a sepsis-related gene list,” and then immediately specify how it was obtained. Do not guess the source. If the provenance of the sepsis gene set is not explicit in the repo, stop and add a placeholder sentence that the gene list source will be inserted once verified, and then point to the run record that contains the actual list of genes for this run in `*_run_info.json`.

### Step 2: Specify the embedding model and embedding dimensionality from run artifacts

Describe the gene embedding model as GeneRAIN embeddings, and state that each gene symbol is mapped to a fixed-length vector. Confirm the embedding dimensionality and the embedding file used by reading `*_run_info.json` and the analysis parameters in `*_preprocess.PARAMS.md`. Report the dimensionality and any gene symbol normalization rules only if they are explicitly stated in those artifacts or in the shared docs.

### Step 3: Describe signature construction as mean and variability in embedding space

Write a concise paragraph explaining that the gene set is summarized as an embedding “signature” consisting of a mean vector and a variability vector across dimensions. Use wording consistent with the pipeline: centroid (mean) and within-set variability (standard deviation). Confirm the exact statistic definitions and any numerical stability constants by reading `sepsis/shared_doc/SHARED_signature.METHOD.md` and `*_signature.PARAMS.md`. Do not include any values computed from the run in this section.

### Step 4: Describe dimension selection as a stability filter, including the rule and threshold

Explain that the analysis restricts similarity calculations to relatively stable embedding dimensions, using a variance-based filter to down-weight noisy dimensions. Identify the dimension filter type used in this run (for example, quantile-based selection) and its threshold by reading `*_signature.PARAMS.md` and `*_run_info.json`. Report the rule and its parameter values exactly as recorded, but do not report the number of selected dimensions here unless it is part of the method definition rather than a run outcome; if the number selected is a computed outcome for this specific gene set, it belongs in Results, not Methods.

### Step 5: Define the reference gene set libraries used for matching, in clinician-friendly terms

State that the query signature was compared against precomputed term signatures derived from curated gene set libraries. In the main manuscript, restrict the primary analyses to clinician-interpretable library types described in `docs/sepsis_manuscript_design.md`, such as pathways, immune cell markers, drug signatures, and disease expression signatures. Confirm the exact library names used in the final manuscript by reading the output folder names under `sepsis/` and by verifying that each included library has a corresponding `sepsis_<LIBRARY>` directory with `*_term_matches.csv`. If the Methods mentions an optional library (for example, LINCS chemical perturbations), include it only if it will also appear in Results.

### Step 6: Describe term signature precomputation without over-technical implementation details

Explain that for each reference term, the associated genes were mapped into the same embedding space and summarized as term-level mean and variability vectors, enabling direct comparison to the query signature. Pull the precise definition from `sepsis/shared_doc/SHARED_term_matching.METHOD.md` and the pipeline description in `docs/embedding_signature_pipeline_plan.md`. Avoid describing internal file naming conventions or script names.

### Step 7: Define similarity metrics and ranking criteria exactly as implemented

Describe that similarity is quantified using a variance-aware distance between mean vectors that prioritizes stable dimensions, and that cosine similarity is also computed as a complementary directional measure. Confirm the formal descriptions in `sepsis/shared_doc/SHARED_term_matching.METHOD.md` and the implementation notes in `docs/embedding_signature_pipeline_plan.md`. If the pipeline computes both metrics but uses one for ranking, state that ranking was based on the recorded primary metric, and specify the metric name as it appears in the `*_term_matches.csv` header.

### Step 8: Describe quality control and reporting outputs without presenting outcomes

Write a paragraph describing that quality control includes gene coverage (fraction of input genes mapped to embeddings) and any checks or diagnostics captured in `*_qc.csv` and `*_qc.RESULTS.md`. In Methods, describe what is assessed and where it is recorded, but do not state the actual coverage values. Briefly state that the pipeline produces tabular outputs and publication-quality figures for transparency and reproducibility, and that each figure has an accompanying data file, as indicated by the presence of `figures/*.data.csv` and `figures/*.data.json` in each library output directory.

### Step 9: Software and reproducibility statement based on recorded metadata

Close Methods with a short reproducibility paragraph that references the existence of run metadata and parameter files. Use `*_run_info.json` and `*_*.PARAMS.md` to populate the analysis date, random seed, and key configuration values if they are explicitly recorded. Do not invent version numbers; include package versions only if they are captured in an artifact or documented in `requirements.txt` and can be cited accurately without inference.

## Final Methods checklist before you move on

Read the completed Methods once, line by line, and confirm that every specific parameter value is traceable to `*_run_info.json` or a `*.PARAMS.md` file, that no Results content has leaked in, and that the language stays accessible to clinicians by avoiding overly granular bioinformatics implementation detail.


