# geneset_embedding_analysis_pipeline

Embedding-based gene-set signature matching using [GeneRAIN](https://doi.org/10.1186/s13059-025-03749-6) representations. This repository contains the code, precomputed reference libraries, and committed outputs for the sepsis susceptibility gene manuscript (*Human Mutation* revision; release tag **`human-mutation-revision-1`**).

## Data

| Resource | Location |
|----------|----------|
| GeneRAIN 200d embeddings | [Zenodo 10408775](https://zenodo.org/records/10408775) — download and place as `data/GeneRAIN-vec.200d.txt` (or set `--embeddings-path`) |
| Reference term signatures | `precomputed/*.signatures.jsonl` (included in repo) |
| Manuscript analysis outputs | `sepsis/` (CSVs, `*_PARAMS.md`, `*_run_info.json`, per-library `figures/`) |

## Environment (UV)

Use the **minimal** dependency set for this paper. The root `requirements.txt` is a legacy full environment (metabolomics tools and private git dependencies) and is **not** required to reproduce the sepsis embedding analysis.

```bash
cd geneset_embedding_analysis_pipeline
uv venv --python 3.10
source .venv/bin/activate   # Windows: .venv\Scripts\activate
uv pip install -r requirements-pipeline-minimal.txt
# or: uv pip install -e .
```

## Reproduce rankings (CLI)

23-gene list (Lu et al., 2019 susceptibility meta-analysis), matching `scripts/run_real_data_example.py`:

```text
TLR1,TLR2,LBP,AGER,NOD2,MBL2,TNF,LTA,IL1B,IL6,CXCL8,IL10,IFNG,CXCL1,CXCL12,MIF,FCGR2A,SERPINE1,ACE,LCE4A,TAGAP,VDR,MIR608
```

Example — one library (Reactome), same defaults as the manuscript (`--random-state 17`, `--q-low 0.3`, `--dims-filter quantile`, `--top-k 20`, no clustering for small gene sets):

```bash
python -m src.main \
  --embeddings-path data/GeneRAIN-vec.200d.txt \
  --entity-id sepsis \
  --associated-item sepsis \
  --input-genes "TLR1,TLR2,LBP,AGER,NOD2,MBL2,TNF,LTA,IL1B,IL6,CXCL8,IL10,IFNG,CXCL1,CXCL12,MIF,FCGR2A,SERPINE1,ACE,LCE4A,TAGAP,VDR,MIR608" \
  --precomputed-terms-jsonl precomputed/Reactome_Pathways_2024.signatures.jsonl \
  --dims-filter quantile \
  --q-low 0.3 \
  --top-k 20 \
  --no-clustering \
  --random-state 17
```

Run the five manuscript libraries by pointing `--precomputed-terms-jsonl` to:

- `precomputed/Reactome_Pathways_2024.signatures.jsonl` → Figure 1E / QC panels
- `precomputed/CellMarker_2024.signatures.jsonl` → Figure 1C–D
- `precomputed/RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO.signatures.jsonl` → Figure 2
- `precomputed/DSigDB.signatures.jsonl` → Figure 3
- `precomputed/RummaGEO_DrugPerturbations_2025.signatures.jsonl` → Figure 4

Outputs are written under `sepsis/sepsis_<LibraryName>/` (or use `--output-folder` to override). **Committed results in `sepsis/` match the submitted manuscript**; re-running should be consistent when using the same embeddings file and random seed.

## Regenerate combined figures

After pipeline runs (or using existing `sepsis/*/figures/`):

| Manuscript figure | Script |
|-------------------|--------|
| Figure 1 (QC + cell + pathway panels) | `scripts/combine_gene_set_embeddings_qc_figures.py`, `scripts/combine_cellmarker_figures.py` (see script headers for paths) |
| Figure 2 (disease GEO) | `scripts/combine_rna_seq_disease_geo_figures.py` |
| Figure 3 (DSigDB) | `scripts/combine_dsigdb_figures.py` |
| Figure 4 (RummaGEO drugs) | `scripts/combine_rummageo_drugperturbations_2025_figures.py` |

```bash
python scripts/combine_gene_set_embeddings_qc_figures.py
python scripts/combine_cellmarker_figures.py
python scripts/combine_rna_seq_disease_geo_figures.py
python scripts/combine_dsigdb_figures.py
python scripts/combine_rummageo_drugperturbations_2025_figures.py
```

## Tests

```bash
uv pip install pytest
pytest tests/ -q
```

## Citation

If you use this pipeline, please cite the GeneRAIN paper and the manuscript under revision. Report issues at https://github.com/suzheng/geneset_embedding_analysis_pipeline/issues.
