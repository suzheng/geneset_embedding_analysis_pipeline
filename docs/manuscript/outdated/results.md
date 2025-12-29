# Results

## Analysis Completeness and Quality Control

The gene set embedding analysis was performed on a clinical query signature consisting of 23 input genes. Quality control metrics recorded in the pipeline output (`sepsis_Reactome_Pathways_2024_geneset_embed_qc.RESULTS.md`) indicated that 17 out of 23 genes (73.9%) were successfully mapped to the embedding space across all six analyzed libraries. The remaining 6 genes were not found in the precomputed embeddings and were excluded from the distance calculations. Clustering analysis was disabled for this run.

## Pathway and Network Matches

The query signature was compared against the Reactome Pathways 2024 library. Terms were ranked by variance-weighted distance (ascending). The highest-ranked pathway matches included "Chemokine Receptors Bind Chemokines" (distance: 3.09), "Signaling by Interleukins" (distance: 3.22), "TNFs Bind Their Physiological Receptors" (distance: 3.43), "Interleukin-20 Family Signaling" (distance: 3.51), and "Cytokine Signaling in Immune System" (distance: 3.73). Figure 1 shows the top 10 matches for this library.

## Immune Cell Type Associations

Comparison with the CellMarker 2024 library identified several cell-type specific signatures. The top-ranked terms by distance included "Cancer Stem Cell Undefined Human" (4.06), "Macrophage Lung Mouse" (4.10), "Basophil Undefined Human" (4.16), "Stage II Neutrophil Bone Marrow Mouse" (4.20), and "Myeloid Cell Bone Marrow Human" (4.25). Other significant associations included "Exhausted CD8+ T Cell Liver Human" (4.25) and "Neutrophil Blood Human" (4.30). These results are visualized in Figure 2.

## Drug Candidate Identification

Drug matches were extracted from three libraries: DSigDB, RummaGEO Drug Perturbations 2025, and LINCS L1000 Chemical Perturbation Consensus Sigs. 

In DSigDB, the top-ranked terms were "lanatoside C MCF7 DOWN" (distance: 10.73), "digoxin MCF7 DOWN" (11.20), "bisacodyl MCF7 DOWN" (11.78), "metoprolol HL60 DOWN" (12.43), and "clobetasol MCF7 DOWN" (16.75). 

RummaGEO Drug Perturbations identified "REMDESIVIR HUMAN DN" (distance: 11.67), "PANOBINOSTAT MOUSE DN" (14.79), "AZACITIDINE HUMAN DN" (18.13), "SULFORAPHANE MOUSE DN" (18.72), and "BLEOMYCIN MOUSE UP" (24.49) among the highest-ranked entries. 

The LINCS L1000 library identified "EX-527 Up" (16.05), "Wortmannin Down" (20.46), "Palbociclib Down" (26.27), and "GR-127935 Down" (13.20) as top matches. Combined drug results are summarized in Figure 3.

## Disease Signature Comparisons

The query signature was compared against disease signatures from GEO. The highest-ranked matches were "Rfx2 Testes Knockout GSE68283 up" (distance: 10.56), "STAU1 HEK293T knockdown GSE52447 down" (17.72), "ZNA MCF7 30uM-3hrs GSE59251 down" (16.72), "FBXL10 Female - Whole embryo Ablation GSE76552 down" (19.60), and "miR-92a-1 Whole embryo Knockout GSE63813 down" (15.17). Figure 4 illustrates these associations.

---

## Figures and Tables

**Figure 1. Top Pathway Matches (Reactome 2024).** Horizontal bar chart showing the 10 most similar terms from the Reactome Pathways 2024 library. Terms are ranked by variance-weighted distance (x-axis, ascending). Source: `sepsis/sepsis_Reactome_Pathways_2024/figures/top_matches_example.svg`.

**Figure 2. Top Cell Type Matches (CellMarker 2024).** Horizontal bar chart showing the 10 most similar terms from the CellMarker 2024 library. Terms are ranked by variance-weighted distance (x-axis, ascending). Source: `sepsis/sepsis_CellMarker_2024/figures/top_matches_example.svg`.

**Figure 3. Drug Signature Matches.** Summary visualization of top drug candidates from DSigDB, RummaGEO, and LINCS libraries. Terms are ranked by ascending distance. Source data: `sepsis/sepsis_DSigDB/sepsis_DSigDB_geneset_embed_sepsis_DSigDB_term_matches.csv`.

**Figure 4. Disease Signature Comparisons.** Top 10 matches from the RNA-Seq Disease Gene and Drug Signatures from GEO library. Source: `sepsis/sepsis_RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO/figures/top_matches_example.svg`.

**Table 1. Summary of Top Ranked Terms Across Selected Libraries.** This table lists the top 5 matches for each primary category. Full results are available in the repository as CSV files (e.g., `*_term_matches.csv`).

| Library | Top Ranked Term | Distance | Cosine | Coverage |
| :--- | :--- | :--- | :--- | :--- |
| Reactome | Chemokine Receptors Bind Chemokines | 3.09 | 0.811 | 1.000 |
| CellMarker | Cancer Stem Cell Undefined Human | 4.06 | 0.640 | 1.000 |
| DSigDB | bisacodyl MCF7 DOWN | 11.78 | -0.205 | 0.922 |
| RummaGEO Drug | BLEOMYCIN MOUSE UP | 24.49 | -0.372 | 1.000 |
| GEO Disease | Rfx2 Testes Knockout GSE68283 up | 10.56 | -0.282 | 0.791 |
