
## Immune Cell Types Associated With the Sepsis Signature

To characterize the immune landscape associated with sepsis, we compared the sepsis genetic signature against the CellMarker reference library [https://doi.org/10.1093/nar/gkac947], which contains curated marker profiles for 1,692 immune and non-immune cell populations. Cell type similarity was assessed using a variance-weighted distance metric in embedding space, where smaller distances indicate stronger functional similarity between the sepsis signature and reference cell type signatures.

The analysis revealed strong enrichment of both adaptive and innate immune cell populations, consistent with the complex immunological dysregulation observed in sepsis. Among adaptive immune cells, B cell signatures showed the strongest similarity to the sepsis signature (distance: 3.250; cosine similarity: 0.809), suggesting activation of humoral immune responses during sepsis. B cells contribute to host defense through antibody production and antigen presentation, and their activation may reflect ongoing pathogen recognition and immune adaptation.

T cell populations were also prominently represented. The sepsis signature showed strong similarity to T helper 17 (Th17) cells (distance: 3.278; cosine similarity: 0.865), a pro-inflammatory T cell subset that promotes neutrophil recruitment and amplifies cytokine signaling. 

Innate immune populations were also strongly represented. Monocyte signatures from bone marrow and peripheral blood compartments showed high similarity to the sepsis signature (Monocyte Bone Marrow Mouse: distance 3.491; Monocyte Peripheral Blood Human: distance 3.744). Macrophage signatures, particularly classically activated M1 macrophages from lung tissue, were also highly ranked (distance: 3.515; cosine similarity: 0.858), supporting the central role of pro-inflammatory myeloid activation in sepsis.

Additional innate immune components included mast cells and basophils, which are known to release histamine and inflammatory mediators during acute immune activation (Mast Cell Bone Marrow Mouse: distance 3.463; cosine similarity 0.811; Basophil Undefined Human: distance 3.749; cosine similarity 0.804). Tissue-level immune signatures were also observed, including lung immune cell signatures, reflecting organ-specific immune responses during systemic infection.

Overall, these findings demonstrate coordinated activation of both innate and adaptive immune compartments during sepsis, consistent with the dual inflammatory and immune-suppressive states observed clinically.

![Cell type analysis results](../../sepsis/sepsis_CellMarker_2024/figures/combined_figure.png)

*Figure 1. Immune cell type association analysis. (A) Cell type similarity ranking based on variance-weighted distance (shorter bars indicate greater similarity). (B) Principal component projection of reference signatures (light blue dots), top matching cell types (green circles), and the sepsis signature (orange star). Clustering of top matches near the sepsis signature indicates functional similarity in embedding space.*

**Table 1. Top matching cell types from CellMarker 2024 library**

| Rank | Cell Type | Distance | Cosine Similarity | Genes Found | Coverage |
|------|-----------|----------|-------------------|-------------|----------|
| 1 | B Cell Undefined Human | 3.250 | 0.809 | 96 | 1.000 |
| 2 | T Helper17 (Th17) Cell Undefined Human | 3.278 | 0.865 | 17 | 0.850 |
| 3 | Mast Cell Bone Marrow Mouse | 3.463 | 0.811 | 30 | 1.000 |
| 4 | Monocyte Bone Marrow Mouse | 3.491 | 0.792 | 68 | 0.986 |
| 5 | M1 Macrophage Lung Human | 3.515 | 0.858 | 23 | 1.000 |
| 6 | Monocyte Peripheral Blood Human | 3.744 | 0.794 | 65 | 0.956 |
| 7 | Basophil Undefined Human | 3.749 | 0.804 | 29 | 1.000 |
| 8 | T Cell Undefined Mouse | 3.775 | 0.795 | 35 | 1.000 |
| 9 | Myeloid Cell Bone Marrow Human | 3.906 | 0.781 | 120 | 1.000 |
| 10 | Immune Cell Lung Human | 3.923 | 0.785 | 23 | 1.000 |

The full ranked list of cell type similarities is provided in the supplementary data file:
[sepsis_CellMarker_2024_term_matches.csv](../../sepsis/sepsis_CellMarker_2024/sepsis_CellMarker_2024_geneset_embed_sepsis_CellMarker_2024_term_matches.csv). The dataset includes cell type name, tissue origin, species, distance metric, cosine similarity, embedding coverage, and marker gene mapping statistics.


