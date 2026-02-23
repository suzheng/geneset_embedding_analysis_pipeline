# Results

## Gene Set, Embeddings and Quality Control

To interpret our sepsis-associated gene set in terms of known biology—such as pathways, cell types, or drug effects—we need a way to compare it with reference gene sets even when the gene lists do not overlap. Traditional overlap-based methods require shared genes; here we use gene embeddings instead. Embeddings are numerical representations of genes learned from large-scale expression data; genes with similar biological roles tend to have similar embedding values. By summarizing our gene set as an embedding signature (a centroid plus a measure of how coherent the genes are in each dimension), we can find reference terms with similar functional profiles. We further focus on dimensions where our genes agree—where they show consistent patterns—and down-weight dimensions where they disagree, which would otherwise add noise to the comparison. This section reports the quality control results for this embedding-based approach.

The analysis was performed on a curated set of 23 sepsis-associated genes representing key mediators of the immune response and inflammatory pathways implicated in sepsis pathophysiology. Gene symbols were normalized through standardized preprocessing procedures, including whitespace trimming and conversion to uppercase, to ensure consistent mapping to the reference embedding space. Following preprocessing, all 23 input genes were processed for embedding retrieval.

Quality control assessment revealed that 21 of the 23 input genes (91.3%) were successfully mapped to the GeneRAIN embedding dictionary, which contains 31,769 genes with precomputed 200-dimensional vector representations. The two genes that could not be mapped were excluded from downstream signature computation. This coverage rate of 91.3% indicates that the vast majority of sepsis-associated genes in the input set were represented in the embedding space, enabling highly reliable signature computation. The successfully mapped genes encompassed key immune mediators including pattern recognition receptors, cytokines, chemokines, and other inflammatory signaling molecules central to sepsis pathogenesis.

The 21 successfully mapped genes were converted into their corresponding 200-dimensional embedding vectors, where each dimension represents a learned biological feature derived from analysis of 410,000 human RNA-seq samples. To summarize the collective functional profile of the sepsis gene set, statistical signatures were computed across all embedding dimensions. For each of the 200 dimensions, two summary statistics were calculated: a mean vector representing the centroid position of the gene set in that dimension, and a standard deviation vector capturing the variability of genes within the set. The mean vector provides a measure of the central tendency of the sepsis genes' functional characteristics, while the standard deviation quantifies the coherence or diversity of the gene set across each semantic dimension.

To focus similarity calculations on the most stable and informative embedding dimensions, a variance-based dimension selection approach was applied. Dimensions were ranked by their standard deviation values, and those with variability at or below the 30th percentile threshold were retained for similarity matching. This quantile-based filtering identified 60 dimensions (30% of 200 total dimensions) that exhibited coherent behavior across the sepsis gene set, where genes showed convergent functional characteristics. The remaining 140 dimensions, which displayed higher variability and potentially introduced noise, were excluded from downstream similarity calculations. This dimension selection strategy enhances the signal-to-noise ratio by prioritizing dimensions where the sepsis genes exhibit consistent functional patterns.

![Dimension filtering visualization](../../sepsis/sepsis_Reactome_Pathways_2024/figures/combined_figure_qc.png)

*Figure 1. Dimension filtering for the sepsis gene set embedding analysis. (A) Distribution of standard deviation values across 200 embedding dimensions, sorted in ascending order. The horizontal dashed line indicates the 30th percentile threshold (approximately 0.071) used for dimension selection. Dimensions with standard deviation values at or below this threshold (60 dimensions) were retained for similarity calculations, while dimensions above the threshold (140 dimensions) were filtered out to reduce noise. (B) Bar chart showing the selection status of embedding dimensions across all 200 dimensions. Selected dimensions (60 dimensions, shown in orange) represent dimensions with low variability where the sepsis gene set exhibits coherent functional characteristics. Filtered dimensions (140 dimensions, shown in blue) represent dimensions with higher variability that were excluded from similarity calculations to enhance signal-to-noise ratio.*

**Table 1. Quality control metrics for the sepsis gene set embedding analysis**

| Metric | Value |
|--------|-------|
| Total input genes | 23 |
| Genes successfully mapped to embeddings | 21 |
| Genes not found in embedding dictionary | 2 |
| Gene coverage rate | 91.3% |
| Embedding dimensions per gene | 200 |
| Dimensions selected for similarity matching | 60 (30%) |
| Dimensions filtered out | 140 (70%) |

The computed embedding signature, consisting of mean and standard deviation vectors across the selected 60 dimensions, provides a quantitative representation of the sepsis gene set's functional profile in the embedding space. This signature serves as the basis for similarity matching against reference gene set libraries, enabling identification of biological pathways, immune cell types, drug signatures, and disease associations that share similar functional characteristics with the sepsis-associated genes. The high gene coverage rate (91.3%) ensures robust and reliable similarity calculations across downstream analyses.