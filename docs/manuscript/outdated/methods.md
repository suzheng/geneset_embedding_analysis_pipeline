# Methods

## Input Gene Set and Data Normalization

The primary analysis unit for this study is a curated sepsis-associated gene set, representing a collection of gene symbols associated with the pathophysiology of sepsis. Prior to analysis, all gene symbols were normalized through a standardized procedure involving whitespace trimming, conversion to uppercase, and removal of empty entries. Genes that could not be mapped to the reference embedding space were excluded from the analysis.

## Gene Embedding Model

Gene signatures were constructed using GeneRAIN embeddings, a high-dimensional semantic representation where each gene symbol is mapped to a fixed-length vector in a 200-dimensional space ($D=200$). This model allows for the representation of complex biological relationships as geometric distances and orientations within the embedding space.

## Embedding Signature Construction

For each gene set, an embedding "signature" was computed as a statistical summary of its constituent gene vectors. This signature consists of a centroid (mean vector, $\mu$) and a variability vector (standard deviation, $\sigma$) across all 200 dimensions. The centroid represents the collective semantic position of the gene set, while the standard deviation captures the semantic coherence or variability of the genes in each dimension. To ensure numerical stability, standard deviation values below $1 \times 10^{-8}$ were clamped to this minimum threshold.

## Stability-Based Dimension Selection

To focus the analysis on the most stable and biologically informative features, a variance-based dimension filter was applied. Embedding dimensions were ranked by their variability ($\sigma$), and only those dimensions with variability at or below the 30th percentile ($q=0.3$) were retained for similarity calculations. This quantile-based selection prioritizes dimensions where the input genes exhibit high semantic convergence, effectively down-weighting noisy or biologically irrelevant dimensions.

## Reference Gene Set Libraries and Precomputation

The query gene signature was compared against precomputed term signatures derived from curated gene set libraries. These libraries include clinical and biological categories such as Reactome Pathways, GO Biological Processes, CellMarker cell-type signatures, and disease-specific expression signatures. For each reference term, the associated genes were mapped into the same 200-dimensional GeneRAIN embedding space and summarized as term-level mean and variability vectors, enabling a direct geometric comparison with the query signature.

## Similarity Metrics and Term Ranking

Similarity between the query signature and reference terms was quantified using a variance-weighted distance that prioritizes dimensions identified as stable in the query set. Cosine similarity was also computed as a complementary measure of directional alignment in the embedding space. Reference terms were ranked primarily by their variance-weighted distance (ascending), with cosine similarity serving as a secondary criterion for ranking.

## Quality Control and Reproducibility

Quality control was performed by assessing the fraction of input genes successfully mapped to the embedding model (gene coverage). All analysis parameters, including the random seed (seed = 17) and software versions (Python 3.10.17, NumPy 2.2.6, scikit-learn 1.7.2), were recorded in run metadata to ensure technical reproducibility. The analysis pipeline produces systematic tabular outputs and data-linked visualizations for all results.
