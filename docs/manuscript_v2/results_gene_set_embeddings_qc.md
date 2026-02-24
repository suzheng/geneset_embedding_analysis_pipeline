## Results

### Gene set embedding analysis and quality control

To investigate sepsis biology, we focused on genes with established genetic associations with the disease. Specifically, we used a curated set of 23 genes derived from a large systematic review and meta-analysis of host genetic variants linked to sepsis susceptibility [https://doi.org/10.1186/s13054-019-2313-0]. That study integrated 349 publications and 405 variants across 172 genes, ultimately identifying 29 variants within 23 genes significantly associated with sepsis risk. These genes represent key components of pathogen sensing, inflammatory signaling, and host defense.

Although these genes have known individual functions, genes often participate in multiple biological processes. To better understand their shared biology, we analyzed the 23 genes as a group rather than individually. Traditional gene set enrichment analysis tests whether disease genes overlap with known functional gene sets [https://doi.org/10.1073/pnas.0506580102]. However, this approach is less effective for small gene sets because biologically related genes may not share identical gene symbols, even when they participate in similar pathways.

To overcome this limitation, we used numerical gene representations in which genes with similar biological functions are placed close to each other in a mathematical space. These representations were generated using GeneRAIN, a deep learning model trained on 410,000 gene expression samples, providing robust and generalizable functional representations [https://doi.org/10.1186/s13059-025-03749-6]. Instead of requiring gene overlap, this approach evaluates whether two gene sets occupy similar regions in biological feature space, even if they share no common genes.

After embedding the sepsis-associated genes, we next identified which embedding dimensions best captured their shared biology. Each gene is represented across 200 dimensions, but not all dimensions contribute equally to defining sepsis biology. Some dimensions reflect shared functional signals, while others capture biological heterogeneity or noise.

To improve interpretability, we applied variance-based dimension filtering. Dimensions were ranked by variability across the sepsis genes; lower variability suggests shared biological behavior. We retained dimensions within the lowest 30% of variability, corresponding to 60 of 200 dimensions (Figure 1). The remaining 140 higher-variability dimensions were excluded to reduce noise and improve similarity detection [https://doi.org/10.1093/bib/bbab354].

![Dimension filtering visualization](../../sepsis/sepsis_Reactome_Pathways_2024/figures/combined_figure_qc.png)

*Figure 1. Dimension filtering for sepsis gene embedding analysis. (A) Distribution of standard deviation across embedding dimensions, sorted in ascending order. The dashed line indicates the 30th percentile threshold. (B) Selection status across dimensions. Retained dimensions (orange) show low variability and stronger functional coherence, while filtered dimensions (blue) show higher variability.*

Of the 23 sepsis-associated genes, 21 were successfully mapped to the embedding space, achieving 91.3% coverage. This high coverage suggests that most genetically supported sepsis susceptibility genes are well represented in the learned biological space, providing a reliable basis for further analysis.

The resulting sepsis embedding signature summarizes the collective functional profile of these genes across the 60 selected dimensions. This signature was then compared with curated libraries of pathways, immune cell types, diseases, and drug perturbation profiles to identify biological processes and therapeutic hypotheses with similar functional patterns [https://doi.org/10.1093/nar/gkad1025, https://doi.org/10.1093/nar/gkac947, https://doi.org/10.1093/bioinformatics/btv313, https://doi.org/10.1016/j.patter.2024.101072, https://doi.org/10.1038/ncomms12846].

By focusing on biologically consistent dimensions and maintaining high gene coverage, we improved signal clarity and ensured that downstream similarity results more likely reflect true biological mechanisms rather than statistical noise.