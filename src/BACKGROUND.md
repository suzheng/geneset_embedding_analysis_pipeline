# Background Summary

### Gene set analysis (GSA)
Interprets omics data by mapping gene lists to pathways, GO terms, diseases, cell markers, drug targets etc.; widely used in biology and medicine.

---

# Limitations of Traditional (Hard-Matching) GSA

### Hard matching issue
Only counts exact overlaps; ignores gene functional similarity.

### Missed biology
- Paralogous genes (e.g., kinases, TF families) treated as unrelated.  
- Functionally redundant genes in parallel pathways ignored.  
- Proteins with similar domains/activities not linked.  

**Impact:** Misses meaningful associations, esp. novel/poorly annotated genes.

### Statistical biases
- **Small sets** → low power, underrepresented.  
- **Large sets** → inflated significance, annotation bias toward well-studied areas.  

**Result:** Size-dependent bias distorts enrichment interpretation.  

**The most important problem:** Possible incomplete interpretation of the input gene set, like a list of genes found associated with a specific phenotype or a disease.

---

# Deep Learning & Transformer Revolution

### Transformers (Vaswani et al., 2017)
Use self-attention to model long-range dependencies; core of BERT, GPT, AlphaFold2.

### Self-supervised learning
Learns from unlabeled data via masked/next prediction.

### Biology applications
- Gene expression modeling (cell identity, regulation, perturbation prediction).  
- Protein/DNA sequence modeling (ESM, DNABERT).  

**Benefit:** Capture complex, contextual biological relationships.

---

# GeneRAIN: Transformer-Based Gene Representation Learning

- Trained on **410k human bulk RNA-seq samples** via self-supervised learning.  
- Learned embeddings encode: protein domains, GO terms, disease links, pathways, TF targets, epigenetics, cell markers, PPIs.  
- Introduced **Binning-By-Gene normalization**; compared BERT vs GPT (**GPT better**).  
- **Evaluation:** 5,000+ classifiers tested; 72.65% (GPT) vs 20.36% (Geneformer) with AUC > 0.6.  
- **lncRNA transfer:** Predicted 62.5M attributes for 13,030 lncRNAs.  
- **Finding:** Distinct embedding dimensions encode different biological attributes; captures regulation & perturbation effects (validated by CRISPRi).

---

# Opportunity: Soft (Embedding-Based) Matching

Embedding-based analysis overcomes hard-matching limits:

- **Gene role similarity** captured by spatial proximity in embedding space.  
- **Annotation gap tolerance** — novel genes inferred via neighbors.  
- **Size-invariant metrics** — use mean/variance signatures instead of counts.  
- **Continuous similarity scores** — degree of relatedness, not binary.  
- **Multifaceted interpretation** — unified analysis across pathways, diseases, cell types, regulation.

---

# Toward Comprehensive Functional Mapping

- Current gene set libraries incomplete; many associations missing/undetected.  
- Hard matching overlooks subtle relationships.  

Embedding-based methods enable:
- Discovery of new gene set–phenotype links.  
- Detection of subtle/complex functional relations.  
- Richer annotation via multifaceted embeddings.  
- Cross-domain prediction (e.g., lncRNAs).  
- Data-driven hypothesis generation.

---

# This Study: Embedding-Based Gene Set Functional Analysis

### Pipeline overview
1. **Gene set embedding:** Retrieve GeneRAIN vectors → compute mean (centroid) + SD (variability).  
2. **Dimension selection:** Keep low-variance, high-signal dimensions.  
3. **Reference library:** Precompute signatures from KEGG, Reactome, GO, OMIM, ClinVar, DisGeNET, LINCS, DSigDB, CRISPR, cell markers, TF targets, protein complexes.  
4. **Similarity matching:** Rank by variance-weighted distance + cosine similarity.  
5. **Interpretation/visualization:** PCA, similarity plots, heatmaps.  
6. **Opposite-term detection:** Identify antagonistic/contrasting perturbations.

### Goal
- Reveal new functional associations via soft, continuous similarity metrics.  
- Integrate AI-learned embeddings into standard GSA workflows.  
- Enable more sensitive, nuanced, comprehensive functional interpretation, advancing gene function, disease mechanism, and therapeutic discovery.
