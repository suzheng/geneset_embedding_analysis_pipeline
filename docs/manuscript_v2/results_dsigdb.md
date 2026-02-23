## Drug Repurposing Analysis Identifies Immunomodulatory and Anti-inflammatory Therapeutic Candidates

To identify potential therapeutic opportunities, the sepsis-associated gene functional signature was compared against the DSigDB drug signature library [https://doi.org/10.1093/bioinformatics/btv313] containing 4,025 compound-related gene expression signatures. Drug matching was performed using embedding-based functional similarity metrics, where smaller distance values and higher cosine similarity indicate stronger biological similarity between drug-induced transcriptional patterns and the sepsis signature.

The top-ranked drug match was POVPC (distance: 1.761, cosine similarity: 0.912), an oxidized phosphatidylcholine lipid mediator with anti-inflammatory properties. This was followed by D-Monapterin (distance: 1.939, cosine similarity: 0.897), a compound involved in pteridine metabolism. Inhaled corticosteroids fluticasone (distance: 2.108, cosine similarity: 0.890) and beclomethasone (distance: 2.192, cosine similarity: 0.882) ranked third and fourth, respectively. Additional high-ranking matches included clodronate disodium (distance: 2.299, cosine similarity: 0.883), an osteoclast activity inhibitor with immunomodulatory effects, and sulfasalazine (distance: 2.380, cosine similarity: 0.875), a disease-modifying anti-rheumatic drug used in inflammatory diseases.

Additional candidate compounds included parthenolide (distance: 2.432, cosine similarity: 0.863), DL-mevalonic acid (distance: 2.495, cosine similarity: 0.871), ferrous sulfate (distance: 2.529, cosine similarity: 0.876), and dinoprostone (distance: 2.548, cosine similarity: 0.871). Overall, the drug matching pattern was dominated by anti-inflammatory, immunomodulatory, and metabolic pathway-related compounds, suggesting that multiple biological processes contribute to sepsis immune dysregulation.

---

### Figure 1. Drug repurposing analysis results from DSigDB library

(A) Drug candidates ranked by variance-weighted distance (horizontal bars; shorter bars indicate greater similarity).
(B) Two-dimensional density distribution of variance-weighted distance (x-axis) and cosine similarity (y-axis). Color intensity represents density across hexagonal bins, and top candidate drugs cluster in regions of low distance and high similarity.
(C) Principal component projection showing global distribution of drug signatures (light blue dots), top matches (green circles), and the sepsis signature (orange star).

---

## Functional Immune Signature Characteristics

Functional mapping of the sepsis signature demonstrated simultaneous representation of pro-inflammatory and immune exhaustion features. The signature showed strong similarity to T helper 17 (Th17)-related immune activation patterns, consistent with Th17-mediated neutrophil recruitment and cytokine amplification during severe infection.

At the same time, signatures associated with exhausted CD4+ and CD8+ T cells were also highly similar to the sepsis signature, indicating impaired adaptive immune function consistent with sepsis-associated immunoparalysis. Innate immune cell populations including monocytes, macrophages, and neutrophils were also strongly represented.

These results suggest that the sepsis signature captures the dual immune state of inflammatory activation and immune suppression characteristic of clinical sepsis.

---

## Molecular Similarity Between Sepsis and Pulmonary Fibrotic Disease

Disease signature mapping identified functional convergence between sepsis and Idiopathic Pulmonary Fibrosis, suggesting shared molecular mechanisms related to epithelial injury and abnormal tissue repair responses.

---

## JAK-STAT Pathway Targeting as a Candidate Therapeutic Strategy

Drug perturbation signature analysis identified ruxolitinib, a JAK1/JAK2 inhibitor, as a candidate therapeutic compound. The result supports potential therapeutic targeting of cytokine-mediated signaling pathways involved in sepsis immune dysregulation.
