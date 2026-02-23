

## Functional Mapping of Sepsis Gene Signature to Drug Perturbation Signatures

To further identify potential therapeutic compounds associated with the sepsis molecular signature, we also compared the sepsis-associated gene signature against 5,344 drug perturbation signatures from the RummaGEO Drug Perturbations library, which contains experimentally derived transcriptional responses to drug exposure. The analysis identified several clinically relevant compounds with strong transcriptional similarity to the sepsis signature, suggesting potential drug repurposing opportunities for immune modulation in sepsis.

The highest-ranked match was glutamine (distance: 3.002, cosine similarity: 0.769), an amino acid that plays a central role in immune cell metabolism and cytokine regulation. This finding is consistent with prior studies demonstrating altered glutamine metabolism during critical illness and sepsis. Clinical trials of glutamine supplementation in critically ill patients have reported heterogeneous outcomes, reflecting the context-dependent nature of immune metabolic support in sepsis.

The second-ranked match was temozolomide (distance: 3.391, cosine similarity: 0.742), an alkylating chemotherapeutic agent primarily used in brain tumor treatment. Although not directly related to sepsis pathogenesis, temozolomide is known to induce immunosuppression through cytokine and lymphocyte regulatory effects, which may explain partial transcriptional similarity with the sepsis signature.

The third-ranked compound was bortezomib (distance: 3.448, cosine similarity: 0.739), a proteasome inhibitor that modulates nuclear factor-kappa B signaling and inflammatory cytokine production. Previous experimental sepsis models have shown context-dependent effects of proteasome inhibition on inflammatory injury, likely reflecting the dual inflammatory and immunosuppressive nature of sepsis.

Calcitriol, the active form of vitamin D, appeared twice among the top-ranked signatures, with similarity scores of distance 3.465 (cosine similarity 0.737) and 3.602 (cosine similarity 0.720) under different experimental conditions. Vitamin D signaling is known to influence innate and adaptive immune responses, and vitamin D deficiency has been associated with increased susceptibility to infection and poorer outcomes in critically ill patients.

Ruxolitinib, a Janus kinase (JAK) 1/2 inhibitor used in myeloproliferative disorders, ranked fifth (distance: 3.496, cosine similarity: 0.720). Ruxolitinib suppresses JAK/STAT-mediated cytokine signaling, including interleukin-6 and interferon-mediated pathways. Preclinical models of hyperinflammatory syndromes have demonstrated that JAK inhibition can reduce cytokine storm responses and attenuate organ injury when administered in appropriate therapeutic windows.

Additional high-ranking matches included curcumin (distance: 3.594, cosine similarity: 0.755), azacitidine (distance: 3.633, cosine similarity: 0.736), decitabine (distance: 3.658, cosine similarity: 0.732), and berberine (distance: 3.712, cosine similarity: 0.722). These compounds represent diverse pharmacological classes but share common biological themes of immune regulation, inflammatory suppression, or antimicrobial activity.

Overall, drug perturbation mapping revealed enrichment of compounds targeting metabolic regulation, epigenetic modification, and cytokine signaling pathways, supporting the hypothesis that sepsis pathophysiology reflects integrated dysregulation across multiple immune regulatory axes.

## Figure 1. Drug Perturbation Signature Landscape

![Combined analysis of drug perturbation signatures](../../sepsis/sepsis_RummaGEO_DrugPerturbations_2025/figures/combined_figure.png)

*Figure 1. Analysis of drug perturbation signatures from the RummaGEO Drug Perturbations 2025 library. (A) Drug perturbation signatures ranked by variance-weighted distance, where shorter bars indicate greater similarity to the sepsis signature. Drug names and experimental conditions are shown on the y-axis. (B) Joint density distribution of variance-weighted distance (x-axis) and cosine similarity (y-axis). Color intensity reflects density in hexagonal bins. Top-ranking perturbations (magenta points) cluster in the low-distance, high-similarity region. (C) Principal component analysis projection of drug perturbation signatures. Light blue dots represent all drug signatures, green circles represent top-ranked matches, and the orange star represents the sepsis signature.*

## Table 1. Top Drug Perturbation Matches

The complete ranked list of drug perturbation matches is provided in the supplementary data file. The table includes drug signature identifiers, experimental descriptions, variance-weighted distance, cosine similarity, embedding coverage, and gene mapping statistics.

The full dataset is available at:
[sepsis_RummaGEO_DrugPerturbations_2025_geneset_embed_sepsis_RummaGEO_DrugPerturbations_2025_term_matches.csv](../../sepsis/sepsis_RummaGEO_DrugPerturbations_2025/sepsis_RummaGEO_DrugPerturbations_2025_geneset_embed_sepsis_RummaGEO_DrugPerturbations_2025_term_matches.csv)
