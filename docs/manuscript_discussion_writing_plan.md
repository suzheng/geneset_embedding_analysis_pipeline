# Discussion: Step-by-Step Writing Plan (Clinician-Facing Interpretation, No New Results)

## Scope and guardrails for this section

This Discussion plan is designed to interpret the Results for clinical physicians while maintaining strict evidence boundaries. The Discussion may explain plausible clinical and biological context, but it must not introduce any new ranked findings, any new term names that were not already reported in Results, or any new quantitative claims. Every specific pathway name, cell type name, drug name, or disease comparison that appears in Discussion must have already appeared in Results and must be traceable to a Results figure or table.

External statements about mechanisms, clinical practice, or prior studies must be supported by citations. If a citation is not available at drafting time, include a placeholder citation marker rather than inventing claims.

## Step-by-step drafting workflow

### Step 1: Begin with a brief, clinician-oriented summary of what the Results reported, without repeating the entire Results

Write the opening Discussion paragraph as a high-level synthesis that mirrors the Results subsection headings, using non-quantitative language. The purpose is to remind the reader of the categories of findings, such as pathways, immune cell types, drug signatures, and disease signature comparisons, without re-listing all terms and without adding new terms.

### Step 2: Interpret the pathway-level findings in clinical and pathophysiologic context

Using only the pathway names already present in Results, write one paragraph explaining how these pathways relate to host immune response and sepsis pathophysiology in clinician-recognizable language. Avoid describing intracellular cascade details unless they matter clinically. If specific cytokine pathways are discussed, anchor them to established sepsis immunology literature and avoid overstating causality.

### Step 3: Interpret the immune cell type findings as hypotheses about cellular contributors

Using only the cell types already reported in Results, write one paragraph that explains what it could mean clinically when a sepsis-associated gene set is similar to specific immune cell marker signatures. Keep the statements framed as hypotheses or associations rather than definitive conclusions, and avoid implying that the analysis measures cell proportions unless the pipeline actually uses single-cell deconvolution, which it does not.

### Step 4: Interpret the drug signature findings as prioritization for further evaluation, not treatment recommendations

Using only the drug names already reported in Results, write one paragraph explaining how similarity to drug perturbation signatures can be used for hypothesis generation and prioritization. Make clear that signature similarity is not evidence of clinical efficacy or safety. If you group drugs by class for readability, ensure the grouping is accurate and supported by pharmacology references, and avoid suggesting off-label use.

### Step 5: Interpret the disease signature comparisons as shared-response hypotheses and differential diagnosis context

Using only the disease or condition signatures already reported in Results, write one paragraph explaining what shared signatures could imply about overlapping host-response programs. Avoid claims that the analysis diagnoses sepsis or distinguishes pathogens unless there is explicit evaluation data, which would need to be presented in Results.

### Step 6: Explain what the embedding-based approach adds relative to traditional enrichment, using clear clinician language

Write a short paragraph that summarizes, in non-technical terms, why embedding-based “soft matching” can complement traditional gene overlap methods. Use the framing from `sepsis/shared_doc/BACKGROUND.md` but remove deep learning implementation details. Emphasize interpretability and continuous similarity while acknowledging that embedding methods are hypothesis-generating and depend on the training data of the embedding model.

### Step 7: State limitations precisely and conservatively

Write a limitations paragraph that includes only limitations that are true for this analysis and can be supported by the pipeline design. Typical limitations include dependency on the input gene list definition, potential bias from embedding training data, the fact that similarity metrics do not establish causality, and the need for experimental or clinical validation. If any run-specific limitation is relevant, it must be supported by QC outputs in `*_qc.csv` or by run metadata.

### Step 8: Provide concise, actionable future directions appropriate for a clinician-facing manuscript

Write a short paragraph describing pragmatic next steps, such as validating top pathway and drug hypotheses in independent transcriptomic cohorts or experimental models, and prioritizing candidates for translational evaluation. Keep future directions general and feasible, and do not propose specific interventions unless the manuscript has sufficient evidence and citations.

### Step 9: Close with a restrained conclusion consistent with the data boundary

Draft a final paragraph that restates the contribution at the level of method and interpretability, without adding new findings. The conclusion should mirror the objective described in the Introduction and remain consistent with the Results.

## Final Discussion checklist before you move on

Confirm that every proper noun in Discussion that refers to a specific result also appears in Results, that no new numbers or rankings appear, and that any mechanistic or clinical claims are cited or marked with a placeholder citation. Confirm that the tone stays clinically useful and does not oversell the analysis.


