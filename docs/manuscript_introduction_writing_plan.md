# Introduction: Step-by-Step Writing Plan (Clinician-Facing, No Results)

## Scope and guardrails for this section

This Introduction plan is designed to produce an Introduction that motivates the study for clinical physicians without including any results, rankings, or outcome statements. The Introduction should focus on the clinical problem of sepsis, the gap in current molecular interpretation approaches, and the rationale for an embedding-based, soft-matching gene set analysis. Any factual claims about sepsis epidemiology or prior evidence must be referenced and must not be invented; if a citation is not available, leave a placeholder citation tag to be filled after literature verification.

The primary repo source for the conceptual background is `sepsis/shared_doc/BACKGROUND.md`, but its language must be rewritten for clinicians and for journal style, and it must be pruned to avoid overly technical transformer details that do not materially help the clinical narrative.

## Step-by-step drafting workflow

### Step 1: Open with the clinical motivation using standard journal framing

Draft the first paragraph as a clinical problem statement describing sepsis as a life-threatening syndrome with high morbidity and mortality and heterogeneous host responses. Keep this paragraph general and avoid any claims that would require a number unless you have a verified citation ready. End the paragraph by stating that improved mechanistic interpretation of sepsis-associated molecular signals could support hypothesis generation for pathophysiology and potential therapeutic strategies.

### Step 2: Introduce the role of gene signatures in sepsis without mentioning this study’s outputs

Write a paragraph explaining that gene expression and curated gene lists are commonly used to summarize host response, compare cohorts, and map signals to pathways, immune cell processes, or drug-related gene sets. Keep the language focused on clinician-understandable concepts such as pathways, immune programs, and drug signatures. Do not reference any specific pathway or drug names, because those would be results.

### Step 3: Explain the limitation of traditional enrichment in clinician-friendly terms

Using `sepsis/shared_doc/BACKGROUND.md` only as a conceptual reference, write a paragraph describing that many standard gene set analyses depend on exact gene overlap, which can miss related biology when different genes reflect similar functions or when annotations are incomplete. Keep the explanation concrete and short. Avoid highly technical examples such as paralog families unless you can express them in a single sentence without jargon.

### Step 4: Present the rationale for embedding-based “soft matching” as an enabling technology

Draft a paragraph that introduces the core idea of embedding-based gene representations as a way to summarize functional similarity beyond exact overlap. Use one short, plain-language sentence defining “embedding” as an AI-learned representation that places functionally related genes closer together in a learned space. Avoid describing transformer architectures, training objectives, or deep learning benchmarks unless they are directly necessary for comprehension; those details belong in Methods only if needed for reproducibility.

### Step 5: Define the study objective and what will be reported, in non-result terms

Write a final Introduction paragraph that states the objective of the manuscript. The objective should align with `docs/sepsis_manuscript_design.md` and should describe the categories of outputs that will be reported in Results without naming specific findings. A suitable objective statement is that the analysis maps a sepsis-associated gene set to clinically interpretable reference libraries, including pathways, immune cell marker signatures, drug-related signatures, and disease expression signatures, using an embedding-based similarity approach.

### Step 6: Place a brief roadmap sentence that matches the eventual section order

Close the Introduction with a single sentence that previews the manuscript structure, for example indicating that Methods describe the embedding signature approach and reference libraries, and Results report ranked matches across clinically relevant libraries. Do not allude to what was found.

## Final Introduction checklist before you move on

Re-read the Introduction and confirm there are no term names, drug names, pathway names, ranks, figures, or numbers that correspond to this study’s outputs. Confirm that the narrative is clinician-facing and that any factual claims that would require external validation are either cited or clearly marked with a placeholder citation to be filled later.


