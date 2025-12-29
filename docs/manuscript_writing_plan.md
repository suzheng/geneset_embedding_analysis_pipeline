# Sepsis Gene Set Embedding Manuscript: Step-by-Step Writing Plan (Index)

This folder contains a section-by-section writing plan for a clinician-facing journal-style manuscript based on the analysis design in `docs/sepsis_manuscript_design.md` and the existing analysis outputs under `sepsis/`. The plan is intentionally procedural and data-bound. Wherever the manuscript needs a number, term name, rank, or parameter, the plan requires you to copy it from the actual output files (for example, `*_term_matches.csv`, `*_qc.csv`, and `*_run_info.json`) rather than inferring it.

The section plans are written in the recommended drafting order, which minimizes the risk of accidentally leaking results into Methods/Introduction and ensures that every statement is supported by an output artifact.

## Methods plan

Write Methods using `docs/manuscript_methods_writing_plan.md`. This plan describes how to extract only reproducibility-relevant details from `sepsis/shared_doc/SHARED_*.METHOD.md`, per-library `*_*.PARAMS.md`, and `*_run_info.json`, while avoiding any result statements.

## Results plan

Write Results using `docs/manuscript_results_writing_plan.md`. This plan is strict about reporting only observed outputs (no interpretation), selecting clinically relevant libraries, and selecting figures and tables by reading the actual CSV/JSON outputs and confirming figure availability in each per-library `figures/` folder. The plan also defines how to decide between single-figure versus multi-panel figures (Figure 1a, 1b, etc.) and how to write publication-standard figure and table legends.

## Introduction plan

Write Introduction using `docs/manuscript_introduction_writing_plan.md`. This plan focuses on the clinical framing of sepsis and the motivation for embedding-based “soft matching,” without including any results.

## Discussion plan

Write Discussion using `docs/manuscript_discussion_writing_plan.md`. This plan explains how to interpret results for clinicians, connect to clinical context and prior literature, and state limitations, while keeping the interpretation grounded in the reported Results and avoiding any new, unreported findings.

## Abstract plan

Write Abstract using `docs/manuscript_abstract_writing_plan.md`. This plan describes how to assemble a concise abstract at the end, using only statements that already appear in Methods/Results/Discussion.


