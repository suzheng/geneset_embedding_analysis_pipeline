# Abstract: Step-by-Step Writing Plan (Written Last, Fully Traceable)

## Scope and guardrails for this section

This Abstract plan is designed to produce a concise, journal-standard abstract that is fully consistent with the manuscript and contains no new information. The abstract must be written only after Methods, Results, and Discussion have been drafted and checked. Any numeric values, term names, or claims in the abstract must already appear in the main text and must be traceable to a table, figure, or run artifact.

## Step-by-step drafting workflow

### Step 1: Choose an abstract structure appropriate for clinical journals

Use a structured abstract if the target journal requires it, with headings such as Background, Methods, Results, and Conclusions, or use an unstructured abstract with the same logical order if the journal prefers narrative abstracts. Keep the style compact and avoid technical jargon. Do not include citations unless the target journal allows them.

### Step 2: Draft the Background in one to two sentences with no study-specific findings

Write a short clinical motivation sentence about sepsis and the need for improved functional interpretation of gene signatures. Then add a second sentence stating the general methodological idea of embedding-based gene set analysis as a soft-matching approach to connect gene sets to clinically interpretable reference knowledge. Do not include any specific pathways, drugs, or cell types here.

### Step 3: Draft the Methods in two to three sentences using only what is already in Methods

Summarize the input as a sepsis-associated gene set and describe, at a high level, computing an embedding signature and matching it to reference libraries. Mention the categories of reference libraries included in the main manuscript, such as pathways, immune cell markers, drug perturbation signatures, and disease expression signatures, but do not list library names unless they are already named in the Methods and are appropriate for the journal.

If the Methods section reports a single primary similarity metric, use the same phrasing in the abstract. Do not add implementation details such as file names or command-line arguments.

### Step 4: Draft the Results in three to five sentences by copying only the key reported outcomes

Select a small number of headline Results statements, typically one per Results subsection, and ensure each statement can be mapped directly to a figure or table in Results. Use the same terminology and ordering as the Results section. Include concrete term names and score summaries only if they are already presented in Results and if the journal’s abstract style supports that level of detail; otherwise describe results at the category level while still remaining faithful to the Results.

Do not interpret why the findings occur and do not recommend therapies. Avoid adjectives that imply interpretation. If any numeric values are included, they must be copied exactly from the reported Results tables or QC outputs.

### Step 5: Draft the Conclusions in one to two sentences, aligned to the Discussion conclusion

Write a restrained conclusion that states what the analysis provides, such as a clinically interpretable mapping of a sepsis gene set to pathways, immune cell types, drug-related signatures, and disease signatures using embedding-based similarity. Add a final sentence noting that the findings are hypothesis-generating and require external validation if this caveat is present in Discussion.

### Step 6: Perform an abstract-to-manuscript consistency audit

After drafting, check each noun phrase and each number in the abstract against the manuscript body. If any term name or value is not present in Results, remove it or add it to Results only if it is supported by the output files. Confirm that the abstract does not contain interpretation that is absent from Discussion, and confirm that the abstract does not include any Methods detail that is not described in Methods.


