# LLM Interpretation Feature

## Overview

The LLM Interpreter module enables automatic interpretation of top matched terms from gene set embedding analysis using advanced language models with web search capability. This feature leverages the OpenAI Responses API to conduct literature reviews and provide evidence-based interpretations of findings.

## Features

- **Web Search Integration**: Automatically searches academic literature to contextualize findings
- **Peer-Reviewed References**: Extracts and citations from peer-reviewed academic publications
- **Configurable Models**: Uses environment variable to specify which LLM model to use
- **Graceful Degradation**: Continues pipeline execution even if LLM interpretation fails
- **Structured Output**: Saves interpretations to separate files for easy review

## File Location

`src/llm_interpreter.py`

## Usage

### In Pipeline

When running `geneset_embedding()` with `show_top5_matches=True`, the LLM interpretation will automatically trigger after displaying the top 5 matches:

```bash
python -m src.main \
    --embeddings-path data/GeneRAIN-vec.200d.txt \
    --entity-id test_run \
    --library-group "Gene_Regulation" \
    --show-top5-matches
```

### Direct Usage

```python
from src.llm_interpreter import interpret_top_matches_with_llm, load_background_text

# Prepare data
top_matches = [
    {"term": "Birinapant Up", "distance": "3.1006"},
    {"term": "AMG-837 Down", "distance": "3.2186"},
    # ... more matches
]

background = load_background_text()

result = interpret_top_matches_with_llm(
    top_matches=top_matches,
    library_name="LINCS_L1000_Chem_Pert_Consensus_Sigs",
    background_text=background,
    entity_id="test_run",
    output_dir="./output"
)
```

## Configuration

### Required Environment Variables

1. **OPENAI_API_KEY**: Your OpenAI API key
   ```bash
   export OPENAI_API_KEY="sk-..."
   ```

2. **ACADEMIC_WEB_SEARCH_LLM_FOR_TOOLS**: The model to use for interpretation
   ```bash
   export ACADEMIC_WEB_SEARCH_LLM_FOR_TOOLS="gpt-5-mini"
   ```

### Optional Environment Variables

- `.env` file support: If a `.env` file exists in the project root, environment variables will be loaded from it automatically

## Output

For each library processed, the LLM interpreter generates:

- **Interpretation File**: `{entity_id}_geneset_embed_{library_name}_llm_interpretation.txt`
  - Contains the full interpretation text with inline citations
  - Cites peer-reviewed publications with URLs

## API Reference

### `interpret_top_matches_with_llm()`

Interprets top matched terms using LLM with web search capability.

**Parameters:**
- `top_matches` (List[Dict[str, Any]]): List of top matched terms with 'term' and 'distance' keys
- `library_name` (str): Name of the reference library being analyzed
- `background_text` (str): Research background context for interpretation
- `entity_id` (str): Unique identifier for the analysis run
- `output_dir` (str): Directory to save interpretation results

**Returns:**
- `str | None`: Interpretation text if successful, None if skipped or failed

**Behavior:**
- Gracefully skips if environment variables not set (warning message printed)
- Catches exceptions and continues without interrupting pipeline
- Automatically saves results to file
- Prints progress messages to stdout

### `load_background_text()`

Loads background research text for context.

**Parameters:**
- `background_file` (str, optional): Path to background file. Defaults to `src/BACKGROUND.md`

**Returns:**
- `str`: Background text content, or empty string if file not found

## Integration with Pipeline

The LLM interpreter is integrated into the main pipeline at `src/pipeline.py`:

1. After top 5 matches are printed to stdout
2. Only triggers if `show_top5_matches=True`
3. Preserves existing pipeline behavior
4. Adds minimal overhead (web search takes ~30-60 seconds per library)

### Code Location

```python
# In pipeline.py, around line 970-987
if show_top5_matches and len(ranked) > 0:
    # ... print top 5 matches ...
    
    # Call LLM interpreter
    top_matches_list = [{"term": t.term, "distance": f"{dist:.4f}"} 
                        for t, dist, cos in ranked[:5]]
    background_text = load_background_text()
    interpret_top_matches_with_llm(
        top_matches=top_matches_list,
        library_name=library_name,
        background_text=background_text,
        entity_id=clean_id,
        output_dir=output_dir,
    )
```

## Example Output

When LLM interpretation runs successfully, you'll see:

```
=== TOP 5 MATCHES - LINCS_L1000_CHEM_PERT_CONSENSUS_SIGS ===
 1. Birinapant Up (distance: 3.1006)
 2. AMG-837 Down (distance: 3.2186)
 3. Ro-15-4513 Down (distance: 3.2836)
 4. MRE-269 Up (distance: 3.2952)
 5. GW-441756 Down (distance: 3.3308)

INFO: Starting LLM interpretation for LINCS_L1000_Chem_Pert_Consensus_Sigs using model: gpt-5-mini
      (This may take a minute as it includes web search...)
✓ LLM interpretation completed for LINCS_L1000_Chem_Pert_Consensus_Sigs
  Saved to: ./output/test_run_geneset_embed_LINCS_L1000_Chem_Pert_Consensus_Sigs_llm_interpretation.txt
```

## Error Handling

The interpreter implements robust error handling:

- **Missing environment variables**: Prints warning and skips interpretation
- **API errors**: Catches exceptions, prints warning, continues pipeline
- **Invalid responses**: Attempts multiple parsing strategies before failing
- **File I/O errors**: Catches and reports, continues execution

This ensures the pipeline is resilient and completes even if LLM interpretation is unavailable.

## Performance Considerations

- **Time per library**: ~30-60 seconds (due to web search)
- **Cost per library**: Depends on model and token usage
- **Best used for**: Targeted analysis of specific library groups, not full batch runs

For batch processing many libraries, consider disabling web search interpretation or running in background/parallel mode.

## Dependencies

Required packages (should already be installed):
- `openai`: OpenAI Python client
- `python-dotenv`: For `.env` file support (optional)

## Troubleshooting

### LLM interpretation is skipped

**Issue**: Seeing "WARNING: ... environment variable not set"

**Solution**: Ensure both `OPENAI_API_KEY` and `ACADEMIC_WEB_SEARCH_LLM_FOR_TOOLS` are exported:
```bash
export OPENAI_API_KEY="sk-..."
export ACADEMIC_WEB_SEARCH_LLM_FOR_TOOLS="gpt-5-mini"
```

Or add to `.env` file:
```
OPENAI_API_KEY=sk-...
ACADEMIC_WEB_SEARCH_LLM_FOR_TOOLS=gpt-5-mini
```

### Model not available

**Issue**: Error message mentioning model not found

**Solution**: Verify the model name is correct and available. Supported models include:
- `gpt-5-mini` (recommended for web search)
- Other models supporting the Responses API

Check OpenAI documentation for the latest available models.

### Partial or no citations in output

**Issue**: LLM output lacks peer-reviewed publication URLs

**Solution**: This may be due to:
- Limited web search results for the specific terms
- Model's search constraints
- Terms being too recent for comprehensive literature coverage

The interpreter will still output analysis even if citations are limited.








