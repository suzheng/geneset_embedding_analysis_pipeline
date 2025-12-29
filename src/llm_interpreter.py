import os
import sys
from typing import List, Optional, Dict, Any

# Add the project root to Python path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

# Load environment variables from .env if present
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

from openai import OpenAI


def interpret_top_matches_with_llm(
    top_matches: List[Dict[str, Any]],
    library_name: str,
    background_text: str,
    associated_item: str,
    entity_id: str,
    output_dir: str,
    query_genes: Optional[List[str]] = None,
) -> Optional[str]:
    """
    Interpret top matched terms using LLM with web search capability.
    
    Args:
        top_matches: List of dicts with 'term', 'distance' keys
        library_name: Name of the reference library
        background_text: Research background context
        associated_item: What the input genes are associated with (e.g., disease, phenotype, event)
        query_genes: Optional list of input query gene symbols to include in prompt
        entity_id: Unique identifier for the run
        output_dir: Directory to save interpretation results
        
    Returns:
        Interpretation text if successful, None otherwise
    """
    # Get model name from environment variable
    model_name = os.getenv("ACADEMIC_WEB_SEARCH_LLM_FOR_TOOLS")
    if not model_name:
        model_name = "gpt-5-mini"
    
    # Get OpenAI API key
    openai_api_key = os.getenv("OPENAI_API_KEY")
    if not openai_api_key:
        print("WARNING: OPENAI_API_KEY environment variable not set. Skipping LLM interpretation.")
        return None
    
    # Initialize OpenAI client
    client = OpenAI(api_key=openai_api_key)
    
    # Format top matches
    top_matches_text = "=== TOP 5 MATCHES ===\n"
    for i, match in enumerate(top_matches[:5], 1):
        distance = match.get('distance', 'N/A')
        term = match.get('term', 'unknown')
        top_matches_text += f"{i}. {term} (distance: {distance})\n"
    
    # Format query gene list (concise)
    genes_section = ""
    if query_genes:
        max_show = 100
        shown = [str(g) for g in list(query_genes)[:max_show]]
        genes_line = ", ".join(shown)
        if len(query_genes) > max_show:
            genes_line += f", ...( +{len(query_genes) - max_show} more )"
        genes_section = f"Query genes (n={len(query_genes)}): {genes_line}\n\n"

    # Build the query (concise, general-purpose)
    query = f"""
I am doing an analysis to use advanced AI and gene embedding to improve the interpretation of the genes that are associated with {associated_item}. 
Below is the analysis result for one library. Please read each found term in the library, and analyze and do literature review if necessary to interpret if the term is a new finding, 
or if it is consistent with findings from previous study if there is any. You may want to search existing study. 
Keep in mind, when you search existing study, make sure you only reference peer-reviewed academic publication. and you have to provide the references.

Output as plain text with inline citations like [URL of peer review publication] at the end of relevant sentences. 
Do not include any other details in the references—just the URL. Keep the output in plain text for easy parsing.

Keep in mind, you have to output [URL of peer review publication]!

Interpret the top matched terms for the query gene list that are associated with {associated_item}.
- For each term, state whether it is non-obvious/novel vs well-known/expected.
- If obvious/common knowledge (e.g., canonical pathways/GO expected by most researchers), do not run web search; just note consistency with established knowledge.
- Only run web search for non-obvious findings; cite peer-reviewed papers with inline [URL] at sentence ends.
- Output plain text only; be concise.

Processing library: {library_name}
Please do research and brief discussion of the finding below.
{genes_section}
{top_matches_text}
Background:
{background_text}
"""
    
    print(f"\nINFO: Starting LLM interpretation for {library_name} using model: {model_name}")
    print("      (This may take a minute as it includes web search...)")
    
    try:
        # Call OpenAI Responses API with web search
        response = client.responses.create(
            model=model_name,
            input=[
                {
                    "role": "user",
                    "content": [
                        {"type": "input_text", "text": query}
                    ],
                }
            ],
            text={"format": {"type": "text"}, "verbosity": "medium"},
            reasoning={"effort": "medium"},
            tools=[
                {
                    "type": "web_search",
                    "user_location": {"type": "approximate"},
                    "search_context_size": "medium",
                }
            ],
            store=False,
        )
        
        # Extract the response text
        result_text = None
        if hasattr(response, 'output_text') and response.output_text:
            result_text = response.output_text
        else:
            # Fallback to parsing the structured response
            output_items = getattr(response, "output", [])
            texts = []
            for item in output_items:
                for block in getattr(item, "content", []):
                    block_type = getattr(block, "type", None)
                    if block_type == "output_text":
                        text_val = getattr(block, "text", None)
                        if text_val:
                            texts.append(str(text_val))
            result_text = "\n\n".join(texts) if texts else None
        
        if result_text:
            print(f"✓ LLM interpretation completed for {library_name}")
            
            # Save the interpretation to a file
            clean_id = entity_id.replace("/", "_").replace("\\", "_")
            lib_name = library_name.replace(".tsv", "").replace(".signatures.jsonl", "")
            interpretation_file = os.path.join(
                output_dir,
                f"{clean_id}_geneset_embed_{lib_name}_llm_interpretation.txt"
            )
            
            with open(interpretation_file, 'w', encoding='utf-8') as f:
                f.write(result_text)
            
            print(f"  Saved to: {interpretation_file}")
            return result_text
        else:
            print(f"WARNING: No text extracted from LLM response for {library_name}")
            return None
            
    except Exception as e:
        print(f"WARNING: LLM interpretation failed for {library_name}: {str(e)}")
        print(f"         This might be due to model availability or API issues. Continuing without interpretation.")
        return None


def load_background_text(background_file: Optional[str] = None) -> str:
    """
    Load background text from file.
    
    Args:
        background_file: Path to background file, or None to use default
        
    Returns:
        Background text as string
    """
    if background_file is None:
        # Default to src/BACKGROUND.md
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        background_file = os.path.join(project_root, "src", "BACKGROUND.md")
    
    if not os.path.exists(background_file):
        print(f"WARNING: Background file not found at {background_file}")
        return ""
    
    with open(background_file, 'r', encoding='utf-8') as f:
        return f.read()

