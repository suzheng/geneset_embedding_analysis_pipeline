from __future__ import annotations

import os
import sys
from pathlib import Path

# Add project root to Python path to import src modules
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from src.main import main as run_main
from scripts.generate_manuscript_prompt import generate_manuscript_prompt


def main() -> None:
    """
    Example script demonstrating how to run the geneset embedding analysis
    for a library group using the updated main.py
    """
    
    # Set up paths
    embeddings_path = str(project_root / "data" / "GeneRAIN-vec.200d.txt")
    
    # Example input gene list (pick a few common ones likely present)
    example_genes = "TLR1,TLR2,LBP,RAGE,NOD2,MBL2,TNFA,LTA,IL1B,IL6,IL8,IL10,IFNG,CXCL1,CXCL12,MIF,FCGR2A,PAI-1,ACE,LCE4A,TAGAP,VDR,MIR608"
    
    # Prepare command-line arguments for main.py
    # Run analysis for the "Gene_Perturbations" group as an example
    sys.argv = [
        "run_real_data_example.py",
        "--embeddings-path", embeddings_path,
        "--entity-id", "sepsis",
        "--associated-item", "sepsis",
        "--input-genes", example_genes,
        "--library-group", "all-but-disease",  # Try: Gene_Ontology, Pathways, Drug_and_Chemical_Perturbations, etc.
        "--dims-filter", "quantile",
        "--q-low", "0.3",
        "--top-k", "20",
        "--no-clustering",
    ]
    
    print("=" * 80)
    print("Running geneset embedding analysis with library group: Pathways")
    print("=" * 80)
    print()
    
    # Call the main function from main.py
    run_main()
    
    # After run: generate consolidated manuscript prompt in the entity-id folder
    def _get_arg_value(flag: str) -> str | None:
        if flag in sys.argv:
            i = sys.argv.index(flag)
            if i + 1 < len(sys.argv):
                return sys.argv[i + 1]
        return None

    entity_id = _get_arg_value("--entity-id") or "run"
    output_folder = _get_arg_value("--output-folder")
    workspace = os.getenv("PROJECT_WORKSPACE_VOLUME_PATH") or os.getcwd()
    if output_folder:
        entity_dir = Path(workspace) / output_folder / entity_id
    else:
        entity_dir = Path(workspace) / entity_id

    prompt_path = generate_manuscript_prompt(str(entity_dir))
    print()
    print("Generated manuscript prompt:")
    print(prompt_path)
    
    print()
    print("=" * 80)
    print("Example completed!")
    print("=" * 80)
    print()
    print("To run other groups, modify the --library-group argument:")
    print("  - Gene_Ontology (3 libraries)")
    print("  - Pathways (2 libraries)")
    print("  - Drug_and_Chemical_Perturbations (7 libraries)")
    print("  - Gene_Perturbations (3 libraries)")
    print("  - Disease_and_Phenotype_Associations (6 libraries)")
    print("  - Rare_Disease_Associations (2 libraries)")
    print("  - Tissue_and_Cell_Type_Markers (3 libraries)")
    print("  - Cellular_Compartments (2 libraries)")
    print("  - Protein_Domains_and_Complexes (3 libraries)")
    print("  - Gene_Regulation (4 libraries)")
    print("  - Miscellaneous (7 libraries)")
    print("  - all (run all groups)")


if __name__ == "__main__":
    main()


