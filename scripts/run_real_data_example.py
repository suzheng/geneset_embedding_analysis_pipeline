from __future__ import annotations

import os
import sys
from pathlib import Path

# Add project root to Python path to import src modules
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from src.main import main as run_main


def main() -> None:
    """
    Example script demonstrating how to run the geneset embedding analysis
    for a library group using the updated main.py
    """
    
    # Set up paths
    embeddings_path = str(project_root / "data" / "GeneRAIN-vec.200d.txt")
    
    # Example input gene list (pick a few common ones likely present)
    example_genes = "TP53,EGFR,BRCA1,BRCA2,PTEN,KRAS,APOE,MYC,ALK,MTHFR,CFTR,TNF"
    
    # Prepare command-line arguments for main.py
    # Run analysis for the "Pathways" group as an example
    sys.argv = [
        "run_real_data_example.py",
        "--embeddings-path", embeddings_path,
        "--entity-id", "example_run",
        "--input-genes", example_genes,
        "--library-group", "Pathways",  # Try: Gene_Ontology, Pathways, Drug_and_Chemical_Perturbations, etc.
        "--dims-filter", "quantile",
        "--q-low", "0.3",
        "--top-k", "25",
        "--no-clustering",
    ]
    
    print("=" * 80)
    print("Running geneset embedding analysis with library group: Pathways")
    print("=" * 80)
    print()
    
    # Call the main function from main.py
    run_main()
    
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


