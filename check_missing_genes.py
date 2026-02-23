#!/usr/bin/env python3

import os
from pathlib import Path

def load_gene_symbols(embeddings_path: str) -> set[str]:
    """Load all gene symbols from the embeddings file."""
    gene_symbols = set()
    with open(embeddings_path, 'r', encoding='utf-8') as f:
        # Skip header line
        next(f)
        for line in f:
            if line.strip():
                parts = line.strip().split()
                if parts:
                    symbol = parts[0].strip().upper()
                    gene_symbols.add(symbol)
    return gene_symbols

def main():
    # Example genes from the script (updated with correct aliases)
    example_genes = [
        "TLR1", "TLR2", "LBP", "AGER", "NOD2", "MBL2", "TNF", "LTA", "IL1B", "IL6", "CXCL8", "IL10",
        "IFNG", "CXCL1", "CXCL12", "MIF", "FCGR2A", "SERPINE1", "ACE", "LCE4A", "TAGAP", "VDR", "MIR608"
    ]

    # Path to embeddings file
    embeddings_path = "data/GeneRAIN-vec.200d.txt"

    if not os.path.exists(embeddings_path):
        print(f"Error: Embeddings file not found at {embeddings_path}")
        return

    # Load all available gene symbols
    available_symbols = load_gene_symbols(embeddings_path)

    print(f"Total genes in embeddings: {len(available_symbols)}")
    print(f"Checking {len(example_genes)} example genes...")
    print()

    found = []
    missing = []

    for gene in example_genes:
        normalized_gene = gene.upper()
        if normalized_gene in available_symbols:
            found.append((gene, normalized_gene))
        else:
            missing.append((gene, normalized_gene))

    print(f"Found: {len(found)}/{len(example_genes)}")
    print(f"Missing: {len(missing)}/{len(example_genes)}")
    print()

    if found:
        print("FOUND genes:")
        for gene, normalized in found:
            print(f"  ✓ {gene} (normalized: {normalized})")

    print()

    if missing:
        print("MISSING genes:")
        for gene, normalized in missing:
            print(f"  ✗ {gene} (normalized: {normalized})")

        print()
        print("Checking for potential aliases...")

        # Check for some common aliases (already applied corrections)
        potential_aliases = {
            "SERPINE1": ["PAI1", "PLANH1"],  # PAI-1 corrected to SERPINE1
            "MIR608": ["MIRN608", "MIMAT0003276", "MI0003621"],  # miRNA naming variations
            "TNF": ["TNFA", "TNF-A"],  # TNFA corrected to TNF
            "LTA": ["TNFB", "TNFSF1"],
            "CXCL12": ["SDF1", "SDF1A"],
            "FCGR2A": ["FCGR2", "CD32"],
            "ACE": ["ACE1", "DCP1"],
            "IFNG": ["IFN-G"],
            "IL1B": ["IL1-B", "IL1F2"],
            "CXCL1": ["GRO1", "GRO"],
            "AGER": ["RAGE"],  # RAGE corrected to AGER
            "CXCL8": ["IL8"],  # IL8 corrected to CXCL8
            "LCE4A": ["LCE4", "LEP8", "SPRL4A"],  # Late cornified envelope protein 4A
        }

        for gene, normalized in missing:
            if gene in potential_aliases:
                print(f"\nChecking aliases for {gene}:")
                aliases_found = []
                aliases_missing = []
                for alias in potential_aliases[gene]:
                    alias_norm = alias.upper()
                    if alias_norm in available_symbols:
                        aliases_found.append(alias)
                    else:
                        aliases_missing.append(alias)

                if aliases_found:
                    print(f"  ✓ Found aliases: {', '.join(aliases_found)}")
                if aliases_missing:
                    print(f"  ✗ Missing aliases: {', '.join(aliases_missing)}")

if __name__ == "__main__":
    main()