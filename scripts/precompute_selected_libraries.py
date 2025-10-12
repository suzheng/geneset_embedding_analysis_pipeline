from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

# Add parent directory to path to import src modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.precompute_terms import precompute_from_file


def main() -> None:
    # Paths
    project_root = Path(__file__).resolve().parent.parent
    scripts_dir = project_root / "scripts"
    selected_libraries_file = scripts_dir / "selected_libraries.txt"
    
    # Embedding path
    embeddings_path = project_root / "data" / "GeneRAIN-vec.200d.txt"
    
    # Base directory for Enrichr TSV files
    enrichr_base_dir = Path("/Users/suzheng/noBackupData/BioBot/backend/data/tools/tool_data/database/gene_enrichment_analysis")
    
    # Output directory for precomputed signatures
    output_dir = project_root / "precomputed"
    output_dir.mkdir(exist_ok=True)
    
    # Read selected libraries
    libraries_to_process = []
    print(f"Reading library list from: {selected_libraries_file}")
    with open(selected_libraries_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            library_filename = row["Library_Filename"]
            group = row.get("Group", "Unknown")
            libraries_to_process.append((library_filename, group))
    
    print(f"Found {len(libraries_to_process)} libraries to process\n")
    
    # Process each library
    successful = 0
    failed = 0
    skipped = 0
    
    for idx, (library_filename, group) in enumerate(libraries_to_process, 1):
        print(f"[{idx}/{len(libraries_to_process)}] Processing: {library_filename} ({group})")
        
        # Construct paths
        tsv_path = enrichr_base_dir / library_filename
        output_name = library_filename.replace(".tsv", ".signatures.jsonl")
        output_path = output_dir / output_name
        
        # Check if TSV file exists
        if not tsv_path.exists():
            print(f"  ⚠️  WARNING: TSV file not found at {tsv_path}")
            skipped += 1
            continue
        
        # Check if already processed
        if output_path.exists():
            print(f"  ℹ️  Already exists: {output_path.name} (skipping)")
            skipped += 1
            continue
        
        try:
            # Precompute signatures
            print(f"  Processing terms from: {tsv_path}")
            precompute_from_file(
                str(embeddings_path),
                str(tsv_path),
                str(output_path),
                min_genes_per_term=3
            )
            print(f"  ✓ Wrote: {output_path.name}")
            successful += 1
        except Exception as e:
            print(f"  ✗ ERROR: {str(e)}")
            failed += 1
        
        print()
    
    # Summary
    print("=" * 80)
    print("PRECOMPUTATION SUMMARY")
    print("=" * 80)
    print(f"Total libraries:  {len(libraries_to_process)}")
    print(f"Successful:       {successful}")
    print(f"Skipped:          {skipped}")
    print(f"Failed:           {failed}")
    print(f"\nOutput directory: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()

