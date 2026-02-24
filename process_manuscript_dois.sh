#!/bin/bash

# Script to process DOI citations in manuscript_v2 markdown files
# This script applies the DOI citation processing script to all markdown files
# in the docs/manuscript_v2 directory and outputs processed files to docs/manuscript_v2_processed/

# Set the project directory (change this if running from a different location)
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

# Activate the virtual environment
source venv/bin/activate

# Create output directory
mkdir -p docs/manuscript_v2_processed

# List of markdown files to process
MANUSCRIPT_FILES=(
    "abstract_and_introduction.md"
    "discussion.md"
    "methods.md"
    "results_cellmarker_2024.md"
    "results_dsigdb.md"
    "results_gene_set_embeddings_qc.md"
    "results_reactome_pathways_2024.md"
    "results_rna_seq_disease_geo.md"
    "results_rummageo_drugperturbations_2025.md"
)

echo "Processing DOI citations in manuscript files..."
echo "=============================================="

# Process each markdown file
for file in "${MANUSCRIPT_FILES[@]}"; do
    input_file="docs/manuscript_v2/$file"
    output_file="docs/manuscript_v2_processed/$file"
    ris_folder="docs/manuscript_v2_processed/ris/${file%.md}"

    echo "Processing $file..."
    python src/process_doi_citations.py \
        --input "$input_file" \
        --output "$output_file" \
        --ris-folder "$ris_folder" \
        --cache "docs/manuscript_v2_processed/doi_cache.json"

    echo "✓ Completed $file"
    echo ""
done

echo "All manuscript files processed successfully!"
echo "Processed files are in: docs/manuscript_v2_processed/"
echo "RIS citation files are in: docs/manuscript_v2_processed/ris/"