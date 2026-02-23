#!/usr/bin/env python3
"""
Script to combine figures for RNA-Seq Disease GEO results.
Combines top_matches_example, distance_cosine_hexbin, and pca_term_space into a multi-panel figure.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from scripts.combine_figures import (
    combine_pdf_figures,
    combine_png_figures,
    svg_to_pdf,
    HAS_CAIROSVG,
    HAS_PIL
)


def main():
    """Combine RNA-Seq Disease GEO figures into multi-panel figure."""
    
    # Define paths
    figures_dir = project_root / "sepsis" / "sepsis_RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO" / "figures"
    
    # Figure files to combine (as mentioned in results_rna_seq_disease_geo.md)
    figure_files = [
        "top_matches_example.svg",
        "distance_cosine_hexbin.svg",
        "pca_term_space.svg"
    ]
    
    # Panel labels
    labels = ["A", "B", "C"]
    
    # Layout: 2 columns in row 1, 1 column in row 2 (centered)
    # This creates a nice layout: 2 panels on top, 1 panel centered on bottom
    cols_per_row = [2, 1]
    
    # Check if figures exist
    figure_paths = []
    for fig_file in figure_files:
        fig_path = figures_dir / fig_file
        if not fig_path.exists():
            print(f"Warning: Figure not found: {fig_path}")
            print("Please run the analysis first to generate figures.")
            return
        figure_paths.append(str(fig_path))
        print(f"Found: {fig_path}")
    
    # Output paths
    output_base = figures_dir / "combined_figure"
    output_pdf = str(output_base) + ".pdf"
    output_png = str(output_base) + ".png"
    output_svg = str(output_base) + ".svg"
    
    print("\n" + "=" * 80)
    print("Combining RNA-Seq Disease GEO figures")
    print("=" * 80)
    print(f"Figures: {figure_files}")
    print(f"Labels: {labels}")
    print(f"Layout: {cols_per_row[0]} columns in row 1, {cols_per_row[1]} column in row 2")
    print()
    
    # Combine as PDF
    print("Generating combined PDF...")
    try:
        # Convert SVG to PDF temporarily if needed
        pdf_paths = []
        temp_pdfs = []
        
        for svg_path in figure_paths:
            temp_pdf = svg_path.replace('.svg', '_temp.pdf')
            if HAS_CAIROSVG:
                svg_to_pdf(svg_path, temp_pdf)
                pdf_paths.append(temp_pdf)
                temp_pdfs.append(temp_pdf)
            else:
                print("Error: cairosvg is required for SVG to PDF conversion")
                print("Install with: pip install cairosvg")
                return
        
        combine_pdf_figures(
            pdf_paths,
            output_pdf,
            labels,
            cols_per_row,
            row_gap=3.0,
            label_offset_x=5.0,
            label_offset_y=15.0,
            label_fontsize=24.0
        )
        
        # Clean up temporary PDFs
        for temp_pdf in temp_pdfs:
            if os.path.exists(temp_pdf):
                os.remove(temp_pdf)
        
        print(f"✓ Combined PDF saved: {output_pdf}")
    except Exception as e:
        print(f"Error generating PDF: {e}")
        import traceback
        traceback.print_exc()
    
    # Combine as PNG
    print("\nGenerating combined PNG...")
    try:
        # Convert SVG to PNG first if needed
        if not HAS_PIL:
            print("Warning: PIL/Pillow not available. Skipping PNG generation.")
            print("Install with: pip install Pillow")
        else:
            from PIL import Image
            import cairosvg
            
            png_paths = []
            temp_pngs = []
            
            for svg_path in figure_paths:
                temp_png = svg_path.replace('.svg', '_temp.png')
                # Convert SVG to PNG using cairosvg
                cairosvg.svg2png(url=svg_path, write_to=temp_png, dpi=300)
                png_paths.append(temp_png)
                temp_pngs.append(temp_png)
            
            combine_png_figures(
                png_paths,
                output_png,
                labels,
                cols_per_row,
                row_gap=30,
                label_offset_x=25,
                label_offset_y=35,
                label_fontsize=42
            )
            
            # Clean up temporary PNGs
            for temp_png in temp_pngs:
                if os.path.exists(temp_png):
                    os.remove(temp_png)
            
            print(f"✓ Combined PNG saved: {output_png}")
    except Exception as e:
        print(f"Error generating PNG: {e}")
        import traceback
        traceback.print_exc()
    
    # Generate SVG from PDF
    print("\nGenerating combined SVG...")
    try:
        import fitz
        # Convert PDF to SVG using PyMuPDF
        pdf_doc = fitz.open(output_pdf)
        page = pdf_doc[0]
        svg_content = page.get_svg_image(matrix=fitz.Identity)
        pdf_doc.close()
        
        # Save SVG
        with open(output_svg, 'w', encoding='utf-8') as f:
            f.write(svg_content)
        
        print(f"✓ Combined SVG saved: {output_svg}")
    except Exception as e:
        print(f"Error generating SVG: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("Figure combination complete!")
    print("=" * 80)
    print(f"\nOutput files:")
    if os.path.exists(output_pdf):
        print(f"  - PDF: {output_pdf}")
    if os.path.exists(output_png):
        print(f"  - PNG: {output_png}")
    if os.path.exists(output_svg):
        print(f"  - SVG: {output_svg}")


if __name__ == "__main__":
    main()
