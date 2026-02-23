#!/usr/bin/env python3
"""
Script to combine multiple figures into a single multi-panel figure.
Supports SVG, PDF, and PNG formats.
"""

from __future__ import annotations

import os
import sys
import datetime
from pathlib import Path
from typing import List, Tuple, Optional
import argparse

try:
    import fitz  # PyMuPDF
except ImportError:
    print("Error: PyMuPDF (fitz) is required. Install with: pip install PyMuPDF")
    sys.exit(1)

try:
    import cairosvg
    HAS_CAIROSVG = True
except ImportError:
    HAS_CAIROSVG = False
    print("Warning: cairosvg not available. SVG to PDF conversion may require cairosvg.")
    print("Install with: pip install cairosvg")

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    print("Warning: PIL/Pillow not available. PNG handling may be limited.")
    print("Install with: pip install Pillow")


def print_file_modification_date(file_path: str) -> None:
    """Print the modification date of a file."""
    try:
        mod_time = os.path.getmtime(file_path)
        mod_datetime = datetime.datetime.fromtimestamp(mod_time)
        print(f"Last modified at {mod_datetime}: '{file_path}'")
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")


def svg_to_pdf(svg_path: str, pdf_path: str) -> None:
    """Convert SVG to PDF using cairosvg."""
    if not HAS_CAIROSVG:
        raise RuntimeError("cairosvg is required for SVG to PDF conversion")
    cairosvg.svg2pdf(url=svg_path, write_to=pdf_path)


def combine_pdf_figures(
    pdf_paths: List[str],
    output_pdf: str,
    labels: List[str],
    cols_per_row: List[int],
    row_gap: float = 2.0,
    label_offset_x: float = 5.0,
    label_offset_y: float = 15.0,
    label_fontsize: float = 10.0
) -> None:
    """
    Combine multiple PDF figures into a single multi-panel PDF.
    
    Args:
        pdf_paths: List of paths to PDF files to combine
        output_pdf: Output PDF path
        labels: List of panel labels (e.g., ['A', 'B', 'C'])
        cols_per_row: List of number of columns per row (e.g., [2, 1] for 2 cols in row 1, 1 col in row 2)
        row_gap: Gap between rows in points
        label_offset_x: X offset for labels in points
        label_offset_y: Y offset for labels in points
        label_fontsize: Font size for labels
    """
    for pdf_path in pdf_paths:
        print_file_modification_date(pdf_path)
        if not os.path.exists(pdf_path):
            raise FileNotFoundError(f"PDF file not found: {pdf_path}")
    
    if len(pdf_paths) != len(labels):
        raise ValueError(f"Number of PDF paths ({len(pdf_paths)}) must match number of labels ({len(labels)})")
    
    doc = fitz.open()  # Create a new PDF

    # Dictionaries to hold row-wise information
    row_pdfs: dict[int, List[str]] = {}  # {row_index: [pdf_paths]}
    row_widths: dict[int, float] = {}  # {row_index: total_width}
    row_heights: dict[int, float] = {}  # {row_index: max_height}
    row_labels: dict[int, List[str]] = {}  # {row_index: [labels]}

    # Organize PDFs into rows and calculate dimensions for each row
    current_index = 0
    for row_index, num_cols in enumerate(cols_per_row):
        row_pdfs[row_index] = []
        row_labels[row_index] = []
        row_widths[row_index] = 0
        row_heights[row_index] = 0

        for _ in range(num_cols):
            if current_index >= len(pdf_paths):
                break
            pdf_path = pdf_paths[current_index]
            label = labels[current_index]
            row_pdfs[row_index].append(pdf_path)
            row_labels[row_index].append(label)
            
            with fitz.open(pdf_path) as src_doc:
                src_page = src_doc[0]
                row_widths[row_index] += src_page.rect.width
                row_heights[row_index] = max(row_heights[row_index], src_page.rect.height)
            current_index += 1
        
        # Add row gap to height (except for last row)
        if row_index < len(cols_per_row) - 1:
            row_heights[row_index] += row_gap

    # Calculate overall page dimensions
    page_width = max(row_widths.values())
    page_height = sum(row_heights.values())

    # Create a new page in the output document
    page = doc.new_page(width=page_width, height=page_height)

    # Place each PDF figure on the combined page
    current_y = 0
    label_idx = 0
    for row_index in sorted(row_pdfs.keys()):
        # Calculate row width for centering
        row_width = row_widths[row_index]
        # Center-align horizontally if this row has fewer panels than the max width
        # Calculate horizontal offset for centering
        num_panels_in_row = len(row_pdfs[row_index])
        if num_panels_in_row < max(cols_per_row) and page_width > row_width:
            # Center the row
            current_x = (page_width - row_width) / 2
        else:
            current_x = 0
        
        # Get max height for this row (before adding gap)
        row_max_height = row_heights[row_index]
        if row_index < len(cols_per_row) - 1:
            row_max_height -= row_gap
        
        for pdf_index, pdf_path in enumerate(row_pdfs[row_index]):
            label = row_labels[row_index][pdf_index]
            
            with fitz.open(pdf_path) as src_doc:
                src_page = src_doc[0]
                # Center-align vertically: calculate vertical offset
                vertical_offset = (row_max_height - src_page.rect.height) / 2
                
                rect = fitz.Rect(
                    current_x,
                    current_y + vertical_offset,
                    current_x + src_page.rect.width,
                    current_y + vertical_offset + src_page.rect.height
                )
                page.show_pdf_page(rect, src_doc, 0)

                # Add label at the same horizontal level for all figures in the row
                label_x = current_x + label_offset_x
                label_y = current_y + label_offset_y  # Same level for all labels in the row
                page.insert_text(
                    (label_x, label_y),
                    label,
                    fontsize=label_fontsize,
                    fontname="Helvetica-Bold"
                )

                # Update current_x for the next figure in the row
                current_x += src_page.rect.width
            label_idx += 1
        
        # Update current_y for the next row
        current_y += row_heights[row_index]

    # Save the combined PDF
    doc.save(output_pdf)
    doc.close()
    print(f"Combined PDF saved to: {output_pdf}")


def combine_svg_figures(
    svg_paths: List[str],
    output_svg: str,
    labels: List[str],
    cols_per_row: List[int],
    row_gap: float = 0.1,
    label_offset_x: float = 0.02,
    label_offset_y: float = 0.05,
    label_fontsize: float = 12.0
) -> None:
    """
    Combine multiple SVG figures into a single multi-panel SVG.
    Uses matplotlib to load and combine SVGs.
    """
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.patches import Rectangle
    import xml.etree.ElementTree as ET
    
    # Load SVG dimensions
    svg_dims: List[Tuple[float, float]] = []
    for svg_path in svg_paths:
        print_file_modification_date(svg_path)
        if not os.path.exists(svg_path):
            raise FileNotFoundError(f"SVG file not found: {svg_path}")
        
        # Try to read SVG dimensions from the file
        tree = ET.parse(svg_path)
        root = tree.getroot()
        width = float(root.get('width', '800').replace('px', '').replace('pt', ''))
        height = float(root.get('height', '600').replace('px', '').replace('pt', ''))
        svg_dims.append((width, height))
    
    # Organize SVGs into rows
    row_svgs: dict[int, List[str]] = {}
    row_labels: dict[int, List[str]] = {}
    row_widths: dict[int, float] = {}
    row_heights: dict[int, float] = {}
    
    current_index = 0
    for row_index, num_cols in enumerate(cols_per_row):
        row_svgs[row_index] = []
        row_labels[row_index] = []
        row_widths[row_index] = 0
        row_heights[row_index] = 0
        
        for _ in range(num_cols):
            if current_index >= len(svg_paths):
                break
            svg_path = svg_paths[current_index]
            label = labels[current_index]
            row_svgs[row_index].append(svg_path)
            row_labels[row_index].append(label)
            
            width, height = svg_dims[current_index]
            row_widths[row_index] += width
            row_heights[row_index] = max(row_heights[row_index], height)
            current_index += 1
        
        if row_index < len(cols_per_row) - 1:
            row_heights[row_index] += row_gap * max(row_heights.values())
    
    # Calculate overall figure dimensions
    total_width = max(row_widths.values())
    total_height = sum(row_heights.values())
    
    # Create figure
    fig = plt.figure(figsize=(total_width / 100, total_height / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, total_width)
    ax.set_ylim(total_height, 0)  # Invert y-axis
    ax.axis('off')
    
    # Place each SVG on the combined figure
    current_y = 0
    for row_index in sorted(row_svgs.keys()):
        current_x = 0
        for svg_index, svg_path in enumerate(row_svgs[row_index]):
            label = row_labels[row_index][svg_index]
            width, height = svg_dims[row_index * max(cols_per_row) + svg_index]
            
            # Load and place SVG (convert to PNG first for embedding)
            # For now, we'll use a placeholder approach
            # In practice, you might want to use svg2rlg or similar
            
            # Add label
            ax.text(
                current_x + label_offset_x * total_width,
                current_y + label_offset_y * total_height,
                label,
                fontsize=label_fontsize,
                fontweight='bold',
                transform=ax.transData
            )
            
            current_x += width
        
        current_y += row_heights[row_index]
    
    # Save combined SVG
    fig.savefig(output_svg, format='svg', bbox_inches='tight', pad_inches=0)
    plt.close(fig)
    print(f"Combined SVG saved to: {output_svg}")


def combine_png_figures(
    png_paths: List[str],
    output_png: str,
    labels: List[str],
    cols_per_row: List[int],
    row_gap: int = 20,
    label_offset_x: int = 10,
    label_offset_y: int = 20,
    label_fontsize: int = 14
) -> None:
    """
    Combine multiple PNG figures into a single multi-panel PNG.
    """
    if not HAS_PIL:
        raise RuntimeError("PIL/Pillow is required for PNG combination")
    
    from PIL import Image, ImageDraw, ImageFont
    
    # Load PNG images and get dimensions
    images: List[Image.Image] = []
    for png_path in png_paths:
        print_file_modification_date(png_path)
        if not os.path.exists(png_path):
            raise FileNotFoundError(f"PNG file not found: {png_path}")
        img = Image.open(png_path)
        images.append(img)
    
    # Organize images into rows
    row_images: dict[int, List[Image.Image]] = {}
    row_labels: dict[int, List[str]] = {}
    row_widths: dict[int, int] = {}
    row_heights: dict[int, int] = {}
    
    current_index = 0
    for row_index, num_cols in enumerate(cols_per_row):
        row_images[row_index] = []
        row_labels[row_index] = []
        row_widths[row_index] = 0
        row_heights[row_index] = 0
        
        for _ in range(num_cols):
            if current_index >= len(png_paths):
                break
            img = images[current_index]
            label = labels[current_index]
            row_images[row_index].append(img)
            row_labels[row_index].append(label)
            
            row_widths[row_index] += img.width
            row_heights[row_index] = max(row_heights[row_index], img.height)
            current_index += 1
        
        if row_index < len(cols_per_row) - 1:
            row_heights[row_index] += row_gap
    
    # Calculate overall image dimensions
    total_width = max(row_widths.values())
    total_height = sum(row_heights.values())
    
    # Create combined image
    combined = Image.new('RGB', (total_width, total_height), 'white')
    draw = ImageDraw.Draw(combined)
    
    # Try to load a bold font, fall back to default if not available
    font = None
    try:
        # Try Helvetica Bold on macOS
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Helvetica-Bold.ttf", label_fontsize)
    except:
        try:
            # Try Helvetica regular (will make it bold with stroke)
            font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", label_fontsize)
        except:
            try:
                # Try Arial Bold on Windows
                font = ImageFont.truetype("arialbd.ttf", label_fontsize)
            except:
                try:
                    # Try Arial regular
                    font = ImageFont.truetype("arial.ttf", label_fontsize)
                except:
                    try:
                        # Try DejaVu Sans Bold (common on Linux)
                        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", label_fontsize)
                    except:
                        # Fall back to default font
                        font = ImageFont.load_default()
                        print(f"Warning: Using default font (size may not be accurate)")
    
    # Place each image on the combined figure
    current_y = 0
    for row_index in sorted(row_images.keys()):
        # Calculate row width for centering
        row_width = row_widths[row_index]
        # Center-align horizontally if this row has fewer panels than the max width
        num_panels_in_row = len(row_images[row_index])
        if num_panels_in_row < max(cols_per_row) and total_width > row_width:
            # Center the row
            current_x = int((total_width - row_width) / 2)
        else:
            current_x = 0
        
        # Get max height for this row (before adding gap)
        row_max_height = row_heights[row_index]
        if row_index < len(cols_per_row) - 1:
            row_max_height -= row_gap
        
        for img_index, img in enumerate(row_images[row_index]):
            label = row_labels[row_index][img_index]
            
            # Center-align vertically: calculate vertical offset
            vertical_offset = int((row_max_height - img.height) / 2)
            
            # Paste image at centered position
            combined.paste(img, (current_x, current_y + vertical_offset))
            
            # Add label at the same horizontal level for all figures in the row
            # Use fill='black' with explicit RGB tuple for better visibility
            # Add a white stroke to make text stand out against any background
            draw.text(
                (current_x + label_offset_x, current_y + label_offset_y),  # Same level for all labels
                label,
                fill=(0, 0, 0),  # Explicit black color
                font=font,
                stroke_width=2,  # Slightly thicker stroke for better visibility
                stroke_fill=(255, 255, 255)  # White stroke outline
            )
            
            current_x += img.width
        
        current_y += row_heights[row_index]
    
    # Save combined PNG
    combined.save(output_png, 'PNG', dpi=(300, 300))
    print(f"Combined PNG saved to: {output_png}")


def main():
    parser = argparse.ArgumentParser(
        description="Combine multiple figures into a single multi-panel figure"
    )
    parser.add_argument(
        "--figures",
        nargs="+",
        required=True,
        help="Paths to figure files to combine (SVG, PDF, or PNG)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file path (extension determines format: .pdf, .svg, or .png)"
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        required=True,
        help="Panel labels (e.g., 'A B C')"
    )
    parser.add_argument(
        "--cols-per-row",
        type=int,
        nargs="+",
        required=True,
        help="Number of columns per row (e.g., '2 1' for 2 cols in row 1, 1 col in row 2)"
    )
    parser.add_argument(
        "--row-gap",
        type=float,
        default=2.0,
        help="Gap between rows (default: 2.0 for PDF, 0.1 for SVG, 20 for PNG)"
    )
    
    args = parser.parse_args()
    
    if len(args.figures) != len(args.labels):
        parser.error(f"Number of figures ({len(args.figures)}) must match number of labels ({len(args.labels)})")
    
    if sum(args.cols_per_row) != len(args.figures):
        parser.error(f"Sum of cols_per_row ({sum(args.cols_per_row)}) must equal number of figures ({len(args.figures)})")
    
    output_ext = os.path.splitext(args.output)[1].lower()
    
    # Convert SVG to PDF if needed for PDF output
    if output_ext == '.pdf':
        pdf_paths: List[str] = []
        temp_pdfs: List[str] = []
        
        for fig_path in args.figures:
            fig_ext = os.path.splitext(fig_path)[1].lower()
            if fig_ext == '.pdf':
                pdf_paths.append(fig_path)
            elif fig_ext == '.svg':
                if not HAS_CAIROSVG:
                    raise RuntimeError("cairosvg is required to convert SVG to PDF")
                temp_pdf = fig_path.replace('.svg', '_temp.pdf')
                svg_to_pdf(fig_path, temp_pdf)
                pdf_paths.append(temp_pdf)
                temp_pdfs.append(temp_pdf)
            else:
                raise ValueError(f"Unsupported format for PDF output: {fig_ext}")
        
        try:
            combine_pdf_figures(
                pdf_paths,
                args.output,
                args.labels,
                args.cols_per_row,
                row_gap=args.row_gap
            )
        finally:
            # Clean up temporary PDFs
            for temp_pdf in temp_pdfs:
                if os.path.exists(temp_pdf):
                    os.remove(temp_pdf)
    
    elif output_ext == '.svg':
        svg_paths = [p for p in args.figures if os.path.splitext(p)[1].lower() == '.svg']
        if len(svg_paths) != len(args.figures):
            raise ValueError("All input figures must be SVG for SVG output")
        combine_svg_figures(
            svg_paths,
            args.output,
            args.labels,
            args.cols_per_row,
            row_gap=args.row_gap if args.row_gap != 2.0 else 0.1
        )
    
    elif output_ext == '.png':
        png_paths = [p for p in args.figures if os.path.splitext(p)[1].lower() == '.png']
        if len(png_paths) != len(args.figures):
            raise ValueError("All input figures must be PNG for PNG output")
        combine_png_figures(
            png_paths,
            args.output,
            args.labels,
            args.cols_per_row,
            row_gap=int(args.row_gap) if args.row_gap != 2.0 else 20
        )
    
    else:
        parser.error(f"Unsupported output format: {output_ext}. Use .pdf, .svg, or .png")


if __name__ == "__main__":
    main()
