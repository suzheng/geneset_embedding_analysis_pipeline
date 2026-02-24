#!/usr/bin/env python3
"""
Script to convert CSV files to Excel format for supplementary tables.
"""

import csv
from pathlib import Path
import sys

try:
    from openpyxl import Workbook
    from openpyxl.styles import Font
except ImportError:
    print("Error: openpyxl is required. Install it with: pip install openpyxl")
    sys.exit(1)

def convert_csv_to_excel(csv_path: Path, excel_path: Path) -> None:
    """Convert a CSV file to Excel format."""
    print(f"Converting {csv_path.name} to {excel_path.name}...")
    
    # Create workbook and worksheet
    wb = Workbook()
    ws = wb.active
    ws.title = "Data"
    
    # Read CSV file and write to Excel
    with open(csv_path, 'r', encoding='utf-8') as csvfile:
        csv_reader = csv.reader(csvfile)
        for row_idx, row in enumerate(csv_reader, start=1):
            for col_idx, value in enumerate(row, start=1):
                cell = ws.cell(row=row_idx, column=col_idx, value=value)
                # Make header row bold
                if row_idx == 1:
                    cell.font = Font(bold=True)
    
    # Save Excel file
    excel_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(excel_path)
    
    print(f"✓ Created {excel_path}")

def main():
    """Main function to convert all supplementary table CSV files to Excel."""
    # Get project root directory
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    
    # Define CSV to Excel mappings
    csv_files = [
        {
            'csv': project_root / 'sepsis' / 'sepsis_CellMarker_2024' / 'sepsis_CellMarker_2024_geneset_embed_sepsis_CellMarker_2024_term_matches.csv',
            'excel_name': 'Supplementary_Table_1_CellMarker_2024.xlsx'
        },
        {
            'csv': project_root / 'sepsis' / 'sepsis_Reactome_Pathways_2024' / 'sepsis_Reactome_Pathways_2024_geneset_embed_sepsis_Reactome_Pathways_2024_term_matches.csv',
            'excel_name': 'Supplementary_Table_2_Reactome_Pathways_2024.xlsx'
        },
        {
            'csv': project_root / 'sepsis' / 'sepsis_DSigDB' / 'sepsis_DSigDB_geneset_embed_sepsis_DSigDB_term_matches.csv',
            'excel_name': 'Supplementary_Table_3_DSigDB.xlsx'
        },
        {
            'csv': project_root / 'sepsis' / 'sepsis_RummaGEO_DrugPerturbations_2025' / 'sepsis_RummaGEO_DrugPerturbations_2025_geneset_embed_sepsis_RummaGEO_DrugPerturbations_2025_term_matches.csv',
            'excel_name': 'Supplementary_Table_4_RummaGEO_DrugPerturbations_2025.xlsx'
        }
    ]
    
    # Create output directory
    output_dir = project_root / 'docs' / 'manuscript_v2' / 'supplementary_tables'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Converting CSV files to Excel format...")
    print("=" * 60)
    
    # Convert each CSV file
    for file_info in csv_files:
        csv_path = file_info['csv']
        excel_path = output_dir / file_info['excel_name']
        
        if not csv_path.exists():
            print(f"⚠ Warning: CSV file not found: {csv_path}")
            continue
        
        try:
            convert_csv_to_excel(csv_path, excel_path)
        except Exception as e:
            print(f"✗ Error converting {csv_path.name}: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    
    print("=" * 60)
    print(f"✓ All Excel files created in: {output_dir}")
    print(f"\nExcel files:")
    for file_info in csv_files:
        excel_path = output_dir / file_info['excel_name']
        if excel_path.exists():
            print(f"  - {excel_path.name}")

if __name__ == '__main__':
    main()
