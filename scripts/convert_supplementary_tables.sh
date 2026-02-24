#!/bin/bash

# Script to convert CSV files to Excel format for supplementary tables

# Set the project directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# Try to use venv Python, fall back to system Python
if [ -f "venv/bin/activate" ]; then
    source venv/bin/activate
    PYTHON_CMD="python"
else
    PYTHON_CMD="python3"
fi

# Use system Python3 if venv Python doesn't work
if ! $PYTHON_CMD -c "import openpyxl" 2>/dev/null; then
    PYTHON_CMD="/usr/local/bin/python3"
    if ! $PYTHON_CMD -c "import openpyxl" 2>/dev/null; then
        echo "Error: openpyxl is required. Installing..."
        $PYTHON_CMD -m pip install openpyxl --quiet
    fi
fi

# Run the conversion script
$PYTHON_CMD scripts/convert_csv_to_excel.py
