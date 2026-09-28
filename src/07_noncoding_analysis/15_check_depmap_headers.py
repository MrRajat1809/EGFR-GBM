# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('outputs/noncoding/DepMap')

import csv

mut_file = project_path('data/raw/noncoding/DepMap/OmicsSomaticMutations.csv')
exp_file = project_path('data/raw/noncoding/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv')

print("====================================================")
print("🔍 INSPECTING DEPMAP DATA STRUCTURE")
print("====================================================\n")

# 1. Inspect Mutation Headers
try:
    with open(mut_file, 'r') as f:
        reader = csv.reader(f)
        mut_headers = next(reader)
        print("🧬 MUTATION FILE HEADERS (First 15):")
        print(mut_headers[:15])
except FileNotFoundError:
    print(f"❌ ERROR: Could not find {mut_file}")

print("\n----------------------------------------------------\n")

# 2. Inspect Expression Headers and Hunt for EGFR
try:
    with open(exp_file, 'r') as f:
        reader = csv.reader(f)
        exp_headers = next(reader)
        print("📊 EXPRESSION FILE HEADERS (First 5):")
        print(exp_headers[:5])
        
        print("\nHunting for 'EGFR' in the 20,000+ gene columns...")
        egfr_cols = [col for col in exp_headers if 'EGFR' in col.upper()]
        
        if egfr_cols:
            print(f"✅ FOUND IT! The exact column name to use is: '{egfr_cols[0]}'")
        else:
            print("❌ WARNING: Could not find 'EGFR' in the column names.")
            
except FileNotFoundError:
    print(f"❌ ERROR: Could not find {exp_file}")
    
print("\n====================================================")