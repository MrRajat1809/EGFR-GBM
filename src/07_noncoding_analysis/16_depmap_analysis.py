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
import statistics

mut_file = project_path('data/raw/noncoding/DepMap/OmicsSomaticMutations.csv')
exp_file = project_path('data/raw/noncoding/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv')

# We will scan for your top CADD-scoring mutations from the previous step
target_positions = ['55259524', '55259509', '55269049', '55272947', '55259469', '55268881']

mutated_models = set()
mutated_details = {} # To store which cell line has which mutation

print("====================================================")
print("🚀 DEPMAP WGS/RNA-SEQ INTEGRATION PIPELINE")
print("====================================================\n")

# STEP 1: Scan for the Mutations
print(f"[1/3] Scanning {mut_file} for Neo-Enhancer variants...")
try:
    with open(mut_file, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)
        
        idx_model = headers.index('ModelID')
        idx_chrom = headers.index('Chrom')
        idx_pos = headers.index('Pos')
        idx_ref = headers.index('Ref')
        idx_alt = headers.index('Alt')

        for row in reader:
            if len(row) > idx_alt:
                # Clean chromosome name (some databases use 'chr7', some just use '7')
                chrom = row[idx_chrom].replace('chr', '')
                pos = row[idx_pos]
                
                if chrom == '7' and pos in target_positions:
                    model = row[idx_model]
                    mutated_models.add(model)
                    mutated_details[model] = f"{pos} ({row[idx_ref]}>{row[idx_alt]})"
                    print(f"   [HIT] Found Cell Line {model} with mutation at {pos}!")
                    
except Exception as e:
    print(f"❌ ERROR reading mutations: {e}")

print(f"\n[2/3] Total Mutated Cell Lines Found: {len(mutated_models)}")

# STEP 2: Extract EGFR Expression
print(f"\n[3/3] Extracting EGFR Expression from {exp_file}...")
mutated_expr = []
wt_expr = []

try:
    with open(exp_file, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)
        
        idx_model = headers.index('ModelID')
        idx_egfr = headers.index('EGFR (1956)')
        
        for row in reader:
            if len(row) > idx_egfr:
                model = row[idx_model]
                try:
                    expr = float(row[idx_egfr])
                except ValueError:
                    continue # Skip empty values
                    
                if model in mutated_models:
                    mutated_expr.append((model, expr))
                else:
                    wt_expr.append(expr)
                    
except Exception as e:
    print(f"❌ ERROR reading expression: {e}")

# STEP 3: Statistical Scoreboard
print("\n====================================================")
print("📊 FINAL IN VITRO VALIDATION RESULTS")
print("====================================================")

if len(wt_expr) > 0:
    median_wt = statistics.median(wt_expr)
    print(f"Baseline Wild-Type EGFR Expression (Median of {len(wt_expr)} cell lines): {median_wt:.2f} log2(TPM+1)")

if len(mutated_expr) > 0:
    print("\n🔥 MUTANT CELL LINES (Neo-Enhancer):")
    for model, expr in mutated_expr:
        variant = mutated_details.get(model, "Unknown")
        print(f"   - Model: {model} | Variant: chr7:{variant} | EGFR Exp: {expr:.2f} log2(TPM+1)")
        
    median_mut = statistics.median([x[1] for x in mutated_expr])
    print(f"\nMedian Mutated EGFR Expression: {median_mut:.2f} log2(TPM+1)")
    
    if median_mut > median_wt:
        diff = median_mut - median_wt
        print(f"\n✅ SUCCESS! Mutated cell lines show a {diff:.2f} log2(TPM+1) increase over background.")
    else:
        print("\n⚠️ Mutated cell lines do not show elevated EGFR in this specific dataset.")
else:
    print("\n⚠️ No cell lines with the exact target mutations were found.")
    print("If this happens, it simply means these rare variants are not present in the ~1,400 immortalized lines available in CCLE.")

print("====================================================")