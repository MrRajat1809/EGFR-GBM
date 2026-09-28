# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('.')

import pandas as pd
import os

# Paths
BASE_CSV = project_path('outputs/variant_prioritization/EGFR_Arm3_CancerCorrelation_Partial.csv')
CSCAPE_TSV = project_path('data/processed/predictors/cscape/EGFR_cscape_results.tsv')
OUTPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm3_CancerCorrelation_Final.csv')

def main():
    print(f"Loading Partial Arm 3 File: {BASE_CSV}")
    df = pd.read_csv(BASE_CSV)
    
    print(f"Loading CScape results: {CSCAPE_TSV}")
    if not os.path.exists(CSCAPE_TSV):
        print("Error: CScape file not found! Did you download it from the web server?")
        return
        
    # Read the CScape TSV
    df_cscape = pd.read_csv(CSCAPE_TSV, sep='\t')
    
    # Just in case it downloaded as comma-separated
    if len(df_cscape.columns) == 1: 
        df_cscape = pd.read_csv(CSCAPE_TSV, sep=',')

    # Hardcoding the exact columns you provided
    chrom_col = '# Chromosome'
    pos_col = 'Position'
    ref_col = 'Ref. Base'
    alt_col = 'Mutant Base'
    score_col = 'Coding Score'

    # Clean the chromosome column (just in case it has 'chr7' instead of '7')
    df_cscape[chrom_col] = df_cscape[chrom_col].astype(str).str.replace('chr', '', case=False)

    # Build the exact matching key (e.g., '7_55233109_G/A')
    df_cscape['Match_Location'] = df_cscape[chrom_col] + '_' + \
                                  df_cscape[pos_col].astype(str) + '_' + \
                                  df_cscape[ref_col].astype(str) + '/' + \
                                  df_cscape[alt_col].astype(str)
    
    # Rename the score column to standard
    df_cscape.rename(columns={score_col: 'CScape_score'}, inplace=True)
    
    # CScape threshold: >= 0.5 is Oncogenic (Driver), < 0.5 is Benign (Passenger)
    df_cscape['CScape_pred'] = df_cscape['CScape_score'].apply(
        lambda x: 'Oncogenic' if pd.to_numeric(x, errors='coerce') >= 0.5 else 'Benign'
    )

    # Merge onto the main DataFrame
    df = pd.merge(df, df_cscape[['Match_Location', 'CScape_score', 'CScape_pred']], 
                  left_on='Genomic_Location', right_on='Match_Location', how='left')
    
    df.drop(columns=['Match_Location'], inplace=True)
    
    # Save the Final Arm 3 dataset
    df.to_csv(OUTPUT_CSV, index=False)
    
    print(f"\nSUCCESS! CScape merged successfully using precise column mapping.")
    print(f"Arm 3 Cancer Correlation is OFFICIALLY COMPLETE!")
    print(f"Saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()