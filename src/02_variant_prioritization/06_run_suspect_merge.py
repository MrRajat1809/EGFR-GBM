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
import sqlite3
import os

# Paths
INPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm2_DiseaseTendency_Partial.csv')
SUSPECT_DB = project_path('data/tools/suspect/suspect_package/data/suspect.db')
OUTPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm2_DiseaseTendency_Final.csv')

# EGFR Canonical UniProt ID
EGFR_UNIPROT = "P00533"

def main():
    print(f"Loading Partial Arm 2 File: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    
    # Initialize columns as object to allow mixed types, or just stick to strings
    df['SuSPect_score'] = '.'
    df['SuSPect_pred'] = '.'
    
    print(f"Connecting to SuSPect Local Database: {SUSPECT_DB}")
    if not os.path.exists(SUSPECT_DB):
        print("Error: suspect.db not found!")
        return
        
    conn = sqlite3.connect(SUSPECT_DB)
    cursor = conn.cursor()
    
    match_count = 0
    
    for i, row in df.iterrows():
        amino_acids = str(row['Amino_acids'])
        pos = str(row['Protein_position'])
        
        if '/' in amino_acids and len(amino_acids.split('/')) == 2:
            wt, mutant = amino_acids.split('/')
            
            query = f"SELECT {mutant} FROM uniprot_201303 WHERE uniprot = ? AND pos = ? AND wt = ?"
            
            try:
                cursor.execute(query, (EGFR_UNIPROT, pos, wt))
                result = cursor.fetchone()
                
                if result and result[0] is not None:
                    # Convert to float for our logic check
                    score_float = float(result[0])
                    
                    # Convert back to string to satisfy pandas strict dtype enforcement
                    df.at[i, 'SuSPect_score'] = str(score_float)
                    
                    if score_float >= 50:
                        df.at[i, 'SuSPect_pred'] = 'Pathogenic'
                    else:
                        df.at[i, 'SuSPect_pred'] = 'Benign'
                        
                    match_count += 1
            except Exception as e:
                print(f"Error querying position {pos} {wt}>{mutant}: {e}")
                
    conn.close()
    
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSUCCESS! Retrieved {match_count} SuSPect scores out of {len(df)} variants.")
    print(f"Arm 2 Disease Tendency is OFFICIALLY COMPLETE!")
    print(f"Saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()