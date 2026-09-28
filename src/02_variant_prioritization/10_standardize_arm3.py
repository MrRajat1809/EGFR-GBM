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
import numpy as np
import os

# Paths
INPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm3_CancerCorrelation_Final.csv')
OUTPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm3_Master_Standardized.csv')

def safe_float(val):
    """Safely convert dbNSFP '.' or missing values to float"""
    try:
        if pd.isna(val) or val == '.' or str(val).strip() == '':
            return np.nan
        return float(val)
    except ValueError:
        return np.nan

def main():
    print(f"Loading raw Arm 3 data: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    
    # 1. FATHMM-XF (> 0.5 is Oncogenic)
    df['FATHMM_Standard'] = df['fathmm-XF_coding_score'].apply(
        lambda x: "Oncogenic" if safe_float(x) > 0.5 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # 2. MetaRNN (> 0.5 is Oncogenic)
    df['MetaRNN_Standard'] = df['MetaRNN_score'].apply(
        lambda x: "Oncogenic" if safe_float(x) > 0.5 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # 3. VEST4 (> 0.5 is Oncogenic)
    df['VEST4_Standard'] = df['VEST4_score'].apply(
        lambda x: "Oncogenic" if safe_float(x) > 0.5 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # 4. BayesDel (> 0.16 is Oncogenic)
    df['BayesDel_Standard'] = df['BayesDel_addAF_score'].apply(
        lambda x: "Oncogenic" if safe_float(x) > 0.16 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # 5. CScape (>= 0.5 is Oncogenic)
    df['CScape_Standard'] = df['CScape_score'].apply(
        lambda x: "Oncogenic" if safe_float(x) >= 0.5 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    print("Calculating Arm 3 Consensus...")
    tools = ['FATHMM_Standard', 'MetaRNN_Standard', 'VEST4_Standard', 'BayesDel_Standard', 'CScape_Standard']
    
    # Calculate how many tools successfully scored the variant
    df['Arm3_Tools_Available'] = df[tools].apply(lambda x: (x != "Unknown").sum(), axis=1)
    
    # Calculate how many tools called it Oncogenic
    df['Arm3_Oncogenic_Count'] = df[tools].apply(lambda x: (x == "Oncogenic").sum(), axis=1)
    
    # ARM 3 CONSENSUS LOGIC: Must pass at least 3 out of 5 cancer tools
    df['Arm3_Consensus'] = df.apply(
        lambda row: "PASS" if row['Arm3_Oncogenic_Count'] >= 3 else "FAIL", 
        axis=1
    )
    
    # Reorder columns to make it look highly professional
    front_cols = ['Genomic_Location', 'HGVSp', 'Existing_variation', 'Arm3_Consensus', 'Arm3_Oncogenic_Count', 'Arm3_Tools_Available']
    standard_cols = tools
    
    # Put essential columns first, then standard predictions, then raw scores at the back
    other_cols = [c for c in df.columns if c not in front_cols + standard_cols]
    final_cols = front_cols + standard_cols + other_cols
    df = df[final_cols]
    
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSUCCESS! Standardized {len(df)} variants for Arm 3.")
    print(f"Total mutations passing Arm 3 (>= 3 tools): {len(df[df['Arm3_Consensus'] == 'PASS'])}")
    print(f"Saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()