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

"""
EGFR Variant Pathogenicity Standardization (Arm 1)
--------------------------------------------------
This script standardizes pathogenicity scores from multiple bioinformatic 
tools (dbNSFP) into a binary 'Deleterious' vs 'Benign' classification 
for EGFR variants.

Methodology (Arm 1 Consensus):
A variant is classified as 'PASS' if at least 4 out of 7 tools predict 
it as deleterious based on the following literature-standard thresholds:
- SIFT: < 0.05
- PolyPhen2 HVAR: > 0.908 (Probably Damaging)
- PROVEAN: < -2.5
- CADD Phred: > 20.0 (Top 1% of human genome)
- REVEL: > 0.5
- MutationTaster: 'D' (Disease causing) or 'A' (Abnormal)
- MutationAssessor: > 1.93

"""

import pandas as pd
import numpy as np

# --- CONFIGURATION ---
INPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm1_Pathogenicity.csv')
OUTPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm1_Master_Standardized.csv')

def safe_float(val):
    """
    Handles missing data and dbNSFP null characters ('.').
    Ensures the pipeline doesn't break on variants missing specific scores.
    """
    try:
        if pd.isna(val) or val == '.' or str(val).strip() == '':
            return np.nan
        return float(val)
    except ValueError:
        return np.nan

def main():
    """
    Main execution block: Loads data, applies thresholds, and calculates consensus.
    """
    print(f"Loading raw Arm 1 data: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    
    # --- 1. TOOL-SPECIFIC STANDARDIZATION ---
    
    # SIFT: Scores range from 0 to 1. Lower scores are more deleterious.
    df['SIFT_Standard'] = df['SIFT_score'].apply(
        lambda x: "Deleterious" if safe_float(x) < 0.05 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # PolyPhen2 HVAR: Scores closer to 1 are more damaging. 
    # 0.908 is the high-confidence threshold for 'Probably Damaging'.
    df['PolyPhen_Standard'] = df['Polyphen2_HVAR_score'].apply(
        lambda x: "Deleterious" if safe_float(x) > 0.908 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # PROVEAN: Scores < -2.5 are considered deleterious for protein function.
    df['PROVEAN_Standard'] = df['PROVEAN_score'].apply(
        lambda x: "Deleterious" if safe_float(x) < -2.5 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # CADD Phred: A rank-based score. 20 means the variant is in the top 1% 
    # of most deleterious substitutions in the human genome.
    df['CADD_Standard'] = df['CADD_phred'].apply(
        lambda x: "Deleterious" if safe_float(x) > 20.0 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # REVEL: An ensemble method for missense variants. 0.5 is the suggested threshold.
    df['REVEL_Standard'] = df['REVEL_score'].apply(
        lambda x: "Deleterious" if safe_float(x) > 0.5 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # MutationTaster: Uses categorical predictions.
    # D = Disease causing; A = Disease causing (automatic).
    def std_mt(row):
        pred = str(row['MutationTaster_pred']).upper()
        if 'D' in pred or 'A' in pred: return "Deleterious"
        if 'N' in pred or 'P' in pred: return "Benign"
        return "Unknown"
    df['MutationTaster_Standard'] = df.apply(std_mt, axis=1)
    
    # MutationAssessor: High/Medium impact (Functional) scores > 1.93.
    df['MutationAssessor_Standard'] = df['MutationAssessor_score'].apply(
        lambda x: "Deleterious" if safe_float(x) > 1.93 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown")
    )
    
    # --- 2. CONSENSUS CALCULATION ---
    
    print("Calculating Arm 1 Consensus...")
    tools = ['SIFT_Standard', 'PolyPhen_Standard', 'PROVEAN_Standard', 'CADD_Standard', 
             'REVEL_Standard', 'MutationTaster_Standard', 'MutationAssessor_Standard']
    
    # Track data density (how many tools actually provided a prediction)
    df['Arm1_Tools_Available'] = df[tools].apply(lambda x: (x != "Unknown").sum(), axis=1)
    
    # Count the number of tools that reached a 'Deleterious' conclusion
    df['Arm1_Deleterious_Count'] = df[tools].apply(lambda x: (x == "Deleterious").sum(), axis=1)
    
    # ARM 1 FINAL CONSENSUS: Majority rule (4/7)
    df['Arm1_Consensus'] = df.apply(
        lambda row: "PASS" if row['Arm1_Deleterious_Count'] >= 4 else "FAIL", 
        axis=1
    )
    
    # --- 3. EXPORT & FORMATTING ---
    
    # Prioritize consensus results in the column order for easy review
    front_cols = ['Genomic_Location', 'HGVSp', 'Existing_variation', 'Arm1_Consensus', 
                  'Arm1_Deleterious_Count', 'Arm1_Tools_Available']
    
    other_cols = [c for c in df.columns if c not in front_cols + tools]
    final_cols = front_cols + tools + other_cols
    df = df[final_cols]
    
    # Save results
    df.to_csv(OUTPUT_CSV, index=False)
    
    # Final Reporting
    print(f"\n--- EXECUTION SUMMARY ---")
    print(f"Total variants processed: {len(df)}")
    print(f"Variants passing Arm 1:  {len(df[df['Arm1_Consensus'] == 'PASS'])}")
    print(f"Output saved to:         {OUTPUT_CSV}")

if __name__ == "__main__":
    main()