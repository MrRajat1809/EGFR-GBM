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

# 1. Paths to your three master standardized tables
BASE_DIR = project_path('outputs/variant_prioritization')
PATH_FILE = f"{BASE_DIR}/EGFR_Arm1_Master_Standardized.csv"
DIS_FILE  = f"{BASE_DIR}/EGFR_Arm2_Master_Standardized.csv"
CAN_FILE  = f"{BASE_DIR}/EGFR_Arm3_Master_Standardized.csv"

# Output file path
OUTPUT_FILE = f"{BASE_DIR}/EGFR_TCGA_Glioblastoma_Master_Consensus.csv"

def main():
    print("1. Loading standardized Arm master tables...")
    # Verify file existence
    for f in [PATH_FILE, DIS_FILE, CAN_FILE]:
        if not os.path.exists(f):
            print(f"Error: Required file not found at {f}")
            return

    # Load tables
    df_path = pd.read_csv(PATH_FILE)
    df_dis  = pd.read_csv(DIS_FILE)
    df_can  = pd.read_csv(CAN_FILE)

    # Standard merge key for this pipeline
    key = 'Genomic_Location'

    print("2. Performing multi-arm intersection (Successive Inner Joins)...")
    
    # We select the consensus and count columns from each arm for a clean master view
    m1 = df_path[[key, 'HGVSp', 'Existing_variation', 'Arm1_Consensus', 'Arm1_Deleterious_Count']]
    m2 = df_dis[[key, 'Arm2_Consensus', 'Arm2_Deleterious_Count']]
    m3 = df_can[[key, 'Arm3_Consensus', 'Arm3_Oncogenic_Count']]

    # Merge into a single master table
    merged = pd.merge(m1, m2, on=key, how='inner')
    final_master = pd.merge(merged, m3, on=key, how='inner')

    print("3. Tagging Triple-Hit variants (High-Confidence Glioblastoma Drivers)...")
    
    # A variant is a "Triple Hit" if it passes the threshold of ALL 3 ARMS
    final_master['Triple_Arm_Consensus'] = (
        (final_master['Arm1_Consensus'] == 'PASS') &
        (final_master['Arm2_Consensus'] == 'PASS') &
        (final_master['Arm3_Consensus'] == 'PASS')
    )

    # Label the classification for easier filtering in ChimeraX later
    final_master['Final_Status'] = final_master['Triple_Arm_Consensus'].map({
        True: 'High-Confidence Driver',
        False: 'Likely Passenger/Neutral'
    })

    # Save the consolidated table
    final_master.to_csv(OUTPUT_FILE, index=False)
    
    # Statistics for the Pipeline Report
    total_variants = len(final_master)
    triple_hits = final_master['Triple_Arm_Consensus'].sum()
    
    print("\n" + "="*45)
    print("      PIPELINE SUMMARY: TRIPLE ARM OVERLAP")
    print("="*45)
    print(f"Total Unique TCGA Mutations:   {total_variants}")
    print(f"High-Confidence Drivers:       {triple_hits}")
    print(f"Likely Passenger Mutations:    {total_variants - triple_hits}")
    print("="*45)
    print(f"Final Deliverable saved to: {OUTPUT_FILE}")

    if triple_hits > 0:
        print("\nTOP 5 HIGH-CONFIDENCE DRIVERS IDENTIFIED:")
        drivers = final_master[final_master['Triple_Arm_Consensus'] == True]
        print(drivers[['HGVSp', 'Existing_variation', 'Arm3_Oncogenic_Count']].head(5).to_string(index=False))

if __name__ == "__main__":
    main()