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

# Your exact WSL paths
MAF_FILE = project_path('data/raw/cbioportal/gbm_tcga_pan_can_atlas_2018/data_mutations.txt')
OUTPUT_DIR = project_path('data/processed/variants')
OUTPUT_CSV = f"{OUTPUT_DIR}/EGFR_TCGA_GBM_Missense_Starting_List.csv"

def main():
    # Ensure the reference directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Loading TCGA Master MAF file: {MAF_FILE}...")
    
    # Check if file exists before trying to read
    if not os.path.exists(MAF_FILE):
        print(f"Error: Could not find the MAF file at {MAF_FILE}. Please check the path!")
        return

    # Read the MAF file (cBioPortal files are tab-separated)
    df = pd.read_csv(MAF_FILE, sep='\t', low_memory=False)
    
    print(f"Total mutations across all genes in GBM cohort: {len(df)}")
    
    # Filter 1: Isolate the EGFR gene
    df_egfr = df[df['Hugo_Symbol'] == 'EGFR']
    print(f"Total EGFR mutations (all types: silent, nonsense, etc.): {len(df_egfr)}")
    
    # Filter 2: Isolate ONLY Missense Mutations
    df_egfr_missense = df_egfr[df_egfr['Variant_Classification'] == 'Missense_Mutation'].copy()
    
    # Select the most critical clinical columns
    cols_to_keep = [
        'Tumor_Sample_Barcode', # The unique patient/tumor ID
        'Hugo_Symbol', 
        'Chromosome', 
        'Start_Position', 
        'Reference_Allele', 
        'Tumor_Seq_Allele2',    # The specific mutated Alternate Allele
        'Variant_Classification',
        'HGVSc',                # The exact DNA change (e.g., c.865G>A)
        'HGVSp_Short'           # The classic protein change (e.g., p.A289T)
    ]
    
    # If the file doesn't have exactly these columns, we'll gracefully grab what is available
    available_cols = [c for c in cols_to_keep if c in df_egfr_missense.columns]
    df_final = df_egfr_missense[available_cols]
    
    # Save the new curated list
    df_final.to_csv(OUTPUT_CSV, index=False)
    
    print(f"\nSUCCESS! Extracted {len(df_final)} somatic EGFR missense mutations.")
    print(f"Saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()