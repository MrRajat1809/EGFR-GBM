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

VEP_TSV = project_path('data/processed/variants/EGFR_TCGA_51_Somatic_VEP.tsv')
OUTPUT_CSV = project_path('data/processed/variants/EGFR_TCGA_51_Base.csv')

def main():
    print(f"Reading VEP output: {VEP_TSV}...")
    
    # 1. Read the TSV, skipping the top meta-information lines (##)
    df_vep = pd.read_csv(VEP_TSV, sep='\t', comment='#', header=None)
    
    # 2. Extract the actual column headers (the line starting with #Uploaded_variation)
    header_line = ""
    with open(VEP_TSV, "r") as f:
        for line in f:
            if line.startswith('#Uploaded_variation'):
                header_line = line.strip()
                break
                
    # Clean up the header line and apply to dataframe
    columns = header_line.lstrip('#').split('\t')
    df_vep.columns = columns
    print(f"Total annotations across all transcripts: {len(df_vep)}")

    # 3. Filter strictly for the Canonical Transcript (ENST00000275493)
    df_canonical = df_vep[df_vep['Feature'].str.contains('ENST00000275493', na=False)].copy()
    print(f"Total somatic variants mapped to Canonical Transcript: {len(df_canonical)}")

    # 4. Optional: Rename 'Uploaded_variation' to 'Location' or 'Variant_ID' for clarity
    df_canonical.rename(columns={'Uploaded_variation': 'Genomic_Location'}, inplace=True)

    # 5. Export to CSV
    df_canonical.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSUCCESS! Your perfect 51-row baseline file is saved at: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()