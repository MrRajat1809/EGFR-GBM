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
BASE_CSV = project_path('data/processed/variants/EGFR_TCGA_51_Base.csv')
OUTPUT_DIR = project_path('data/processed/predictors/phdsnp')
OUTPUT_FILE = f"{OUTPUT_DIR}/EGFR_TCGA_51_PHDSNP_Input.vcf"

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Loading Base File: {BASE_CSV}")
    df = pd.read_csv(BASE_CSV)
    
    # Lists to build our VCF-style format
    chroms, poses, ids, refs, alts = [], [], [], [], []
    
    for _, row in df.iterrows():
        # 1. Parse 'Genomic_Location' (Format: 7_55233109_G/A)
        loc_parts = str(row['Genomic_Location']).split('_')
        chrom = loc_parts[0]
        pos = loc_parts[1]
        ref, alt = loc_parts[2].split('/')
        
        # 2. Parse 'Existing_variation' for the ID (rsID or COSMIC)
        var_id = str(row['Existing_variation'])
        if pd.isna(var_id) or var_id == 'nan' or var_id.strip() == '':
            var_id = '.'
        else:
            # If multiple IDs exist (e.g., rs123,COSM456), grab the first one
            var_id = var_id.split(',')[0]
            
        chroms.append(chrom)
        poses.append(pos)
        ids.append(var_id)
        refs.append(ref)
        alts.append(alt)

    # Build the final DataFrame
    vcf_df = pd.DataFrame({
        '#CHROM': chroms,
        'POS': poses,
        'ID': ids,
        'REF': refs,
        'ALT': alts
    })

    # Save as a tab-separated VCF-style file
    vcf_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    
    print(f"\nSUCCESS! Created PHD-SNPg input file for exactly {len(vcf_df)} variants.")
    print(f"Saved to: {OUTPUT_FILE}")
    print("You can upload this single file directly to the PHD-SNPg web server!")

if __name__ == "__main__":
    main()