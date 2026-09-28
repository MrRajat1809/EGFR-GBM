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
import pysam
import os

# Your exact WSL paths
BASE_CSV = project_path('data/processed/variants/EGFR_TCGA_51_Base.csv')
DBNSFP_PATH = project_path('data/reference/dbnsfp_5.3.1a/dbNSFP5.3.1a_grch37.gz')
OUTPUT_DIR = project_path('outputs/variant_prioritization')
OUTPUT_CSV = f"{OUTPUT_DIR}/EGFR_Arm1_Pathogenicity.csv"

# The exact column names we found in your dbNSFP header
TARGET_COLS = [
    'SIFT_score', 'SIFT_pred',
    'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred',
    'PROVEAN_score', 'PROVEAN_pred',
    'CADD_phred', 'CADD_raw',
    'REVEL_score',
    'MutationTaster_score', 'MutationTaster_pred',
    'MutationAssessor_score', 'MutationAssessor_pred'
]

def main():
    # 1. Create the results directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Loading Base File: {BASE_CSV}")
    df = pd.read_csv(BASE_CSV)
    
    # Initialize empty columns for our tools
    for col in TARGET_COLS:
        df[col] = '.'
        
    print(f"Opening dbNSFP Database: {DBNSFP_PATH}")
    try:
        tb = pysam.TabixFile(DBNSFP_PATH)
    except Exception as e:
        print(f"Error opening dbNSFP file (Check your path or pysam installation): {e}")
        return

    # Extract dbNSFP header to map the exact column indices dynamically
    header = list(tb.header)[0].lstrip('#').split('\t')
    col_indices = {col: header.index(col) for col in TARGET_COLS if col in header}
    
    # Get indices for the coordinate columns to ensure a perfect biological match
    idx_chr = header.index('chr')
    idx_pos = header.index('pos(1-based)')
    idx_ref = header.index('ref')
    idx_alt = header.index('alt')

    print("Teleporting into dbNSFP to extract Arm 1 scores for 51 variants...")
    match_count = 0

    # 2. Iterate over the 51 variants
    for i, row in df.iterrows():
        # Genomic_Location is formatted as '7_55233109_G/A'
        loc_parts = str(row['Genomic_Location']).split('_')
        chrom = loc_parts[0]
        pos = int(loc_parts[1])
        ref, alt = loc_parts[2].split('/')
        
        # 3. Query the Tabix file (uses 0-based half-open intervals: pos-1 to pos)
        try:
            records = tb.fetch(chrom, pos - 1, pos)
            for record in records:
                fields = record.split('\t')
                
                # 4. Verify the Reference and Alternate alleles match exactly!
                if fields[idx_ref] == ref and fields[idx_alt] == alt:
                    match_count += 1
                    
                    # Extract our targeted scores
                    for col, col_idx in col_indices.items():
                        raw_val = fields[col_idx]
                        
                        # dbNSFP sometimes lists multiple scores separated by ';' for different transcripts.
                        # We grab the first valid score we find.
                        valid_vals = [v for v in raw_val.split(';') if v != '.' and v != '']
                        df.at[i, col] = valid_vals[0] if valid_vals else '.'
                        
                    break # We found our exact match, move to the next variant
                    
        except ValueError:
            pass # Standard Tabix exception if a coordinate is entirely missing from the database

    # 5. Save the enriched dataset
    df.to_csv(OUTPUT_CSV, index=False)
    
    print(f"\nSUCCESS! Matched {match_count} out of {len(df)} variants in dbNSFP.")
    print(f"Arm 1 Pathogenicity data safely stored at: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()