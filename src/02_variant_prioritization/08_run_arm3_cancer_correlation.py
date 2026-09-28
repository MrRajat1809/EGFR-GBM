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

# Paths
BASE_CSV = project_path('data/processed/variants/EGFR_TCGA_51_Base.csv')
DBNSFP_PATH = project_path('data/reference/dbnsfp_5.3.1a/dbNSFP5.3.1a_grch37.gz')
OUTPUT_DIR = project_path('outputs/variant_prioritization')
OUTPUT_CSV = f"{OUTPUT_DIR}/EGFR_Arm3_CancerCorrelation_Partial.csv"

# The specific Arm 3 tools found in your dbNSFP header
TARGET_COLS = [
    'fathmm-XF_coding_score', 'fathmm-XF_coding_pred',
    'MetaRNN_score', 'MetaRNN_pred',
    'VEST4_score',
    'BayesDel_addAF_score', 'BayesDel_addAF_pred'
]

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(BASE_CSV)
    
    # Initialize empty columns
    for col in TARGET_COLS:
        df[col] = '.'
        
    print(f"Teleporting into dbNSFP for Arm 3 tools...")
    try:
        tb = pysam.TabixFile(DBNSFP_PATH)
    except Exception as e:
        print(f"Error opening dbNSFP file: {e}")
        return

    header = list(tb.header)[0].lstrip('#').split('\t')
    col_indices = {col: header.index(col) for col in TARGET_COLS if col in header}
    
    idx_chr = header.index('chr')
    idx_pos = header.index('pos(1-based)')
    idx_ref = header.index('ref')
    idx_alt = header.index('alt')

    match_count = 0

    for i, row in df.iterrows():
        loc_parts = str(row['Genomic_Location']).split('_')
        chrom = loc_parts[0]
        pos = int(loc_parts[1])
        ref, alt = loc_parts[2].split('/')
        
        try:
            records = tb.fetch(chrom, pos - 1, pos)
            for record in records:
                fields = record.split('\t')
                
                # Strict matching
                if fields[idx_ref] == ref and fields[idx_alt] == alt:
                    match_count += 1
                    for col, col_idx in col_indices.items():
                        valid_vals = [v for v in fields[col_idx].split(';') if v != '.' and v != '']
                        df.at[i, col] = valid_vals[0] if valid_vals else '.'
                    break
        except ValueError:
            pass 

    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSUCCESS! Retrieved Arm 3 scores for {match_count} variants.")
    print(f"Arm 3 Partial Data safely stored at: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()