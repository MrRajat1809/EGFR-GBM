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
PHDSNP_PATH = project_path('data/processed/predictors/phdsnp/PHD-SNPg_output.txt')
OUTPUT_DIR = project_path('outputs/variant_prioritization')
OUTPUT_CSV = f"{OUTPUT_DIR}/EGFR_Arm2_DiseaseTendency_Partial.csv"

TARGET_COLS = [
    'MetaSVM_score', 'MetaSVM_pred',
    'MetaLR_score', 'MetaLR_pred',
    'AlphaMissense_score', 'AlphaMissense_pred',
    'VARITY_R_score', 'VARITY_ER_score'
]

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(BASE_CSV)
    
    # 1. Initialize empty dbNSFP columns
    for col in TARGET_COLS:
        df[col] = '.'
        
    # 2. Extract from dbNSFP
    print(f"Teleporting into dbNSFP for Arm 2 tools...")
    tb = pysam.TabixFile(DBNSFP_PATH)
    header = list(tb.header)[0].lstrip('#').split('\t')
    col_indices = {col: header.index(col) for col in TARGET_COLS if col in header}
    
    idx_chr, idx_pos, idx_ref, idx_alt = header.index('chr'), header.index('pos(1-based)'), header.index('ref'), header.index('alt')

    for i, row in df.iterrows():
        loc_parts = str(row['Genomic_Location']).split('_')
        chrom, pos = loc_parts[0], int(loc_parts[1])
        ref, alt = loc_parts[2].split('/')
        
        try:
            records = tb.fetch(chrom, pos - 1, pos)
            for record in records:
                fields = record.split('\t')
                if fields[idx_ref] == ref and fields[idx_alt] == alt:
                    for col, col_idx in col_indices.items():
                        valid_vals = [v for v in fields[col_idx].split(';') if v != '.' and v != '']
                        df.at[i, col] = valid_vals[0] if valid_vals else '.'
                    break
        except ValueError:
            pass

    # 3. Merge PHD-SNPg Data
    print(f"Merging PHD-SNPg results from: {PHDSNP_PATH}")
    if os.path.exists(PHDSNP_PATH):
        # Read PHD-SNPg (skipping the meta-info header starting with ##)
        df_phd = pd.read_csv(PHDSNP_PATH, sep='\t', comment='#', header=None)
        
        # Manually extract the column header from the file
        with open(PHDSNP_PATH, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    phd_cols = line.strip().lstrip('#').split('\t')
                    break
        df_phd.columns = phd_cols
        
        # Rename columns to avoid conflicts and make merging clean
        df_phd.rename(columns={'SCORE': 'PHDSNPg_score', 'PREDICTION': 'PHDSNPg_pred'}, inplace=True)
        
        # Create a matching string (e.g., '7_55233109_G/A')
        df_phd['Match_Location'] = df_phd['CHROM'].astype(str) + '_' + df_phd['POS'].astype(str) + '_' + df_phd['REF'] + '/' + df_phd['ALT']
        
        # Merge exactly on our unique identifier
        df = pd.merge(df, df_phd[['Match_Location', 'PHDSNPg_score', 'PHDSNPg_pred']], 
                      left_on='Genomic_Location', right_on='Match_Location', how='left')
        df.drop(columns=['Match_Location'], inplace=True)
        print("PHD-SNPg successfully merged!")
    else:
        print(f"Warning: Could not find PHD-SNPg file at {PHDSNP_PATH}")

    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSUCCESS! Arm 2 (Partial) saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()