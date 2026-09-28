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
import re

# Paths
INPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm2_DiseaseTendency_Final.csv')
METASNP_FILE = project_path('data/processed/predictors/meta predictors/MetaSNP_output.txt')
OUTPUT_CSV = project_path('outputs/variant_prioritization/EGFR_Arm2_Master_Standardized.csv')

def safe_float(val):
    try:
        if pd.isna(val) or val == '.' or str(val).strip() == '':
            return np.nan
        return float(val)
    except ValueError:
        return np.nan

def parse_metasnp(filepath):
    """Parses the Meta-SNP output text format with score on the subsequent line"""
    results = {}
    with open(filepath, 'r') as f:
        content = f.read()
        # Find mutation blocks
        # Pattern: MutationName, then Neutral/Disease, then score on next line
        pattern = r"([A-Z]\d+[A-Z])\s+(?:Neutral|Disease|NA)\s+(?:Neutral|Disease|NA)\s+(?:Neutral|Disease|NA)\s+(?:Neutral|Disease|NA)\s+(Neutral|Disease|NA).*?\n\s+([\d\.-]+)"
        matches = re.findall(pattern, content, re.DOTALL)
        for mut, pred, score in matches:
            results[mut] = (pred, score)
    return results

def main():
    print(f"Loading Arm 2 Partial data: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    
    # 1. Parse and Merge Meta-SNP
    print("Parsing Meta-SNP output...")
    metasnp_data = parse_metasnp(METASNP_FILE)
    
    df['MetaSNP_pred_raw'] = '.'
    df['MetaSNP_score'] = '.'
    
    aa_map = {'Ala':'A','Cys':'C','Asp':'D','Glu':'E','Phe':'F','Gly':'G','His':'H','Ile':'I','Lys':'K','Leu':'L','Met':'M','Asn':'N','Pro':'P','Gln':'Q','Arg':'R','Ser':'S','Thr':'T','Val':'V','Trp':'W','Tyr':'Y'}

    for i, row in df.iterrows():
        hgvsp = str(row['HGVSp'])
        match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvsp)
        if match:
            wt_3, pos, mut_3 = match.groups()
            if wt_3 in aa_map and mut_3 in aa_map:
                key = f"{aa_map[wt_3]}{pos}{aa_map[mut_3]}"
                if key in metasnp_data:
                    df.at[i, 'MetaSNP_pred_raw'] = str(metasnp_data[key][0])
                    df.at[i, 'MetaSNP_score'] = str(metasnp_data[key][1])

    # 2. Standardization Logic
    print("Standardizing all Arm 2 tools...")
    
    df['MetaSVM_Standard'] = df['MetaSVM_pred'].apply(lambda x: "Deleterious" if str(x) == 'D' else ("Benign" if str(x) == 'T' else "Unknown"))
    df['MetaLR_Standard'] = df['MetaLR_pred'].apply(lambda x: "Deleterious" if str(x) == 'D' else ("Benign" if str(x) == 'T' else "Unknown"))
    df['AlphaMissense_Standard'] = df['AlphaMissense_score'].apply(lambda x: "Deleterious" if safe_float(x) > 0.56 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown"))
    df['SuSPect_Standard'] = df['SuSPect_score'].apply(lambda x: "Deleterious" if safe_float(x) >= 50 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown"))
    df['PHDSNPg_Standard'] = df['PHDSNPg_pred'].apply(lambda x: "Deleterious" if str(x) == 'Pathogenic' else ("Benign" if str(x) == 'Benign' else "Unknown"))
    df['MetaSNP_Standard'] = df['MetaSNP_pred_raw'].apply(lambda x: "Deleterious" if str(x) == 'Disease' else ("Benign" if str(x) == 'Neutral' else "Unknown"))
    df['VARITY_Standard'] = df['VARITY_R_score'].apply(lambda x: "Deleterious" if safe_float(x) > 0.5 else ("Benign" if not np.isnan(safe_float(x)) else "Unknown"))

    # 3. Consensus Calculation
    tools = ['MetaSVM_Standard', 'MetaLR_Standard', 'AlphaMissense_Standard', 'SuSPect_Standard', 'PHDSNPg_Standard', 'MetaSNP_Standard', 'VARITY_Standard']
    df['Arm2_Tools_Available'] = df[tools].apply(lambda x: (x != "Unknown").sum(), axis=1)
    df['Arm2_Deleterious_Count'] = df[tools].apply(lambda x: (x == "Deleterious").sum(), axis=1)
    
    df['Arm2_Consensus'] = df.apply(lambda row: "PASS" if row['Arm2_Deleterious_Count'] >= 4 else "FAIL", axis=1)

    front_cols = ['Genomic_Location', 'HGVSp', 'Existing_variation', 'Arm2_Consensus', 'Arm2_Deleterious_Count', 'Arm2_Tools_Available']
    final_cols = front_cols + tools + [c for c in df.columns if c not in front_cols + tools]
    df[final_cols].to_csv(OUTPUT_CSV, index=False)
    
    print(f"\nSUCCESS! Arm 2 Standardized.")
    print(f"Total mutations passing Arm 2 (>= 4 tools): {len(df[df['Arm2_Consensus'] == 'PASS'])}")

if __name__ == "__main__":
    main()