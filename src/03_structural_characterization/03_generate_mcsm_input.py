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
import re
import os

# Paths
INPUT_CSV = project_path('outputs/variant_prioritization/EGFR_TCGA_Glioblastoma_Master_Consensus.csv')
OUTPUT_TXT = project_path('outputs/structural_analysis/thermodynamic_stability_ranking/EGFR_42_mCSM_Input.txt')

# 3-letter to 1-letter amino acid dictionary
aa_map = {
    'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Cys':'C',
    'Gln':'Q', 'Glu':'E', 'Gly':'G', 'His':'H', 'Ile':'I',
    'Leu':'L', 'Lys':'K', 'Met':'M', 'Phe':'F', 'Pro':'P',
    'Ser':'S', 'Thr':'T', 'Trp':'W', 'Tyr':'Y', 'Val':'V'
}

def main():
    if not os.path.exists(INPUT_CSV):
        print(f"Error: Could not find {INPUT_CSV}")
        return

    df = pd.read_csv(INPUT_CSV)
    
    # Filter for the 42 winners
    drivers = df[df['Final_Status'] == 'High-Confidence Driver'].copy()
    
    chain_id = 'A'
    mcsm_lines = []

    for _, row in drivers.iterrows():
        hgvsp = str(row['HGVSp'])
        
        # Match p.Ala289Val
        match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvsp)
        
        if match:
            wt_3, pos, mut_3 = match.groups()
            
            if wt_3 in aa_map and mut_3 in aa_map:
                wt_1 = aa_map[wt_3]
                mut_1 = aa_map[mut_3]
                
                # EXACT mCSM FORMAT:
                # Column1 = Chain
                # Column2 = WT+Position+Mut
                mcsm_lines.append(f"{chain_id} {wt_1}{pos}{mut_1}")

    # Remove duplicates (optional but recommended)
    mcsm_lines = sorted(set(mcsm_lines))

    with open(OUTPUT_TXT, 'w') as f:
        for line in mcsm_lines:
            f.write(line + "\n")

    print(f"\nSUCCESS! Formatted {len(mcsm_lines)} variants for mCSM.")
    print(f"File saved to: {OUTPUT_TXT}")

if __name__ == "__main__":
    main()