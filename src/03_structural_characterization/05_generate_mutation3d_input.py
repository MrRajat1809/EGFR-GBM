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
OUTPUT_TXT = project_path('outputs/structural_analysis/3d_spatial_clustering/EGFR_42_Mutation3D_Input.txt')

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
    drivers = df[df['Final_Status'] == 'High-Confidence Driver'].copy()
    
    mut_list = []

    for _, row in drivers.iterrows():
        hgvsp = str(row['HGVSp'])
        match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvsp)
        
        if match:
            wt_3, pos, mut_3 = match.groups()
            if wt_3 in aa_map and mut_3 in aa_map:
                mut_list.append(f"{aa_map[wt_3]}{pos}{aa_map[mut_3]}")

    # Mutation3D Batch Format: [Identifier] [Mut1] [Mut2] [Mut3]...
    batch_line = f"EGFR {' '.join(mut_list)}"

    with open(OUTPUT_TXT, 'w') as f:
        f.write(batch_line + "\n")

    print(f"\nSUCCESS! Formatted {len(mut_list)} variants for Mutation3D Batch Upload.")
    print(f"File saved to: {OUTPUT_TXT}")
    print(f"Preview: {batch_line[:75]}...")

if __name__ == "__main__":
    main()