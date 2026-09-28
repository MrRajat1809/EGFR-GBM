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

INPUT_CSV = project_path('outputs/variant_prioritization/EGFR_TCGA_Glioblastoma_Master_Consensus.csv')
OUTPUT_TXT = project_path('outputs/structural_analysis/thermodynamic_stability_ranking/EGFR_PremPS_Input.txt')

aa_map = {'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Cys':'C', 'Gln':'Q', 
          'Glu':'E', 'Gly':'G', 'His':'H', 'Ile':'I', 'Leu':'L', 'Lys':'K', 
          'Met':'M', 'Phe':'F', 'Pro':'P', 'Ser':'S', 'Thr':'T', 'Trp':'W', 
          'Tyr':'Y', 'Val':'V'}

def main():
    if not os.path.exists(INPUT_CSV):
        print(f"Error: Could not find {INPUT_CSV}")
        return

    df = pd.read_csv(INPUT_CSV)
    drivers = df[df['Final_Status'] == 'High-Confidence Driver'].copy()
    
    # Add the exact header PremPS is asking for
    premps_list = ["#Chain to Mutate\tResidue\tMutant Residue"]
    
    for _, row in drivers.iterrows():
        match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', str(row['HGVSp']))
        if match:
            wt_3, pos, mut_3 = match.groups()
            if wt_3 in aa_map and mut_3 in aa_map:
                wt_1 = aa_map[wt_3]
                mut_1 = aa_map[mut_3]
                
                # Tab-separated format: Chain(A) \t WT+Pos(C620) \t Mut(Y)
                premps_list.append(f"A\t{wt_1}{pos}\t{mut_1}")

    with open(OUTPUT_TXT, 'w') as f:
        f.write("\n".join(premps_list))

    print(f"SUCCESS! PremPS input list generated with {len(premps_list)-1} variants.")
    print(f"Saved to: {OUTPUT_TXT}")
    print("\nFile Preview:")
    print("\n".join(premps_list[:4]))

if __name__ == "__main__":
    main()