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

import os
import pandas as pd
import re

# ==========================================
# 1. FILE PATHS & SETUP
# ==========================================
s3_path = project_path('outputs/variant_prioritization/S3.xlsx') 
mutant_dir = project_path('outputs/structural_analysis/mutants')
cxc_output = project_path('outputs/structural_analysis/run_mutations.cxc')

os.makedirs(mutant_dir, exist_ok=True)

# ==========================================
# 2. LOAD MUTATION DATA
# ==========================================
df = pd.read_excel(s3_path, skiprows=3, engine='openpyxl') 
pass_muts = df[df['Consensus_Stability_ddG'] == 'PASS']['Mutation'].dropna().unique().tolist()

print(f"Found {len(pass_muts)} unique PASS mutations. Generating ChimeraX automation script...")

aa_map = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 
          'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 
          'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', 
          'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}

# ==========================================
# 3. GENERATE CHIMERAX SCRIPT
# ==========================================
# Use project-relative paths so the commands work on Windows and Linux.
# Set ChimeraX's working directory to the project root before opening the file.
wt_pdb = 'outputs/structural_analysis/EGFR_mature_WT.pdb'
mutant_relative_dir = 'outputs/structural_analysis/mutants'

with open(cxc_output, 'w') as f:
    for mut in pass_muts:
        match = re.match(r'([A-Za-z])(\d+)([A-Za-z])', mut)
        if match:
            res_num = match.group(2)
            new_aa = match.group(3).upper()
            new_aa_3 = aa_map.get(new_aa)
            
            if new_aa_3:
                f.write(f'open "{wt_pdb}"\n')
                f.write(f"swapaa :{res_num} {new_aa_3} log false\n")
                f.write(f'save "{mutant_relative_dir}/{mut}.pdb" format pdb\n')
                f.write("close session\n\n")

print(f"\nSUCCESS! ChimeraX script saved at: {cxc_output}. Run it from the project root.")
