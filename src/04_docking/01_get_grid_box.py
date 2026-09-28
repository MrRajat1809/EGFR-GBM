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

from Bio.PDB import PDBParser
import numpy as np

pdb_path = project_path('outputs/structural_analysis/docking/targets/WT_TKD_clean.pdb')

# These are the classic ATP-pocket residues where Erlotinib binds
atp_pocket_residues = [718, 726, 743, 745, 790, 793, 844]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("TKD", pdb_path)

coords = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[1] in atp_pocket_residues and 'CA' in residue:
                coords.append(residue['CA'].get_coord())

coords = np.array(coords)
center = coords.mean(axis=0)

print(f"--- VINA GRID BOX COORDINATES ---")
print(f"center_x = {center[0]:.3f}")
print(f"center_y = {center[1]:.3f}")
print(f"center_z = {center[2]:.3f}")
print(f"size_x = 20.0")
print(f"size_y = 20.0")
print(f"size_z = 20.0")
