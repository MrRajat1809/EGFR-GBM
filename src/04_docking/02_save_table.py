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

data = [
    {"Target": "WT_TKD", "Ligand": "Erlotinib", "Mode": 1, "Affinity (kcal/mol)": -6.595, "RMSD_lb": 0.0, "RMSD_ub": 0.0},
    {"Target": "WT_TKD", "Ligand": "Erlotinib", "Mode": 2, "Affinity (kcal/mol)": -6.462, "RMSD_lb": 1.643, "RMSD_ub": 2.114},
    {"Target": "WT_TKD", "Ligand": "Erlotinib", "Mode": 3, "Affinity (kcal/mol)": -6.390, "RMSD_lb": 2.838, "RMSD_ub": 5.180},
    {"Target": "L861Q_TKD", "Ligand": "Erlotinib", "Mode": 1, "Affinity (kcal/mol)": -6.591, "RMSD_lb": 0.0, "RMSD_ub": 0.0},
    {"Target": "L861Q_TKD", "Ligand": "Erlotinib", "Mode": 2, "Affinity (kcal/mol)": -6.568, "RMSD_lb": 1.490, "RMSD_ub": 1.951},
    {"Target": "L861Q_TKD", "Ligand": "Erlotinib", "Mode": 3, "Affinity (kcal/mol)": -6.433, "RMSD_lb": 4.420, "RMSD_ub": 8.100},
]

df = pd.DataFrame(data)
df.to_csv(project_path('outputs/variant_prioritization/S5_Erlotinib_Docking.csv'), index=False)
print("Table saved as S5_Erlotinib_Docking.csv")
