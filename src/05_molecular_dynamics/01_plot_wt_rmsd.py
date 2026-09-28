# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('outputs/gromacs/WT_Erlotinib')

import matplotlib.pyplot as plt

def read_xvg(filename):
    x, y = [], []
    with open(filename, 'r') as f:
        for line in f:
            # Ignore headers and comments in the xvg file
            if not line.startswith(('@', '#')):
                parts = line.split()
                if len(parts) >= 2:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))
    return x, y

# Read the data
time_prot, rmsd_prot = read_xvg('rmsd_protein_200.xvg')
time_lig, rmsd_lig = read_xvg('rmsd_ligand_200.xvg')

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(time_prot, rmsd_prot, label='Protein Backbone', color='blue', linewidth=1.5, alpha=0.8)
plt.plot(time_lig, rmsd_lig, label='Erlotinib (Ligand)', color='red', linewidth=1.5, alpha=0.8)

# Format the graph
plt.title('RMSD of WT EGFR + Erlotinib (100 ns)', fontsize=14, fontweight='bold')
plt.xlabel('Time (ns)', fontsize=12)
plt.ylabel('RMSD (nm)', fontsize=12)
plt.ylim(0, 0.5) # Sets the Y-axis from 0 to 0.5 nm for a clean look
plt.legend(loc='upper right', fontsize=11)
plt.grid(True, linestyle='--', alpha=0.5)

# Save and show
plt.tight_layout()
plt.savefig('WT_Erlotinib_RMSD.png', dpi=300)
print("Graph saved successfully as WT_Erlotinib_RMSD.png!")