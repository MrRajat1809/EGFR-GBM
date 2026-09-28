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
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
import re

# ==========================================
# 1. FILE PATHS & SETUP
# ==========================================
s3_path = project_path('outputs/variant_prioritization/S3.xlsx')
# USE THE NEW MATURE WT FILE
wt_pdb_path = project_path('outputs/structural_analysis/EGFR_mature_WT.pdb')
output_csv = project_path('outputs/variant_prioritization/S4_Clustering_Table_Final.csv')

print("1. Loading Mutation Data...")
df = pd.read_excel(s3_path, skiprows=3, engine='openpyxl')
pass_df = df[df['Consensus_Stability_ddG'] == 'PASS'].copy()

# Extract numbering EXACTLY as it appears in the clinical data (No +24 offset!)
pass_df['PDB_Residue'] = pass_df['Mutation'].apply(lambda x: int(re.search(r'\d+', x).group()))
pass_df['PremPS_ddG'] = pd.to_numeric(pass_df['PremPS_ddG'], errors='coerce')

# ==========================================
# 2. EXTRACT 3D COORDINATES FROM MATURE PDB
# ==========================================
print("2. Mapping 3D Coordinates...")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("EGFR", wt_pdb_path)
chain = structure[0]['A']

coord_dict = {}
for res in chain:
    if res.id[0] == ' ' and 'CA' in res:
        coord_dict[res.id[1]] = res['CA'].get_coord()

coords, valid_indices = [], []
for idx, row in pass_df.iterrows():
    if row['PDB_Residue'] in coord_dict:
        coords.append(coord_dict[row['PDB_Residue']])
        valid_indices.append(idx)
    else:
        print(f"Warning: Residue {row['PDB_Residue']} not found in PDB.")

pass_df = pass_df.loc[valid_indices].copy()

# ==========================================
# 3. AVERAGE LINKAGE SPATIAL CLUSTERING
# ==========================================
print("3. Performing Average-Linkage Clustering (25 Angstroms)...")
dist_matrix = pdist(np.array(coords))
Z = linkage(dist_matrix, method='average')
pass_df['Raw_Cluster_ID'] = fcluster(Z, t=25, criterion='distance')

# ==========================================
# 4. DOMAIN ANNOTATION (Based on 1-1210 Precursor Numbering)
# ==========================================
def get_domain(resnum):
    if 25 <= resnum <= 189: return "Domain I (ECD)"
    elif 190 <= resnum <= 336: return "Domain II (ECD Tether)"
    elif 337 <= resnum <= 505: return "Domain III (ECD)"
    elif 506 <= resnum <= 645: return "Domain IV (ECD Tether)"
    elif 646 <= resnum <= 668: return "Transmembrane"
    elif 669 <= resnum <= 711: return "Juxtamembrane"
    elif 712 <= resnum <= 979: return "Tyrosine Kinase Domain"
    else: return "C-Terminal Tail"

# ==========================================
# 5. FORMAT FINAL TABLE
# ==========================================
print("4. Formatting Final Table S4...")
cluster_data = []
cluster_counter = 1

for cid in sorted(pass_df['Raw_Cluster_ID'].unique()):
    c_df = pass_df[pass_df['Raw_Cluster_ID'] == cid]
    
    if len(c_df) >= 2:
        muts_with_ddg = [f"{row['Mutation']} ({row['PremPS_ddG']:.2f})" for _, row in c_df.iterrows()]
        rep_mut = c_df.loc[c_df['PremPS_ddG'].idxmax()]['Mutation']
        
        regions = [get_domain(r) for r in c_df['PDB_Residue']]
        primary_region = max(set(regions), key=regions.count)
        
        if "Domain II" in regions and "Domain IV" in regions:
            primary_region = "Domain II/IV Tether Interface"
            
        cluster_data.append({
            'Cluster_ID': f"Cluster {cluster_counter}",
            'Primary_Region': primary_region,
            'Total_Variants': len(c_df),
            'Variants_and_ddG': ", ".join(muts_with_ddg),
            'Representative_Target (Max ddG)': rep_mut
        })
        cluster_counter += 1

final_df = pd.DataFrame(cluster_data)
final_df = final_df.sort_values(by='Total_Variants', ascending=False)
final_df.to_csv(output_csv, index=False)

print(f"\nSUCCESS! Identified {len(final_df)} distinct structural hotspots.")
pd.set_option('display.max_colwidth', 100)
print(final_df.to_string(index=False))