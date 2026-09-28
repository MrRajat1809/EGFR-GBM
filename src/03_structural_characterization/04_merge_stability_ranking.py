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

# --- PATH CONFIGURATION ---
GENOMIC_BASE = project_path('outputs/variant_prioritization')
STRUCTURAL_BASE = project_path('outputs/structural_analysis/thermodynamic_stability_ranking')
OUTPUT_FILE = os.path.join(GENOMIC_BASE, "EGFR_Final_Genomic_Structural_Integrated.csv")

def extract_pos(val):
    """Extracts numerical protein position for sorting."""
    if pd.isna(val) or val == '': return 0
    val = str(val).strip()
    m = re.search(r'(\d+)', val)
    return int(m.group(1)) if m else 0

def extract_mutation(val):
    """Normalizes variant strings to 1-letter format."""
    if pd.isna(val) or val == '': return None
    val = str(val).strip()
    
    # Handle MAESTRO format: C620.A{Y}
    if '{' in val:
        m = re.search(r'([A-Z])(\d+).+\{([A-Z])\}', val)
        if m: return f"{m.group(1)}{m.group(2)}{m.group(3)}"

    # Standard pattern (handles p.Ala289Val, Ala289Val, or A289V)
    m = re.search(r'([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})', val)
    if m:
        aa_map = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H',
                  'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
        wt, pos, mut = m.groups()
        wt = aa_map.get(wt.capitalize(), wt) if len(wt) > 1 else wt
        mut = aa_map.get(mut.capitalize(), mut) if len(mut) > 1 else mut
        return f"{wt.upper()}{pos}{mut.upper()}"
    return None

def main():
    print("🚀 Starting Protein-Position Sorted Integration...")

    # 1. LOAD & SORT MASTER GENOMIC TABLE
    df_master = pd.read_csv(os.path.join(GENOMIC_BASE, "EGFR_TCGA_Glioblastoma_Master_Consensus.csv"))
    df_master['prot_pos'] = df_master['HGVSp'].apply(extract_pos)
    df_master['mut_id'] = df_master['HGVSp'].apply(extract_mutation)
    # Primary sort by protein position
    df_master = df_master.sort_values(by='prot_pos').reset_index(drop=True)

    # 2. LOAD & SORT PREMPS
    print("Processing PremPS...")
    premps = pd.read_csv(os.path.join(STRUCTURAL_BASE, "PremPS_results.csv"))
    # In your file: 'Mutated chain' has the residue (e.g. E45K), 'Mutation' has the score
    premps['mut_id'] = premps['Mutated chain'].apply(extract_mutation)
    premps_subset = premps[['mut_id', 'Mutation']].rename(columns={'Mutation': 'PremPS_ddG'})

    # 3. LOAD & SORT MAESTRO
    print("Processing MAESTRO...")
    maestro = pd.read_csv(os.path.join(STRUCTURAL_BASE, "MAESTRO_results.csv"), sep=';')
    maestro['mut_id'] = maestro['substitution'].apply(extract_mutation)
    maestro_subset = maestro[['mut_id', 'ddG_pred']].rename(columns={'ddG_pred': 'MAESTRO_ddG'})

    # 4. LOAD & SORT mCSM (Normalizing -ve to +ve)
    print("Processing mCSM...")
    mcsm = pd.read_csv(os.path.join(STRUCTURAL_BASE, "mCSM_results.txt"), sep='\t')
    mcsm['mut_id'] = mcsm['WILD_RES'] + mcsm['RES_POS'].astype(str) + mcsm['MUT_RES']
    mcsm['mCSM_ddG_norm'] = mcsm['PRED_DDG'] * -1
    mcsm_subset = mcsm[['mut_id', 'mCSM_ddG_norm']]

    # --- MERGING ---
    print("Merging structural data into sorted Master Table...")
    final_df = df_master.merge(premps_subset, on='mut_id', how='left')
    final_df = final_df.merge(maestro_subset, on='mut_id', how='left')
    final_df = final_df.merge(mcsm_subset, on='mut_id', how='left')

    # Remove duplicates
    final_df = final_df.drop_duplicates(subset=['Genomic_Location', 'HGVSp'])
    
    # CALCULATE CONSENSUS SCORE
    ddg_cols = ['PremPS_ddG', 'MAESTRO_ddG', 'mCSM_ddG_norm']
    final_df['Consensus_Stability_ddG'] = final_df[ddg_cols].mean(axis=1)

    # FINAL SORT: Ascending order of protein position (No ranking by ddG)
    final_df = final_df.sort_values(by='prot_pos', ascending=True)

    # Cleanup temp columns
    final_df = final_df.drop(columns=['mut_id', 'prot_pos'])

    # --- SAVE ---
    try:
        final_df.to_csv(OUTPUT_FILE, index=False)
        print(f"✅ Success! Table sorted by Protein Position.")
        print(f"Output: {OUTPUT_FILE}")
    except PermissionError:
        print(f"❌ ERROR: Close '{OUTPUT_FILE}' in Excel and rerun.")

if __name__ == "__main__":
    main()