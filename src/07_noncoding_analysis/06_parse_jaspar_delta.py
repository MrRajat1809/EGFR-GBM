# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('outputs/noncoding')

import pandas as pd

print("Loading JASPAR results from copy/paste...")
try:
    # Read the pasted TSV file
    df = pd.read_csv(project_path('data/raw/noncoding/jaspar_results.tsv'), sep='\t')
    
    # Strip whitespace and make columns lowercase
    df.columns = df.columns.str.strip().str.lower()
    
    # Auto-detect the right columns
    seq_col = [c for c in df.columns if 'sequence' in c and 'id' in c][0] if any('sequence' in c and 'id' in c for c in df.columns) else 'sequence id'
    name_col = [c for c in df.columns if 'name' in c][0] if any('name' in c for c in df.columns) else 'name'
    score_col = [c for c in df.columns if 'relative' in c][0] if any('relative' in c for c in df.columns) else 'relative score'
    
    # Clean up the sequence names to pair WT and Mutated
    df['variant'] = df[seq_col].apply(lambda x: str(x).replace('_WildType', '').replace('_Mutated', '').strip())
    df['state'] = df[seq_col].apply(lambda x: 'WT' if 'WildType' in str(x) else 'Mutated')
    
    df[score_col] = pd.to_numeric(df[score_col], errors='coerce')
    df = df.dropna(subset=[score_col])
    
    # THE FIX: If the max score is around 1.0, convert the decimal to a percentage
    if df[score_col].max() <= 2.0:
        df[score_col] = df[score_col] * 100.0
    
    # Pivot the data
    pivot_df = df.pivot_table(index=['variant', name_col], 
                              columns='state', 
                              values=score_col, 
                              aggfunc='max').reset_index()
    
    pivot_df = pivot_df.fillna(0)
    pivot_df['delta'] = pivot_df['Mutated'] - pivot_df['WT']
    
    # Filter for massive swings (>5% shift or complete appearance/disappearance)
    significant = pivot_df[(abs(pivot_df['delta']) >= 5.0) | 
                           ((pivot_df['WT'] == 0) & (pivot_df['Mutated'] > 85.0)) |
                           ((pivot_df['Mutated'] == 0) & (pivot_df['WT'] > 85.0))]
    
    significant = significant.sort_values(by='delta', key=abs, ascending=False)
    
    print("\n" + "="*60)
    print("🔥 TOP MOTIF DISRUPTIONS (GAIN OR LOSS OF BINDING) 🔥")
    print("="*60)
    
    count = 0
    for _, row in significant.head(20).iterrows():
        var = row['variant'].replace('>', '') 
        
        # Clean up the TF name (e.g., MA2542.1.SCAND3 -> SCAND3)
        raw_tf = str(row[name_col])
        tf = raw_tf.split('.')[-1] if '.' in raw_tf else raw_tf
        
        wt_score = row['WT']
        mut_score = row['Mutated']
        delta = row['delta']
        
        if delta > 0:
            status = "🟢 GAINED BINDING"
        else:
            status = "🔴 LOST BINDING  "
            
        print(f"Variant: {var}")
        print(f"TF Name: {tf}")
        print(f"Change:  {status} ({delta:+.2f}%)  [WT: {wt_score:.2f}% -> Mut: {mut_score:.2f}%]")
        print("-" * 40)
        count += 1

    if count == 0:
        print("No massive swings >5% found. The mutations might be subtly shifting affinity rather than destroying motifs.")

    significant.to_csv('jaspar_significant_deltas.csv', index=False)
    print(f"\nSaved {len(significant)} significant delta events to 'jaspar_significant_deltas.csv'")

except Exception as e:
    print(f"Error parsing file: {e}")
