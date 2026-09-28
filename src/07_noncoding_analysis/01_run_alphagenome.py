# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
import os
use_working_directory('outputs/noncoding')

import pandas as pd
import time
from alphagenome.data import genome
from alphagenome.models import dna_client

print("Initializing AlphaGenome client... ", flush=True)
API_KEY = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
if not API_KEY:
    raise RuntimeError("Set ALPHAGENOME_API_KEY in the environment before running this script.")
model = dna_client.create(API_KEY)

print("Loading vep_non_coding_filtered.csv...", flush=True)
df = pd.read_csv('vep_non_coding_filtered.csv')
unique_vars = df['Uploaded_variation'].unique()
ATAC_SCREEN_THRESHOLD = 0.5

results_list = []

for var in unique_vars:
    try:
        chrom, pos_str, alleles = var.split('_')
        pos = int(pos_str)
        ref, alt = alleles.split('/')
        
        if '-' in ref or '-' in alt:
            continue

        print(f"Scoring {var}... ", end="", flush=True)
        
        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt
        )
        
        # Use the absolute minimum supported context window
        half_window = 16384 // 2
        interval = genome.Interval(
            chromosome=chrom, 
            start=pos - half_window, 
            end=pos + half_window
        )
        
        outputs = model.predict_variant(
            interval=interval,
            variant=variant,
            requested_outputs=[dna_client.OutputType.RNA_SEQ, dna_client.OutputType.ATAC],
            ontology_terms=None 
        )
        
        ref_atac = outputs.reference.atac.values
        alt_atac = outputs.alternate.atac.values
        
        if ref_atac.size > 0 and alt_atac.size > 0:
            max_atac_delta = float(abs(alt_atac - ref_atac).max())
            print(f"Success! Max Global Delta: {max_atac_delta:.4f}", flush=True)
        else:
            max_atac_delta = None
            print("Warning: No ATAC data returned.", flush=True)
        
        results_list.append({
            'Uploaded_variation': var,
            'Status': 'Success',
            'Max_ATAC_Disruption_Global': max_atac_delta,
            'ATAC_Screen_Pass': (
                max_atac_delta >= ATAC_SCREEN_THRESHOLD
                if max_atac_delta is not None else None
            )
        })
        
    except Exception as e:
        print(f"Failed: {e}", flush=True)
        results_list.append({
            'Uploaded_variation': var,
            'Status': f'Failed: {e}',
            'Max_ATAC_Disruption_Global': None,
            'ATAC_Screen_Pass': None
        })

results_df = pd.DataFrame(results_list)
results_df.to_csv('alphagenome_scores.csv', index=False)
print("Done! Scores saved to alphagenome_scores.csv", flush=True)
