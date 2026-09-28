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

import urllib.request
import json
import time

variants = [
    {"id": "chr7_55259524_T/A", "region": "7:55259524-55259524"},
    {"id": "chr7_55240707_G/A", "region": "7:55240707-55240707"},
    {"id": "chr7_55238084_C/T", "region": "7:55238084-55238084"},
    {"id": "chr7_55238076_A/T", "region": "7:55238076-55238076"},
    {"id": "chr7_55259509_T/G", "region": "7:55259509-55259509"}
]

print("Querying Ensembl GRCh37 Regulatory Build for TF Motifs...\n")

for var in variants:
    print(f"Checking {var['id']}...")
    
    # We ask the API specifically for overlapping 'regulatory' and 'motif' features
    url = f"http://grch37.rest.ensembl.org/overlap/region/human/{var['region']}?feature=regulatory;feature=motif"
    
    try:
        req = urllib.request.Request(url, headers={'Content-Type': 'application/json'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            
            if not data:
                print("  -> No Ensembl regulatory or motif features found (Pure Dark Matter!).")
            else:
                for feature in data:
                    # Check for general enhancers/promoters
                    if feature['feature_type'] == 'regulatory':
                        desc = feature.get('description', 'Regulatory Region')
                        print(f"  -> Regulatory Element: {desc}")
                    
                    # Check for exact Transcription Factor binding motifs
                    elif feature['feature_type'] == 'motif':
                        matrix = feature.get('binding_matrix', 'Unknown TF')
                        score = feature.get('score', 'N/A')
                        print(f"  -> TF Motif Overlap: {matrix} (Binding Match Score: {score:.3f})")
                        
    except Exception as e:
        print(f"  -> API Error: {e}")
        
    time.sleep(1) # Be polite to the server
    print("-" * 40)