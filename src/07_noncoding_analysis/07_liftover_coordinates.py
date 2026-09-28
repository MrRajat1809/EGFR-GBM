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

# Our exact hg19/GRCh37 targets
targets = {
    "Window Start": 55000000,
    "EGFR Promoter (Approx)": 55080000,
    "Neo-Enhancer 1": 55259509,
    "Neo-Enhancer 2": 55259524,
    "Window End": 55350000
}

print("Performing Assembly Liftover (hg19 -> hg38) via Ensembl API...\n")

for name, pos in targets.items():
    # The API endpoint maps a specific point from GRCh37 to GRCh38
    url = f"https://rest.ensembl.org/map/human/GRCh37/7:{pos}..{pos}/GRCh38?content-type=application/json"
    
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            
            if data['mappings']:
                # Extract the new hg38 position
                hg38_pos = data['mappings'][0]['mapped']['start']
                print(f"{name}:")
                print(f"  hg19: chr7:{pos}  -->  hg38: chr7:{hg38_pos}\n")
            else:
                print(f"{name}: Could not map to hg38.\n")
                
    except Exception as e:
        print(f"Error mapping {name}: {e}")
        
    # Polite pause for the API
    time.sleep(0.5)