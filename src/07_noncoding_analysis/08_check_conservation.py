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
import urllib.error
import json
import time

# Our exact Neo-Enhancer mutations formatted for MyVariant (hg19)
variants = {
    "chr7:g.55259509T>G": "Locus 55259509 (Repressor Loss)",
    "chr7:g.55259524T>A": "Locus 55259524 (OLIG2 E-box Gain)"
}

print("Pivoting to MyVariant.info API for Conservation Scores...\n")

for hgvs, name in variants.items():
    # Ask the API specifically for the PhyloP and PhastCons scores
    url = f"https://myvariant.info/v1/variant/{hgvs}?fields=cadd.phylop,cadd.phast_cons"
    
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Python/3.10'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            
            print("========================================")
            print(f"📊 CONSERVATION RESULT: {name}")
            print("========================================")
            
            if 'cadd' in data:
                # Safely extract PhyloP
                p_score = data['cadd'].get('phylop', 0)
                if isinstance(p_score, dict):
                    p_score = list(p_score.values())[0] if p_score else 0
                
                # Safely extract PhastCons
                phast = data['cadd'].get('phast_cons', 0)
                if isinstance(phast, dict):
                    phast = list(phast.values())[0] if phast else 0

                print(f"PhyloP Score:     {p_score}")
                print(f"PhastCons Score:  {phast}\n")
                
                if float(p_score) > 1.5:
                    print("🔥 VERDICT: HIGHLY CONSERVED.")
                    print("Evolution strictly protected this exact base pair across 100 species.")
                    print("The Glioblastoma systematically attacked a critical biological anchor.")
                elif float(p_score) > 0:
                    print("VERDICT: Moderately Conserved.")
                else:
                    print("VERDICT: Low Conservation / Rapidly Evolving.")
            else:
                print("No conservation data found.")
            print("========================================\n")
            
    except urllib.error.URLError as e:
        print(f"Error querying {hgvs}: {e}")
        
    time.sleep(1) # Polite pause for the API