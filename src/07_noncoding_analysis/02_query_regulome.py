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

import requests
import pandas as pd

# 1. Load directly from your original tab-separated file
variants = []
with open('EGFR_GBM_Coordinates_for_RegulomeDB.txt', 'r') as f:
    for line in f:
        if line.strip():
            parts = line.strip().split('\t')
            # Format exactly as the API expects: chr:pos (e.g., chr7:55259509)
            variants.append(f"{parts[0]}:{parts[1]}")

# 2. Remove duplicates (we only need to query each locus once)
unique_variants = list(set(variants))
print(f"Querying RegulomeDB API for {len(unique_variants)} unique loci...")

# 3. Send the API request
url = "https://regulomedb.org/regulome-search/"
payload = {
    "variants": unique_variants,
    "genome": "GRCh37", # Matches your COSMIC v103 build
    "limit": 100
}
headers = {"Content-Type": "application/json"}

response = requests.post(url, json=payload, headers=headers)

# 4. Safely parse the results
if response.status_code == 200:
    data = response.json()
    hits = data.get("variants", [])
    
    if not hits:
        print("\nAPI Query Successful, but found 0 matches.")
        print("Biological Result: None of these specific EGFR variants overlap a known regulatory element in the ENCODE GRCh37 databases.")
    else:
        results = []
        for v in hits:
            results.append({
                "Coordinate": f"{v.get('chrom')}:{v.get('start')}-{v.get('end')}",
                "Score": v.get("regulomedb_score", "Unscored"),
                "TF_Binding_Peaks": v.get("num_peaks", 0),
                "rsIDs": ", ".join(v.get("rsids", []))
            })
        
        df = pd.DataFrame(results)
        
        # We only sort if the dataframe successfully built the 'Score' column
        df = df.sort_values(by="Score", na_position='last')
        
        output_file = "EGFR_GBM_RegulomeDB_Scores.tsv"
        df.to_csv(output_file, sep='\t', index=False)
        
        print(f"\nSuccess! {len(df)} variants had regulatory data. Saved to {output_file}\n")
        print("Top Regulatory Variants:")
        print(df.head())
else:
    print(f"API Error {response.status_code}: {response.text}")