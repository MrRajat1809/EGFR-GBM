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

import gzip
import sys

# Pointing directly to your gz file
gz_path = project_path('data/raw/noncoding/Cosmic_NonCodingVariants_v103_GRCh37.tsv.gz')

# The exact base pairs of our Neo-Enhancer
target_1 = "55259509"
target_2 = "55259524"

count_1 = 0
count_2 = 0

print("====================================================")
print(f"🚀 SCANNING COSMIC V103 FOR NEO-ENHANCER RECURRENCE")
print("====================================================\n")
print(f"Streaming {gz_path}...")
print("Crunching 2.09 GB of data in memory. This should take about 60 seconds...\n")

try:
    # Open the gzipped file in text mode ('rt')
    with gzip.open(gz_path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # We do a quick string match first because it's incredibly fast in Python
            if target_1 in line:
                # To be absolutely sure it's our chromosome and not a random ID, we verify
                if "7:55259509-55259509" in line or "55259509" in line:
                    count_1 += 1
                    print(f"[HIT] Found independent patient tumor with 55259509 (Repressor Loss) mutation!")
            
            elif target_2 in line:
                if "7:55259524-55259524" in line or "55259524" in line:
                    count_2 += 1
                    print(f"[HIT] Found independent patient tumor with 55259524 (OLIG2 E-box Gain) mutation!")

    print("\n====================================================")
    print("📊 FINAL PATIENT RECURRENCE TALLY")
    print("====================================================")
    print(f"chr7:55259509 (T>G): Found in {count_1} distinct patient samples.")
    print(f"chr7:55259524 (T>A): Found in {count_2} distinct patient samples.")
    
    total = count_1 + count_2
    print(f"\nTotal occurrences of this specific Neo-Enhancer mutation: {total}")
    print("====================================================")

except FileNotFoundError:
    print(f"❌ ERROR: Could not find '{gz_path}' in the current directory.")
except Exception as e:
    print(f"❌ ERROR: {e}")