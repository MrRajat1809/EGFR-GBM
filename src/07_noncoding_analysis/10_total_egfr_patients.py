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
import os

nc_file = project_path('data/raw/noncoding/Cosmic_NonCodingVariants_v103_GRCh37.tsv.gz')
class_file_gz = project_path('data/raw/noncoding/Cosmic_Classification_v103_GRCh37.tsv.gz')
class_file_tsv = project_path('data/raw/noncoding/Cosmic_Classification.tsv')

cns_coso_ids = set()
unique_patients_all = set()
unique_patients_cns = set()

print("====================================================")
print("🚀 COSMIC V103 RELATIONAL DATABASE SCANNER")
print("====================================================\n")

# ---------------------------------------------------------
# STEP 1: Build the CNS Phenotype Dictionary
# ---------------------------------------------------------
print("[1/2] Mapping CNS Phenotype IDs from Classification table...")

class_path = class_file_gz if os.path.exists(class_file_gz) else class_file_tsv

try:
    # Use gzip.open if it's zipped, standard open if unzipped
    open_func = gzip.open if class_path.endswith('.gz') else open
    
    with open_func(class_path, 'rt', encoding='utf-8', errors='ignore') as f:
        next(f) # Skip header
        for line in f:
            cols = line.split('\t')
            if len(cols) > 2:
                coso_id = cols[0].strip()
                primary_site = cols[1].strip().lower()
                
                # If it's a brain/CNS tissue, memorize the COSO ID
                if "central_nervous_system" in primary_site or "brain" in primary_site:
                    cns_coso_ids.add(coso_id)
                    
    print(f"   ✅ Memorized {len(cns_coso_ids)} unique CNS Phenotype IDs.\n")
except Exception as e:
    print(f"❌ ERROR: Could not read Classification file: {e}")
    print("Please make sure the Cosmic_Classification file is in the directory!")
    exit()

# ---------------------------------------------------------
# STEP 2: Stream the Non-Coding Database
# ---------------------------------------------------------
print(f"[2/2] Streaming 2.09 GB Non-Coding Database ({nc_file})...")
print("Hunting for EGFR mutations mapped to CNS IDs...\n")

try:
    with gzip.open(nc_file, 'rt', encoding='utf-8', errors='ignore') as f:
        next(f) # skip header
        line_count = 0
        
        for line in f:
            line_count += 1
            if line_count % 5000000 == 0:
                print(f"   ⏳ [PROGRESS] Scanned {line_count:,} lines...")

            if "EGFR" in line:  
                cols = line.split('\t')
                
                # Based on your awk command: 1=Gene, 6=SampleID, 7=PhenotypeID
                if len(cols) > 6:
                    gene = cols[0].strip()
                    
                    if gene == "EGFR":
                        sample_id = cols[5].strip()
                        coso_id = cols[6].strip()
                        
                        unique_patients_all.add(sample_id)
                        
                        if coso_id in cns_coso_ids:
                            if sample_id not in unique_patients_cns:
                                unique_patients_cns.add(sample_id)
                                print(f"✨ [HIT] Found CNS patient! (Total unique CNS so far: {len(unique_patients_cns)})")
                            
except FileNotFoundError:
    print(f"❌ ERROR: Could not find {nc_file}.")

print("\n====================================================")
print("📊 FINAL EGFR PATIENT COUNTS")
print("====================================================")
print(f"Total unique patients with ANY non-coding EGFR mutation (All Cancers): {len(unique_patients_all)}")
print(f"Total unique patients with ANY non-coding EGFR mutation (CNS/Brain):   {len(unique_patients_cns)}")
print("====================================================")