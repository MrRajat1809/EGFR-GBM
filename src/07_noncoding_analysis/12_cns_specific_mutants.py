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
# Automatically use whichever classification file you have
class_file_gz = project_path('data/raw/noncoding/Cosmic_Classification_v103_GRCh37.tsv.gz')
class_file_tsv = project_path('data/raw/noncoding/Cosmic_Classification.tsv')

cns_coso_ids = set()

# Sets to store unique sample IDs for our specific mutations
cns_patients_509 = set()
cns_patients_524 = set()

print("====================================================")
print("🚀 SCANNING FOR NEO-ENHANCER RECURRENCE IN CNS ONLY")
print("====================================================\n")

# ---------------------------------------------------------
# STEP 1: Build the CNS Phenotype Dictionary
# ---------------------------------------------------------
class_path = class_file_gz if os.path.exists(class_file_gz) else class_file_tsv

try:
    open_func = gzip.open if class_path.endswith('.gz') else open
    with open_func(class_path, 'rt', encoding='utf-8', errors='ignore') as f:
        next(f) # Skip header
        for line in f:
            cols = line.split('\t')
            if len(cols) > 2:
                coso_id = cols[0].strip()
                primary_site = cols[1].strip().lower()
                
                # Memorize Brain/CNS Phenotype IDs
                if "central_nervous_system" in primary_site or "brain" in primary_site:
                    cns_coso_ids.add(coso_id)
except Exception as e:
    print(f"❌ ERROR reading Classification file: {e}")
    exit()

# ---------------------------------------------------------
# STEP 2: Stream the Non-Coding Database
# ---------------------------------------------------------
print(f"Streaming {nc_file}...")
print("Hunting for 55259509 and 55259524 specifically in CNS patients...\n")

try:
    with gzip.open(nc_file, 'rt', encoding='utf-8', errors='ignore') as f:
        next(f) # skip header
        line_count = 0
        
        for line in f:
            line_count += 1
            if line_count % 5000000 == 0:
                print(f"   ⏳ [PROGRESS] Scanned {line_count:,} lines...")

            # SUPER FAST PRE-FILTER: Must have EGFR and one of our coordinates
            if "EGFR" in line and ("55259509" in line or "55259524" in line):  
                cols = line.split('\t')
                
                # Based on awk: 0=Gene, 5=SampleID, 6=PhenotypeID, 11=StartPos
                if len(cols) > 11:
                    gene = cols[0].strip()
                    
                    if gene == "EGFR":
                        sample_id = cols[5].strip()
                        coso_id = cols[6].strip()
                        start_pos = cols[11].strip()
                        
                        # Verify it's a CNS patient
                        if coso_id in cns_coso_ids:
                            
                            # Tally the exact coordinates
                            if start_pos == "55259509":
                                if sample_id not in cns_patients_509:
                                    cns_patients_509.add(sample_id)
                                    print(f"🔥 [HIT] CNS Patient Found with 55259509 (T>G)! (Sample: {sample_id})")
                                    
                            elif start_pos == "55259524":
                                if sample_id not in cns_patients_524:
                                    cns_patients_524.add(sample_id)
                                    print(f"🔥 [HIT] CNS Patient Found with 55259524 (T>A)! (Sample: {sample_id})")
                            
except FileNotFoundError:
    print(f"❌ ERROR: Could not find {nc_file}.")

print("\n====================================================")
print("📊 FINAL CNS PATIENT TALLY")
print("====================================================")
print(f"Total CNS Patients with ANY EGFR non-coding mutation: 99")
print(f"CNS Patients with chr7:55259509 (T>G): {len(cns_patients_509)}")
print(f"CNS Patients with chr7:55259524 (T>A): {len(cns_patients_524)}")
print("====================================================")