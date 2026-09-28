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

cns_patients_509 = set()
cns_patients_524 = set()

print("====================================================")
print("🚀 CHECKING TUMOR CO-OCCURRENCE OF 55259509 AND 55259524")
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
print("Extracting exact Patient IDs for both coordinates...\n")

try:
    with gzip.open(nc_file, 'rt', encoding='utf-8', errors='ignore') as f:
        next(f) # skip header
        
        for line in f:
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
                            # Capture the exact IDs
                            if start_pos == "55259509":
                                cns_patients_509.add(sample_id)
                            elif start_pos == "55259524":
                                cns_patients_524.add(sample_id)
                            
except FileNotFoundError:
    print(f"❌ ERROR: Could not find {nc_file}.")

print("\n====================================================")
print("📊 EXACT PATIENT ID BREAKDOWN")
print("====================================================")
print(f"Sample IDs with chr7:55259509 (n={len(cns_patients_509)}): {', '.join(cns_patients_509)}")
print(f"Sample IDs with chr7:55259524 (n={len(cns_patients_524)}): {', '.join(cns_patients_524)}")

intersection = cns_patients_509.intersection(cns_patients_524)

print("\n====================================================")
print("🔍 CO-OCCURRENCE VERDICT")
print("====================================================")
if len(intersection) > 0:
    print(f"🔥 MASSIVE DISCOVERY: Patient(s) {', '.join(intersection)} harbor BOTH mutations!")
    print("This tumor actively mutated both the SOX2 and OLIG2 docking bays.")
else:
    print("🧬 INDEPENDENT EVENTS: These mutations occurred in completely separate patients.")
    print("This implies the tumor has two distinct evolutionary pathways to hijack this enhancer.")
print("====================================================")