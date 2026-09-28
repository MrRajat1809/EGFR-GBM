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

# Our top 5 targets
variants = [
    {"id": "chr7_55259524_T_A", "chrom": "7", "pos": 55259524, "ref": "T", "alt": "A"},
    {"id": "chr7_55240707_G_A", "chrom": "7", "pos": 55240707, "ref": "G", "alt": "A"},
    {"id": "chr7_55238084_C_T", "chrom": "7", "pos": 55238084, "ref": "C", "alt": "T"},
    {"id": "chr7_55238076_A_T", "chrom": "7", "pos": 55238076, "ref": "A", "alt": "T"},
    {"id": "chr7_55259509_T_G", "chrom": "7", "pos": 55259509, "ref": "T", "alt": "G"}
]

# We want +/- 15 base pairs around the mutation
WINDOW = 15 

fasta_output = []

print("Extracting genomic sequences from Ensembl (GRCh37)...\n")

for var in variants:
    start = var['pos'] - WINDOW
    end = var['pos'] + WINDOW
    
    url = f"http://grch37.rest.ensembl.org/sequence/region/human/{var['chrom']}:{start}..{end}:1?content-type=application/json"
    
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            wt_seq = data['seq']
            
            # Verify the reference base matches Ensembl to ensure we are on the right strand
            actual_ref = wt_seq[WINDOW]
            if actual_ref != var['ref']:
                print(f"Warning for {var['id']}: Expected REF {var['ref']}, but Ensembl shows {actual_ref}. Strand issue?")
            
            # Create the mutated sequence
            mut_seq = wt_seq[:WINDOW] + var['alt'] + wt_seq[WINDOW+1:]
            
            # Format as FASTA
            fasta_output.append(f">{var['id']}_WildType")
            fasta_output.append(wt_seq)
            fasta_output.append(f">{var['id']}_Mutated")
            fasta_output.append(mut_seq)
            
            print(f"Extracted {var['id']}")
            
    except Exception as e:
        print(f"Error fetching {var['id']}: {e}")
        
    time.sleep(1) # Be polite to Ensembl

# Save to file
with open("jaspar_motif_scan_input.fasta", "w") as f:
    f.write("\n".join(fasta_output) + "\n")

print("\nDone! Sequences saved to 'jaspar_motif_scan_input.fasta'")