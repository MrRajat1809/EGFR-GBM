# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('outputs/noncoding/CADD')

import gzip
import csv

input_file = project_path('data/raw/noncoding/CADD/GRCh37-v1.7_anno_21a52e3ce66a0dba879293aea9ea2591.tsv.gz')
output_file = "CADD_Ranked_Variants.csv"

print(f"🚀 Ripping open {input_file}...")

# Helper function to bypass the browser auto-unzip trap
def open_file_safely(filepath):
    # Read the first 2 bytes to check the 'magic number' of the file
    with open(filepath, 'rb') as test_f:
        is_gzipped = test_f.read(2) == b'\x1f\x8b'
    
    if is_gzipped:
        return gzip.open(filepath, 'rt')
    else:
        print("   (Detected uncompressed text despite the .gz extension. Adjusting...)")
        return open(filepath, 'rt', encoding='utf-8')

extracted_data = []

try:
    with open_file_safely(input_file) as f:
        for line in f:
            # Skip the CADD meta-information comments and header row
            if line.startswith("##") or line.startswith("#Chrom"):
                continue

            # Process the actual data rows
            cols = line.strip().split('\t')
            
            # Ensure the row actually has data (CADD output has over 100 columns)
            if len(cols) > 20:
                chrom = cols[0]
                pos = cols[1]
                ref = cols[2]
                alt = cols[3]
                
                # In CADD, column 7 is AnnoType and column 8 is Consequence
                anno_type = cols[6]
                consequence = cols[7]
                
                # PHRED score is ALWAYS the very last column
                try:
                    phred = float(cols[-1])
                except ValueError:
                    phred = 0.0

                extracted_data.append({
                    "Coordinate": f"chr{chrom}:{pos}",
                    "Mutation": f"{ref}>{alt}",
                    "Region_Type": anno_type,
                    "Biological_Impact": consequence,
                    "PHRED_Score": phred
                })

    # Sort the list of dictionaries by PHRED score (Highest to Lowest)
    extracted_data.sort(key=lambda x: x["PHRED_Score"], reverse=True)

    # Write everything to a clean CSV
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ["Coordinate", "Mutation", "Region_Type", "Biological_Impact", "PHRED_Score"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in extracted_data:
            writer.writerow(row)

    print(f"✅ BOOM. Extracted {len(extracted_data)} annotations.")
    print(f"✅ Data sorted and saved cleanly to: {output_file}")

except FileNotFoundError:
    print(f"❌ ERROR: Could not find {input_file}. Make sure the filename is exact!")