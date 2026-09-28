# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('.')

import pandas as pd
import os

INPUT_CSV = project_path('data/processed/variants/EGFR_TCGA_GBM_Missense_Starting_List.csv')
OUTPUT_VCF = project_path('data/processed/variants/EGFR_TCGA_51_Somatic.vcf')

def main():
    print(f"Loading TCGA Missense List...")
    df = pd.read_csv(INPUT_CSV)
    
    # 1. Isolate the unique genomic coordinates (we only need to run each unique mutation once)
    # TCGA MAF columns: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2
    unique_muts = df[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']].drop_duplicates()
    
    print(f"Found {len(unique_muts)} unique somatic mutations to format for VEP.")
    
    # 2. Open the VCF file and write the standard header
    with open(OUTPUT_VCF, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=TCGA_GBM_cBioPortal\n")
        vcf.write("##reference=GRCh37\n") # TCGA PanCancer uses GRCh37/hg19 usually. Change to GRCh38 if you use hg38!
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # 3. Write each mutation as a standard VCF row
        for _, row in unique_muts.iterrows():
            chrom = str(row['Chromosome']).replace('chr', '') # Ensure it's just '7', not 'chr7'
            pos = row['Start_Position']
            ref = row['Reference_Allele']
            alt = row['Tumor_Seq_Allele2']
            
            # ID is '.' because we don't know the rsID yet! VEP will fill this in.
            vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n"
            vcf.write(vcf_line)
            
    print(f"\nSUCCESS! Created {OUTPUT_VCF}")
    print("You can now feed this VCF directly into your Ensembl VEP container!")

if __name__ == "__main__":
    main()