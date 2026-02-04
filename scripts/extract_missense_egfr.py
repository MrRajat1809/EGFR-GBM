import csv

input_file = "EGFR_GRCh37_vep_raw.tsv"
output_file = "EGFR_All_Missense_nsSNPs_GRCh37.csv"

with open(input_file) as f, open(output_file, "w", newline="") as out:
    writer = csv.writer(out)

    # Write CSV Header
    writer.writerow([
        "rsID", "Chr", "Pos_GRCh37", "Ref", "Alt", 
        "Gene", "Transcript_ID", "Protein_ID", 
        "Consequence", "HGVS_c", "HGVS_p", 
        "Protein_position", "AA_change", "Codon_change"
    ])

    for line in f:
        # Skip VEP headers and comments
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        
        # Adjust this number based on how you formatted your TSV
        if len(parts) < 6:
            continue

        # Logic for a custom TSV where CSQ is the last column
        rsid, chrom, pos, ref, alt, csq = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]

        # Split the VEP CSQ fields
        csq_fields = csq.split("|")

        # Safety: Ensure the CSQ string has all the fields you requested in VEP
        if len(csq_fields) < 11:
            continue

        consequence = csq_fields[5]
        if "missense_variant" not in consequence:
            continue

        # Extract fields (Mapping based on your VEP --fields flag)
        gene = csq_fields[3]
        transcript = csq_fields[4]
        hgvsc = csq_fields[6]
        hgvsp = csq_fields[7]
        protein_pos = csq_fields[8]
        aa_change = csq_fields[9]
        codon_change = csq_fields[10]

        # Only keep records that actually resulted in a protein change
        if not hgvsp or hgvsp == "":
            continue

        writer.writerow([
            rsid, chrom, pos, ref, alt, 
            gene, transcript, "P00533", # Hardcoded EGFR UniProt
            consequence, hgvsc, hgvsp, 
            protein_pos, aa_change, codon_change
        ])